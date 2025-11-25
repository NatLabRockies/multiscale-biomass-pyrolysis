#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Prob.H>
#include <TranspReact.H>
#include <Tagging.H>
#include <Chemistry.H>
#include <ProbParm.H>
#include <VarDefines.H>

using namespace amrex;

ProbParm* TranspReact::h_prob_parm = nullptr;
ProbParm* TranspReact::d_prob_parm = nullptr;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
TranspReact::TranspReact()
{
    ReadParameters();
    h_prob_parm = new ProbParm{};
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));
    amrex_probinit(*h_prob_parm, *d_prob_parm);

    allvarnames.resize(NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; i++)
    {
        allvarnames[i] = tr_species::specnames[i];
    }

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev)
    {
        nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    for(int lev=0;lev<nlevs_max;lev++)
    {
        dt[lev]=fixed_dt;
    }

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
    
    ParmParse pp("tr");
    amrex::Vector<int> def_bc_lo{ROBINBC, ROBINBC, ROBINBC};
    amrex::Vector<int> def_bc_hi{ROBINBC, ROBINBC, ROBINBC};
    
    pp.getarr("default_bc_lo", def_bc_lo, 0, AMREX_SPACEDIM);
    pp.getarr("default_bc_hi", def_bc_hi, 0, AMREX_SPACEDIM);
    
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        //make it size 3 all the time
        all_bcs_lo[sp].resize(3, ROBINBC);
        all_bcs_hi[sp].resize(3, ROBINBC);
        
        for(int n=0;n<AMREX_SPACEDIM;n++)
        {
           all_bcs_lo[sp][n]=def_bc_lo[n];
           all_bcs_hi[sp][n]=def_bc_hi[n];
        }
        
        pp.queryarr((allvarnames[sp]+"_bc_lo").c_str(), all_bcs_lo[sp], 0, AMREX_SPACEDIM);
        pp.queryarr((allvarnames[sp]+"_bc_hi").c_str(), all_bcs_hi[sp], 0, AMREX_SPACEDIM);
    }

    //foextrap all states as bcs imposed
    //through linear solver
    bcspec.resize(NUM_SPECIES);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {

        int bctype=(geom[0].isPeriodic(idim))?BCType::int_dir:BCType::foextrap;

        for (int sp=0; sp < NUM_SPECIES; sp++) 
        {
            bcspec[sp].setLo(idim, bctype);
            bcspec[sp].setHi(idim, bctype);
        }
    }
}

TranspReact::~TranspReact()
{
    delete h_prob_parm;
    The_Arena()->free(d_prob_parm);
}
// initializes multilevel data
void TranspReact::InitData()
{
    ProbParm* localprobparm = d_prob_parm;

    if (restart_chkfile == "")
    {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        
        if(transform_vars)
        {
            Vector<MultiFab> Sborder(finest_level+1);
            for(int lev=0;lev<=finest_level;lev++)
            {
                Sborder[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), ngrow_for_fillpatch);
                Sborder[lev].setVal(0.0);

                FillPatch(lev, time, Sborder[lev], 0, Sborder[lev].nComp());
            }
            transform_variables(Sborder,time);
        }

        AverageDown();

        if (chk_int > 0 || chk_time > 0.0)
        {
            WriteCheckpointFile(0);
        }

    } 
    else
    {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0 || plot_time > 0.0)
    {
        WritePlotFile(amrex::Math::floor
                      (amrex::Real(istep[0])/amrex::Real(plot_int)));
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void TranspReact::ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        ParmParse pp("tr");
        if (pp.contains("tagged_vars"))
        {
            int nvars = pp.countval("tagged_vars");
            refine_phi.resize(nvars);
            refine_phigrad.resize(nvars);
            refine_phi_comps.resize(nvars);
            std::string varname;
            for (int i = 0; i < nvars; i++)
            {
                pp.get("tagged_vars", varname, i);
                pp.get((varname + "_refine").c_str(), refine_phi[i]);
                pp.get((varname + "_refinegrad").c_str(), refine_phigrad[i]);

                int spec_id = tr_species::find_id(varname);
                if (spec_id == -1)
                {
                    Print() << "Variable name:" << varname << " not found for tagging\n";
                    amrex::Abort("Invalid tagging variable");
                }
                else
                {
                    refine_phi_comps[i] = spec_id;
                }
            }
        }
    }

    if (refine_phi.size() == 0) return;

    const int tagval = TagBox::SET;

    const MultiFab& state = phi_new[lev];
    MultiFab Sborder(grids[lev], dmap[lev], state.nComp(), 1);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto statefab = Sborder.array(mfi);
            const auto tagfab = tags.array(mfi);

            amrex::Real* refine_phi_dat = refine_phi.data();
            amrex::Real* refine_phigrad_dat = refine_phigrad.data();
            int* refine_phi_comps_dat = refine_phi_comps.data();
            int ntagged_comps = refine_phi_comps.size();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                state_based_refinement(i, j, k, tagfab, statefab, refine_phi_dat, refine_phi_comps_dat, ntagged_comps, tagval);
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                stategrad_based_refinement(i, j, k, tagfab, statefab, refine_phigrad_dat, refine_phi_comps_dat, ntagged_comps, tagval);
            });
        }
    }
}

// read in some parameters from inputs file
void TranspReact::ReadParameters()
{
    {
        ParmParse pp; // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.
        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        plot_int_old=plot_int;
        pp.query("plot_int_old", plot_int_old);
        pp.query("plot_time", plot_time);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        chk_int_old=chk_int;
        pp.query("chk_int_old", chk_int_old);
        pp.query("chk_time", chk_time);
        pp.query("restart", restart_chkfile);
    }

    {
        ParmParse pp("tr");
        pp.query("dt", fixed_dt);

        pp.query("linsolve_reltol",linsolve_reltol);
        pp.query("linsolve_abstol",linsolve_abstol);
        pp.query("linsolve_bot_reltol",linsolve_bot_reltol);
        pp.query("linsolve_bot_abstol",linsolve_bot_abstol);

        pp.query("linsolve_num_pre_smooth",linsolve_num_pre_smooth);
        pp.query("linsolve_num_post_smooth",linsolve_num_post_smooth);
        pp.query("linsolve_num_final_smooth",linsolve_num_final_smooth);
        pp.query("linsolve_num_bottom_smooth",linsolve_num_bottom_smooth);

        pp.query("linsolve_maxiter",linsolve_maxiter);
        pp.query("linsolve_max_coarsening_level",linsolve_max_coarsening_level);
        pp.query("min_species_conc",min_species_conc);
        pp.query("hyp_order",hyp_order);
        pp.query("do_reactions",do_reactions);
        pp.query("do_transport",do_transport);
        pp.query("do_advection",do_advection);
        pp.query("transform_vars",transform_vars);

        Vector<int> steady_specid_list;
        Vector<int> unsolved_specid_list;
        pp.queryarr("steady_species_ids", steady_specid_list);
        pp.queryarr("unsolved_species_ids", unsolved_specid_list);

        for(unsigned int i=0;i<steady_specid_list.size();i++)
        {
            steadyspec[steady_specid_list[i]]=1;    
        }

        for(unsigned int i=0;i<unsolved_specid_list.size();i++)
        {
            unsolvedspec[unsolved_specid_list[i]]=1;    
        }


        if(hyp_order==1) //first order upwind
        {
            ngrow_for_fillpatch=1;
        }
        else
        {
            ngrow_for_fillpatch=3;
        }

        //in case we need to set it manually
        pp.query("ngrow",ngrow_for_fillpatch);
        pp.query("num_timestep_correctors",num_timestep_correctors);
        pp.query("num_split_correctors",num_split_correctors);
        pp.query("adaptive_dt",adaptive_dt);
        pp.query("advective_cfl",advective_cfl);
        pp.query("diffusive_cfl",diffusive_cfl);
        pp.query("split_chemistry",split_chemistry);
        pp.query("dt_min",dt_min);
        pp.query("dt_max",dt_max);
        
        pp.query("do_reactions",do_reactions);

#ifdef AMREX_USE_HYPRE
        pp.query("use_hypre",use_hypre);
#endif
        //local cleanup
        steady_specid_list.clear();
        unsolved_specid_list.clear();

        if(split_chemistry)
        {
           ParmParse pp_int("integration");
           pp_int.query("type",integration_type);
           ParmParse pp_int_sd("integration.sundials");
           pp_int_sd.query("type",integration_sd_type);
        }

    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void TranspReact::GetData(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_new[lev]);
    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(&phi_old[lev]);
        datatime.push_back(t_old[lev]);
    } else
    {
        data.push_back(&phi_old[lev]);
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}
