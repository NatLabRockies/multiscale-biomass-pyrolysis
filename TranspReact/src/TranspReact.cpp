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
#include <UserBCs.H>

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

    allvarnames.resize(NVAR);
    for (int i = 0; i < NUM_SPECIES; i++)
    {
        allvarnames[i] = tr_species::specnames[i];
    }
    allvarnames[CMASK_ID]="cellmask";

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
    bcspec.resize(NVAR);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {

        int bctype=(geom[0].isPeriodic(idim))?BCType::int_dir:BCType::foextrap;

        for (int sp=0; sp < NVAR; sp++) 
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
        smooth_cellmask();
        
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
        pp.query("linsolve_use_prvs_soln",linsolve_use_prvs_soln);

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
        pp.query("interface_update_maxiter",interface_update_maxiter);


        Vector<int> steady_specid_list;
        pp.queryarr("steady_species_ids", steady_specid_list);
        for(unsigned int i=0;i<steady_specid_list.size();i++)
        {
            steadyspec[steady_specid_list[i]]=1;    
        }
        
        Vector<int> unsolved_specid_list;
        pp.queryarr("unsolved_species_ids", unsolved_specid_list);
        for(unsigned int i=0;i<unsolved_specid_list.size();i++)
        {
            unsolvedspec[unsolved_specid_list[i]]=1;    
        }
        
        Vector<int> conjsolve_specid_list;
        pp.queryarr("conjugate_solve_species_ids", conjsolve_specid_list);
        pp.query("conjsolve_maxiter",conjsolve_maxiter);
        for(unsigned int i=0;i<conjsolve_specid_list.size();i++)
        {
            conjugate_solve[conjsolve_specid_list[i]]=1;    
        } 
        
        Vector<int> under_relax_specid_list;
        Vector<amrex::Real> under_relax_fac_list;
        Vector<int> under_relax_maxiter_list;
        Vector<Real> under_relax_tol_list;
        int enable_under_relaxation=0;
        pp.queryarr("under_relax_species_ids", under_relax_specid_list);
        pp.queryarr("under_relax_fac",under_relax_fac_list);
        pp.queryarr("under_relax_maxiter",under_relax_maxiter_list);
        pp.queryarr("under_relax_tol",under_relax_tol_list);

        if(under_relax_specid_list.size()>0)
        {
           enable_under_relaxation=1;
        }

        //set default values
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
          under_relax[sp]=0;
          relaxfac[sp]=0.0;
          under_relax_maxiter[sp]=1;
          under_relax_tol[sp]=1e-5;
        }

        if(enable_under_relaxation)
        {
            if(under_relax_fac_list.size()!=under_relax_specid_list.size() ||
               under_relax_maxiter_list.size()!=under_relax_specid_list.size() ||
               under_relax_maxiter_list.size()!=under_relax_fac_list.size() ||
               under_relax_specid_list.size()!=under_relax_tol_list.size())
            {
                //there may be other false conditions..
                amrex::Print()<<"under_relax_species_ids, under_relax_fac, " 
                             <<"under_relax_maxiter, and under_relax_tol should "
                             <<"have same length \n";
                amrex::Abort("Under relaxation parameter lengths dont match up.");
            }

            for(unsigned int i=0;i<under_relax_specid_list.size();i++)
            {
                under_relax[under_relax_specid_list[i]]=1;   
                relaxfac[under_relax_specid_list[i]]=under_relax_fac_list[i];
                under_relax_maxiter[under_relax_specid_list[i]]=under_relax_maxiter_list[i];
                under_relax_tol[under_relax_specid_list[i]]=under_relax_tol_list[i];
            }
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

        pp.query("using_ib",using_ib);
        if(using_ib)
        {
            ngrow_for_fillpatch=3;
        }
        pp.query("cmask_smoothing_iters",nsmoothiters);

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

//IB functions
void TranspReact::null_bcoeff_at_ib(int ilev, Array<MultiFab, 
                                    AMREX_SPACEDIM>& face_bcoeff, 
                                    MultiFab& Sborder,int conjsolve)
{
    int captured_conjsolve=conjsolve;

    for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array<Box,AMREX_SPACEDIM> face_boxes;
        face_boxes[0] = mfi.nodaltilebox(0);
#if AMREX_SPACEDIM > 1
        face_boxes[1] = mfi.nodaltilebox(1);
#if AMREX_SPACEDIM == 3
        face_boxes[2] = mfi.nodaltilebox(2);
#endif
#endif

        Array4<Real> sb_arr = Sborder.array(mfi);
        GpuArray<Array4<Real>, AMREX_SPACEDIM> 
        face_bcoeff_arr{AMREX_D_DECL(face_bcoeff[0].array(mfi), 
                                     face_bcoeff[1].array(mfi), face_bcoeff[2].array(mfi))};

        for(int idim=0;idim<AMREX_SPACEDIM;idim++)
        {
            amrex::ParallelFor(face_boxes[idim], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                IntVect face{AMREX_D_DECL(i,j,k)};
                IntVect lcell{AMREX_D_DECL(i,j,k)};
                IntVect rcell{AMREX_D_DECL(i,j,k)};

                lcell[idim]-=1;

                int mask_L,mask_R;

                if(!captured_conjsolve)
                {
                    mask_L=int(sb_arr(lcell,CMASK_ID));
                    mask_R=int(sb_arr(rcell,CMASK_ID));
                }
                else
                {
                    mask_L=int(1.0-sb_arr(lcell,CMASK_ID));
                    mask_R=int(1.0-sb_arr(rcell,CMASK_ID));
                }

                //1 when both mask_L and mask_R are 0
                int covered_interface=(!mask_L)*(!mask_R);
                //1 when both mask_L and mask_R are 1
                int regular_interface=(mask_L)*(mask_R);
                //1-0 or 0-1 interface
                int cov_uncov_interface=(mask_L)*(!mask_R)+(!mask_L)*(mask_R);

                if(cov_uncov_interface) //1*0 0*1 cases
                {
                    face_bcoeff_arr[idim](face)=0.0;
                }
                else if(covered_interface) //0*0 case
                {
                    //keeping bcoeff non zero in dead cells just in case
                    face_bcoeff_arr[idim](face)=1.0;
                }
                else
                {
                    //do nothing
                }
            });
        }

    }
}

void TranspReact::set_explicit_fluxes_at_ib(int ilev, MultiFab& rhs,
                                            MultiFab& acoeff,
                                            MultiFab& Sborder,
                                            Real time,int compid,int conjsolve)
{
    Real captured_time=time;
    int solved_comp=compid;
    int captured_conjsolve=conjsolve;
    ProbParm const* localprobparm = d_prob_parm;

    for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto dx = geom[ilev].CellSizeArray();
        auto prob_lo = geom[ilev].ProbLoArray();
        auto prob_hi = geom[ilev].ProbHiArray();
        const Box& domain = geom[ilev].Domain();
        const int* domlo_arr = geom[ilev].Domain().loVect();
        const int* domhi_arr = geom[ilev].Domain().hiVect();

        GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
        GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

        Array4<Real> sb_arr = Sborder.array(mfi);
        Array4<Real> rhs_arr = rhs.array(mfi);
        Array4<Real> acoeff_arr = acoeff.array(mfi);

        Array<Box,AMREX_SPACEDIM> face_boxes;
        face_boxes[0] = mfi.nodaltilebox(0);
#if AMREX_SPACEDIM > 1
        face_boxes[1] = mfi.nodaltilebox(1);
#if AMREX_SPACEDIM == 3
        face_boxes[2] = mfi.nodaltilebox(2);
#endif
#endif
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) { 

            int cmask=(captured_conjsolve==0)?
            int(sb_arr(i,j,k,CMASK_ID)):int(1.0-sb_arr(i,j,k,CMASK_ID));
            if(!cmask) //for dead cells
            {
                rhs_arr(i,j,k)=0.0;
            }
        });

        for(int idim=0;idim<AMREX_SPACEDIM;idim++)
        {

            amrex::ParallelFor(face_boxes[idim], [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                IntVect face{AMREX_D_DECL(i,j,k)};
                IntVect lcell{AMREX_D_DECL(i,j,k)};
                IntVect rcell{AMREX_D_DECL(i,j,k)};
                lcell[idim]-=1;

                int mask_L,mask_R;

                if(!captured_conjsolve)
                {
                    mask_L=int(sb_arr(lcell,CMASK_ID));
                    mask_R=int(sb_arr(rcell,CMASK_ID));
                }
                else
                {
                    mask_L=int(1.0-sb_arr(lcell,CMASK_ID));
                    mask_R=int(1.0-sb_arr(rcell,CMASK_ID));
                }

                int cov_uncov_interface=(mask_L)*(!mask_R)+(!mask_L)*(mask_R);

                if(cov_uncov_interface) //1-0 or 0-1 interface
                {
                    int sgn=(int(sb_arr(lcell,CMASK_ID))==1)?1:-1;

                    tr_boundaries::bc_ib(face,idim,sgn,solved_comp,sb_arr,acoeff_arr,rhs_arr,
                                         domlo,domhi,prob_lo,prob_hi,dx,captured_time,*localprobparm,
                                         captured_conjsolve);
                }
            });
        }
    }
}

void TranspReact::set_solver_mask(Vector<iMultiFab>& solvermask,
                                  Vector<MultiFab>& neg_solvermask,
                                  Vector<MultiFab>& Sborder,int conjsolve)
{
    int captured_conjsolve=conjsolve;
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        for (MFIter mfi(Sborder[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Array4<int> smask_arr = solvermask[ilev].array(mfi);
            Array4<Real> neg_smask_arr=neg_solvermask[ilev].array(mfi);
            const Box& bx = mfi.tilebox();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                smask_arr(i,j,k)=(captured_conjsolve==0)?
                int(sb_arr(i,j,k,CMASK_ID)):(int(1.0-sb_arr(i,j,k,CMASK_ID)));
                neg_smask_arr(i,j,k)=1.0-smask_arr(i,j,k);
            });
        }
    }
}

void TranspReact::null_field_in_covered_cells(Vector<MultiFab>& fld,
                                              Vector<MultiFab>& Sborder,int startcomp,int numcomp)
{

    //multiply syntax
    //Multiply (FabArray<FAB>& dst, FabArray<FAB> const& src, 
    //int srccomp, int dstcomp, int numcomp, int nghost)

    for(int lev=0;lev<=finest_level;lev++)
    {
        for(int c=startcomp;c<(startcomp+numcomp);c++)
        {
            amrex::MultiFab::Multiply(fld[lev],Sborder[lev],CMASK_ID, c, 1, 0);
        }
    }
}

void TranspReact::smooth_cellmask()
{
    const Real time=0.0;
    amrex::Print()<<"finest_level:"<<finest_level<<"\t"<<phi_new.size()<<"\n";
    Vector<MultiFab> Sborder(finest_level+1);
    for(int lev=0;lev<=finest_level;lev++)
    {
        Sborder[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), ngrow_for_fillpatch);
        Sborder[lev].setVal(0.0);

        FillPatch(lev, time, Sborder[lev], 0, Sborder[lev].nComp());
    }

    for(int siter=0;siter<nsmoothiters;siter++)
    {
        for(int lev=0;lev<=finest_level;lev++)
        {
            const auto dx = geom[lev].CellSizeArray();
            auto prob_lo = geom[lev].ProbLoArray();
            auto prob_hi = geom[lev].ProbHiArray();
            const int* domlo_arr = geom[lev].Domain().loVect();
            const int* domhi_arr = geom[lev].Domain().hiVect();

            GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
            GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

            for (MFIter mfi(phi_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx, 1);

                Array4<Real> sborder_arr = Sborder[lev].array(mfi);
                Array4<Real> phi_arr = phi_new[lev].array(mfi);

                // update residual
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    //laplace smoothing
                    Real nadds=0.0;
                    Real neighborsum=0.0;
#if AMREX_SPACEDIM == 3
                    for(int kk=-1;kk<=1;kk++)
#endif
                    {
#if AMREX_SPACEDIM > 1
                        for(int jj=-1;jj<=1;jj++)
#endif
                        {
                            for(int ii=-1;ii<=1;ii++)
                            {
                                if(!(ii==0 && jj==0 && kk==0))
                                {
                                    nadds=nadds+1.0;
                                    neighborsum += sborder_arr(i+ii,j+jj,k+kk,CMASK_ID);
                                }
                            }
                        }
                    }

                    phi_arr(i,j,k,CMASK_ID)=neighborsum/nadds;
                });
            }

            //fill patching with new transformed variables        
            FillPatch(lev, time, Sborder[lev], 0, Sborder[lev].nComp());
        }

    }

    ProbParm* localprobparm = d_prob_parm;
    for(int lev=0;lev<=finest_level;lev++)
    {
        MultiFab& state = phi_new[lev];
        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            Array4<Real> fab = state[mfi].array();
            GeometryData geomData = geom[lev].data();
            const Box& box = mfi.validbox();

            amrex::launch(box, [=] AMREX_GPU_DEVICE(Box const& tbx) {
                    initdomaindata(tbx, fab, geomData, localprobparm);
                    });
        }
    
        amrex::MultiFab::Copy(phi_old[lev], phi_new[lev], 0, 0, phi_new[lev].nComp(), 0);
    }
}
