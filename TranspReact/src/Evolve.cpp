#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <ProbParm.H>
#include <TranspReact.H>
#include <Species.H>
#include <AMReX_MLABecLaplacian.H>
#include <Utils.H>

// advance solution to final time
void TranspReact::Evolve_coupled()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    Real plottime = 0.0;
    Real chktime = 0.0;

    //there is a slight issue when restart file is not a multiple
    //a plot file may get the same number with an "old" file generated
    //note: if the user changes the chk_int and plot_int, they have to
    //manually set the old values for chk_int and plot_int,chk_int_old 
    //and plt_int_old, in the inputs, so that the offsets are correct 
    int plotfilenum=amrex::Math::floor(amrex::Real(istep[0])/amrex::Real(plot_int_old));
    int chkfilenum=amrex::Math::floor(amrex::Real(istep[0])/amrex::Real(chk_int_old));
    if(plot_time > 0.0) plotfilenum=amrex::Math::floor(amrex::Real(cur_time)/amrex::Real(plot_time));
    if(chk_time > 0.0) chkfilenum=amrex::Math::floor(amrex::Real(cur_time)/amrex::Real(chk_time));

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Real strt_time = amrex::second();
        amrex::Print() << "\nCoarse STEP " << step + 1 << " starts ..." << std::endl;
        amrex::Real dt_diff = std::numeric_limits<Real>::max();
        amrex::Real dt_adv = std::numeric_limits<Real>::max();
        amrex::Real dt_diff_lev,dt_adv_lev;

        for(int lev=0;lev<=finest_level;lev++)
        {
            find_transp_timescales(lev,cur_time,dt_diff_lev,dt_adv_lev);
            amrex::Print()<<"diffusion and adv time:"<<lev<<"\t"<<
            dt_diff_lev<<"\t"<<dt_adv_lev<<"\n";

            if(dt_diff_lev < dt_diff)
            {
                dt_diff = dt_diff_lev;
            }
            if(dt_adv_lev < dt_adv)
            {
                dt_adv = dt_adv_lev;
            }
        }

        ComputeDt(cur_time,dt_diff,dt_adv);

        if (max_level > 0 && regrid_int > 0)  // We may need to regrid
        {
            if (istep[0] % regrid_int == 0)
            {
                regrid(0, cur_time);
            }
        }

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
            amrex::Print() << "ADVANCE with time = " << t_new[lev]
            << " dt = " << dt[0] << std::endl;
        }
        amrex::Real dt_common=dt[0]; //no subcycling

        //ngrow fillpatch set in TranspReact.cpp
        //depending on hyperbolic order
        int num_grow=ngrow_for_fillpatch; 

        // Solution and sources MFs
        Vector<MultiFab> rxn_src(finest_level+1);
        Vector<MultiFab> Sborder(finest_level+1);
        Vector<MultiFab> Sborder_old(finest_level+1);
        Vector<MultiFab> phi_tmp(finest_level+1);
        Vector<MultiFab> adv_src(finest_level+1);

        //copy new to old and update time
        for(int lev=0;lev<=finest_level;lev++)
        {
            phi_tmp[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), phi_new[lev].nGrow());
            phi_tmp[lev].setVal(0.0);
            amrex::MultiFab::Copy(phi_old[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            amrex::MultiFab::Copy(phi_tmp[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            t_old[lev] = t_new[lev];
            t_new[lev] += dt_common;
        }

        //allocate flux, adv_src, Sborder
        for(int lev=0;lev<=finest_level;lev++)
        {
            Sborder[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder[lev].setVal(0.0);

            Sborder_old[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder_old[lev].setVal(0.0);

            FillPatch(lev, cur_time, Sborder_old[lev], 0, Sborder_old[lev].nComp());

            rxn_src[lev].define(grids[lev], dmap[lev], NUM_SPECIES, 0);
            rxn_src[lev].setVal(0.0);

            adv_src[lev].define(grids[lev], dmap[lev], 1, 0);
            adv_src[lev].setVal(0.0);
        }


        for(int niter=0;niter<num_timestep_correctors;niter++)
        {
            //for second order accuracy in mid point method
            amrex::Real time_offset=(niter>0)?0.5*dt_common:0.0;

            //reset all
            for(int lev=0;lev<=finest_level;lev++)
            {
                Sborder[lev].setVal(0.0);

                //grab phi_new all the time
                //at first iter phi new and old are same
                FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                rxn_src[lev].setVal(0.0);
            }
            
            //transform any variables
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }

            for(int lev=0;lev<=finest_level;lev++)
            {
                update_rxnsrc_at_level(lev, Sborder[lev], rxn_src[lev], cur_time+time_offset);
            }

            for(unsigned int ind=0;ind<NUM_SPECIES;ind++)
            {
                if(!unsolvedspec[ind])
                {
                    amrex::Print()<<"Solving species:"<<allvarnames[ind]<<"\n";
                    if(!conjugate_solve[ind])
                    {
                            if(do_advection)
                            {
                                update_advsrc_at_all_levels(ind, Sborder, adv_src, cur_time+time_offset,0);
                            }
                            implicit_solve_scalar(cur_time+time_offset, dt_common, ind, Sborder, Sborder_old, rxn_src, adv_src,0);
                    
                    }
                    else
                    {
                        for(int it=0;it<conjsolve_maxiter;it++)
                        {
                            for(int csolve=0;csolve<2;csolve++)
                            {
                                amrex::Print()<<"conjugate solve begins for spec "<<ind<<":"<<csolve<<"==============\n";
                                for(int lev=0;lev<=finest_level;lev++)
                                {
                                    Sborder[lev].setVal(0.0);

                                    //grab phi_new all the time
                                    //at first iter phi new and old are same
                                    FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                                }
                                
                                //update interface cells with neighbor averages
                                for(int iter=0;iter<interface_update_maxiter;iter++)
                                {
                                    for (int lev = 0; lev <= finest_level; lev++)
                                    {
                                        for (MFIter mfi(phi_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                                        {
                                            const Box& bx = mfi.tilebox();
                                            Array4<Real> phi_arr = phi_new[lev].array(mfi);
                                            Array4<Real> sb_arr = Sborder[lev].array(mfi);
                                            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                                update_interface_cells(i,j,k,ind,sb_arr,phi_arr);
                                            });
                                        }
                                    }
                                    for(int lev=0;lev<=finest_level;lev++)
                                    {
                                        Sborder[lev].setVal(0.0);

                                        //grab phi_new all the time
                                        //at first iter phi new and old are same
                                        FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                                    }
                                }
                                
                                if(do_advection)
                                {
                                    update_advsrc_at_all_levels(ind, Sborder, adv_src, cur_time+time_offset,csolve);
                                }
                                implicit_solve_scalar(cur_time+time_offset, dt_common, ind, Sborder, Sborder_old, rxn_src, adv_src,csolve);
                                amrex::Print()<<"conjugate solve for spec "<<ind<<" ends:"<<csolve<<"==============\n";
                            }
                        }
                    }
                }
            }

            if(niter<num_timestep_correctors-1)
            {
                //copy new to old and update time
                for(int lev=0;lev<=finest_level;lev++)
                {
                    amrex::Print()<<"averaging state at iter:"<<niter<<"\n";
                    MultiFab::LinComb(phi_tmp[lev], 0.5, phi_old[lev], 0, 0.5, 
                                      phi_new[lev], 0, 0, phi_new[lev].nComp(), 0);

                    //phi_{k+1}=0.5*(phi_n+phi_k)
                    amrex::MultiFab::Copy(phi_new[lev], phi_tmp[lev], 
                                          0, 0, phi_new[lev].nComp(), 0);
                }
            }
            amrex::Print()<<"\n================== Finished timestep iter:"<<niter+1<<" ================\n";
        }

        AverageDown ();

        for (int lev = 0; lev <= finest_level; lev++)
            ++istep[lev];

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }

        cur_time += dt_common;
        plottime += dt_common;
        chktime += dt_common;
        Real run_time = amrex::second() - strt_time;

        amrex::Print() << "Coarse STEP " << step + 1 << " ends."
        << " TIME = " << cur_time << " DT = " << dt_common << std::endl;
        amrex::Print()<<"Time step wall clock time:"<<run_time<<"\n";

        if (plot_time > 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            if(plottime > plot_time)
            {
                last_plot_file_step = step + 1;
                plotfilenum++;
                WritePlotFile(plotfilenum);
                plottime = 0.0;
            }
        }
        else if (plot_int > 0 && (step + 1) % plot_int == 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            last_plot_file_step = step + 1;
            plotfilenum++;
            WritePlotFile(plotfilenum);
        }

        if(chk_time > 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            if(chktime > chk_time)
            {
                chkfilenum++;
                amrex::Print()<<"writing chk file\n";
                WriteCheckpointFile(chkfilenum);
                chktime = 0.0;
            }
        }
        else if (chk_int > 0 && (step + 1) % chk_int == 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            amrex::Print()<<"writing chk file 1\n";
            chkfilenum++;
            WriteCheckpointFile(chkfilenum);
        }

        if (cur_time >= stop_time - 1.e-6 * dt_common) break;


        //local cleanup
        Sborder.clear();
        Sborder_old.clear();
        phi_tmp.clear();
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step)
    {
        plotfilenum++;
        WritePlotFile(plotfilenum);
    }
}

void TranspReact::Evolve_split()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    Real plottime = 0.0;
    Real chktime = 0.0;

    int plotfilenum=amrex::Math::floor(amrex::Real(istep[0])/amrex::Real(plot_int_old));
    int chkfilenum=amrex::Math::floor(amrex::Real(istep[0])/amrex::Real(chk_int_old));
    if(plot_time > 0.0) plotfilenum=amrex::Math::floor(amrex::Real(cur_time)/amrex::Real(plot_time));
    if(chk_time > 0.0) chkfilenum=amrex::Math::floor(amrex::Real(cur_time)/amrex::Real(chk_time));

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Real strt_time = amrex::second();
        amrex::Print() << "\nCoarse STEP " << step + 1 << " starts ..." << std::endl;
        amrex::Real dt_diff = std::numeric_limits<Real>::max();
        amrex::Real dt_adv = std::numeric_limits<Real>::max();
        amrex::Real dt_diff_lev,dt_adv_lev;

        for(int lev=0;lev<=finest_level;lev++)
        {
            find_transp_timescales(lev,cur_time,dt_diff_lev,dt_adv_lev);
            amrex::Print()<<"diffusion and adv time:"<<lev<<"\t"<<
            dt_diff_lev<<"\t"<<dt_adv_lev<<"\n";

            if(dt_diff_lev < dt_diff)
            {
                dt_diff = dt_diff_lev;
            }
            if(dt_adv_lev < dt_adv)
            {
                dt_adv = dt_adv_lev;
            }
        }

        ComputeDt(cur_time,dt_diff,dt_adv);

        if (max_level > 0 && regrid_int > 0)  // We may need to regrid
        {
            if (istep[0] % regrid_int == 0)
            {
                regrid(0, cur_time);
            }
        }

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
            amrex::Print() << "ADVANCE with time = " << t_new[lev]
            << " dt = " << dt[0] << std::endl;
        }
        amrex::Real dt_common=dt[0]; //no subcycling
        amrex::Real dt_common_inv=1.0/dt_common;

        //ngrow fillpatch set in TranspReact.cpp
        //depending on hyperbolic order
        int num_grow=ngrow_for_fillpatch; 

        // Solution and sources MFs
        Vector<MultiFab> rxn_src(finest_level+1);
        Vector<MultiFab> rxn_src_steady(finest_level+1);
        Vector<MultiFab> Sborder(finest_level+1);
        Vector<MultiFab> Sborder_old(finest_level+1);
        Vector<MultiFab> phi_tmp(finest_level+1);
        Vector<MultiFab> adv_src(finest_level+1);
        Vector<MultiFab> advdiff_src(finest_level+1);

        //copy new to old and update time
        for(int lev=0;lev<=finest_level;lev++)
        {
            phi_tmp[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), phi_new[lev].nGrow());
            phi_tmp[lev].setVal(0.0);
            amrex::MultiFab::Copy(phi_old[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            amrex::MultiFab::Copy(phi_tmp[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            t_old[lev] = t_new[lev];
            t_new[lev] += dt_common;
        }

        //allocate flux, adv_src, Sborder
        for(int lev=0;lev<=finest_level;lev++)
        {
            Sborder[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder[lev].setVal(0.0);

            Sborder_old[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder_old[lev].setVal(0.0);

            FillPatch(lev, cur_time, Sborder_old[lev], 0, Sborder_old[lev].nComp());

            rxn_src[lev].define(grids[lev], dmap[lev], NUM_SPECIES, 0);
            rxn_src[lev].setVal(0.0);
            
            rxn_src_steady[lev].define(grids[lev], dmap[lev], NUM_SPECIES, 0);
            rxn_src_steady[lev].setVal(0.0);

            adv_src[lev].define(grids[lev], dmap[lev], 1, 0);
            adv_src[lev].setVal(0.0);

            advdiff_src[lev].define(grids[lev], dmap[lev], NUM_SPECIES, 0);
            advdiff_src[lev].setVal(0.0);
        }



        for(int niter=0;niter<num_timestep_correctors;niter++)
        {
            //for second order accuracy in mid point method
            amrex::Real time_offset=(niter>0)?0.5*dt_common:0.0;


            //reset all
            for(int lev=0;lev<=finest_level;lev++)
            {
                Sborder[lev].setVal(0.0);

                //grab phi_new all the time
                //at first iter phi new and old are same
                FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                rxn_src[lev].setVal(0.0);
                rxn_src_steady[lev].setVal(0.0);
            }
            
            for(int inner_iter=0;inner_iter<num_split_correctors;inner_iter++)
            {
                for(int lev=0;lev<=finest_level;lev++)
                {
                    chemistry_advance(lev,cur_time,dt_common,advdiff_src[lev],phi_old[lev],phi_tmp[lev]);
                }
                for(int lev=0;lev<=finest_level;lev++)
                {
                    MultiFab::LinComb(rxn_src[lev], dt_common_inv, phi_tmp[lev], 0, -dt_common_inv, 
                                      phi_old[lev], 0, 0, NUM_SPECIES, 0);

                    amrex::MultiFab::Saxpy(rxn_src[lev], -1.0, advdiff_src[lev], 0, 0, NUM_SPECIES, 0);
                }
            
                
                //transform any variables
                if(transform_vars)
                {
                    //sborder old is already with phi_new
                    transform_variables(Sborder,cur_time,dt_common);
                }

                for(unsigned int ind=0;ind<NUM_SPECIES;ind++)
                {
                    if(!unsolvedspec[ind])
                    {
                        if(!steadyspec[ind])
                        {
                            if(!conjugate_solve[ind])
                            {
                                if(do_advection)
                                {
                                    update_advsrc_at_all_levels(ind, Sborder, adv_src, cur_time+time_offset,0);
                                }
                                implicit_solve_scalar(cur_time+time_offset, dt_common, ind, Sborder, Sborder_old, rxn_src, adv_src,0);
                            }
                            else
                            {
                                for(int it=0;it<conjsolve_maxiter;it++)
                                {
                                    for(int csolve=0;csolve<2;csolve++)
                                    {
                                        for(int lev=0;lev<=finest_level;lev++)
                                        {
                                            Sborder[lev].setVal(0.0);

                                            //grab phi_new all the time
                                            //at first iter phi new and old are same
                                            FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                                        }

                                        //update interface cells with neighbor averages
                                        for(int iter=0;iter<interface_update_maxiter;iter++)
                                        {
                                            for (int lev = 0; lev <= finest_level; lev++)
                                            {
                                                for (MFIter mfi(phi_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                                                {
                                                    const Box& bx = mfi.tilebox();
                                                    Array4<Real> phi_arr = phi_new[lev].array(mfi);
                                                    Array4<Real> sb_arr = Sborder[lev].array(mfi);
                                                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                                        update_interface_cells(i,j,k,ind,sb_arr,phi_arr);
                                                    });
                                                }
                                            }
                                            for(int lev=0;lev<=finest_level;lev++)
                                            {
                                                Sborder[lev].setVal(0.0);

                                                //grab phi_new all the time
                                                //at first iter phi new and old are same
                                                FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                                            }
                                        }
                                        if(do_advection)
                                        {
                                            update_advsrc_at_all_levels(ind, Sborder, adv_src, cur_time+time_offset,csolve);
                                        }
                                        implicit_solve_scalar(cur_time+time_offset, dt_common, ind, Sborder, Sborder_old, rxn_src, adv_src,csolve);
                                    }
                                }
                            }
                        }
                    }
                }
                for(int lev=0;lev<=finest_level;lev++)
                {
                    MultiFab::LinComb(advdiff_src[lev], dt_common_inv, phi_new[lev], 0, -dt_common_inv, 
                                      phi_old[lev], 0, 0, NUM_SPECIES, 0);

                    amrex::MultiFab::Saxpy(advdiff_src[lev], -1.0, rxn_src[lev], 0, 0, NUM_SPECIES, 0);
                }
            }

            for(unsigned int ind=0;ind<NUM_SPECIES;ind++)
            {
                if(!unsolvedspec[ind] && steadyspec[ind])
                {
                    amrex::Print()<<"Solving species:"<<allvarnames[ind]<<"\n";
                    if(!conjugate_solve[ind])
                    {
                        if(do_advection)
                        {
                            update_advsrc_at_all_levels(ind, Sborder, adv_src, cur_time+time_offset,0);
                        }
                        
                        for(int lev=0;lev<=finest_level;lev++)
                        {
                            int only_steady_sources=1;
                            update_rxnsrc_at_level(lev, Sborder[lev], rxn_src_steady[lev], cur_time+time_offset, only_steady_sources);
                        }
                        implicit_solve_scalar(cur_time+time_offset, dt_common, ind, Sborder, Sborder_old, rxn_src_steady, adv_src,0);
                    }
                    else
                    {
                        for(int it=0;it<conjsolve_maxiter;it++)
                        {
                            for(int csolve=0;csolve<2;csolve++)
                            {
                                amrex::Print()<<"conjugate solve begins for spec "<<ind<<":"<<csolve<<"==============\n";
                                for(int lev=0;lev<=finest_level;lev++)
                                {
                                    Sborder[lev].setVal(0.0);

                                    //grab phi_new all the time
                                    //at first iter phi new and old are same
                                    FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                                }

                                //update interface cells with neighbor averages
                                for(int iter=0;iter<interface_update_maxiter;iter++)
                                {
                                    for (int lev = 0; lev <= finest_level; lev++)
                                    {
                                        for (MFIter mfi(phi_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                                        {
                                            const Box& bx = mfi.tilebox();
                                            Array4<Real> phi_arr = phi_new[lev].array(mfi);
                                            Array4<Real> sb_arr = Sborder[lev].array(mfi);
                                            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                                update_interface_cells(i,j,k,ind,sb_arr,phi_arr);
                                            });
                                        }
                                    }
                                    for(int lev=0;lev<=finest_level;lev++)
                                    {
                                        Sborder[lev].setVal(0.0);

                                        //grab phi_new all the time
                                        //at first iter phi new and old are same
                                        FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                                    }
                                }
                                if(do_advection)
                                {
                                    update_advsrc_at_all_levels(ind, Sborder, adv_src, cur_time+time_offset,csolve);
                                }
                                for(int lev=0;lev<=finest_level;lev++)
                                {
                                    int only_steady_sources=1;
                                    update_rxnsrc_at_level(lev, Sborder[lev], rxn_src_steady[lev], cur_time+time_offset, only_steady_sources);
                                }
                                implicit_solve_scalar(cur_time+time_offset, dt_common, ind, Sborder, Sborder_old, rxn_src, adv_src,csolve);
                                amrex::Print()<<"conjugate solve for spec "<<ind<<" ends:"<<csolve<<"==============\n";
                            }
                        }
                    }
                }
            }

            if(niter<num_timestep_correctors-1)
            {
                //copy new to old and update time
                for(int lev=0;lev<=finest_level;lev++)
                {
                    amrex::Print()<<"averaging state at iter:"<<niter<<"\n";
                    MultiFab::LinComb(phi_tmp[lev], 0.5, phi_old[lev], 0, 0.5, 
                                      phi_new[lev], 0, 0, phi_new[lev].nComp(), 0);

                    //phi_{k+1}=0.5*(phi_n+phi_k)
                    amrex::MultiFab::Copy(phi_new[lev], phi_tmp[lev], 
                                          0, 0, phi_new[lev].nComp(), 0);
                }
            }
            amrex::Print()<<"\n================== Finished timestep iter:"<<niter+1<<" ================\n";
        }

        AverageDown ();

        for (int lev = 0; lev <= finest_level; lev++)
            ++istep[lev];

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }

        cur_time += dt_common;
        plottime += dt_common;
        chktime += dt_common;
        Real run_time = amrex::second() - strt_time;

        amrex::Print() << "Coarse STEP " << step + 1 << " ends."
        << " TIME = " << cur_time << " DT = " << dt_common << std::endl;
        amrex::Print()<<"Time step wall clock time:"<<run_time<<"\n";

        if (plot_time > 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            if(plottime > plot_time)
            {
                last_plot_file_step = step + 1;
                plotfilenum++;
                WritePlotFile(plotfilenum);
                plottime = 0.0;
            }
        }
        else if (plot_int > 0 && (step + 1) % plot_int == 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            last_plot_file_step = step + 1;
            plotfilenum++;
            WritePlotFile(plotfilenum);
        }

        if(chk_time > 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            if(chktime > chk_time)
            {
                chkfilenum++;
                amrex::Print()<<"writing chk file\n";
                WriteCheckpointFile(chkfilenum);
                chktime = 0.0;
            }
        }
        else if (chk_int > 0 && (step + 1) % chk_int == 0)
        {
            if(transform_vars)
            {
                //sborder old is already with phi_new
                transform_variables(Sborder,cur_time,dt_common);
            }
            amrex::Print()<<"writing chk file 1\n";
            chkfilenum++;
            WriteCheckpointFile(chkfilenum);
        }

        if (cur_time >= stop_time - 1.e-6 * dt_common) break;


        //local cleanup
        Sborder.clear();
        Sborder_old.clear();
        phi_tmp.clear();
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step)
    {
        plotfilenum++;
        WritePlotFile(plotfilenum);
    }

}
