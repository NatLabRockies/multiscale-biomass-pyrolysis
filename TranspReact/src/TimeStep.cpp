#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Prob.H>
#include <Transport.H>
#include <TranspReact.H>

// a wrapper for EstTimeStep
void TranspReact::ComputeDt(amrex::Real cur_time,amrex::Real dt_diff, amrex::Real dt_adv)
{
    if(adaptive_dt)
    {
        amrex::Real adaptive_dt = std::min(dt_adv*advective_cfl, dt_diff*diffusive_cfl);
        if(adaptive_dt > dt_max) 
        {
            adaptive_dt = dt_max;
        }
        if(adaptive_dt < dt_min) 
        {
            adaptive_dt = dt_min;
        }
        dt[0] = adaptive_dt;
    }
    else
    {
        dt[0]=fixed_dt;
    }
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        dt[lev] = dt[lev - 1];
    }
}


void TranspReact::find_transp_timescales(int lev,amrex::Real cur_time,
                                         amrex::Real &dt_diff,amrex::Real &dt_adv)
{
    BL_PROFILE("Vidyut::EstTimeStep()");
    //1e-5 is to accomodate a cfl mulitply
    dt_diff = 1e-5*std::numeric_limits<Real>::max(); 
    dt_adv = 1e-5*std::numeric_limits<Real>::max();
    ProbParm const* localprobparm = d_prob_parm;
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();

    GpuArray<int,NUM_SPECIES> steady_species;
    GpuArray<int,NUM_SPECIES> unsolved_species;

    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        steady_species[sp]=steadyspec[sp];
        unsolved_species[sp]=unsolvedspec[sp];
    }

    MultiFab& S_new = phi_new[lev];
    Real time=cur_time;
    int do_adv=do_advection;

    MultiFab vel(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab diff(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> state_array = S_new.array(mfi);
            Array4<Real> diff_array = diff.array(mfi);
            Array4<Real> vel_array = vel.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                amrex::Real maxdcoeff=0.0;
                amrex::Real maxvel=0.0;
                for(int c=0;c<NUM_SPECIES;c++)
                {
                    amrex::Real dcoeff=0.0;
                    amrex::Real specvel=0.0;

                    if(!steady_species[c] && !unsolved_species[c])
                    {
                        dcoeff=tr_transport::specDiff(i,j,k,c,state_array,
                                                      dx,prob_lo,prob_hi,time,*localprobparm); 

                        specvel=tr_transport::specDiff(i,j,k,c,state_array,
                                                       dx,prob_lo,prob_hi,time,*localprobparm); 
                    }

                    if(amrex::Math::abs(dcoeff)>maxdcoeff)
                    {
                        maxdcoeff=amrex::Math::abs(dcoeff);
                    }
                    if(do_adv)
                    {
                        IntVect iv_cell{AMREX_D_DECL(i, j, k)};
                        amrex::Real vx=0.0;
                        amrex::Real vy=0.0;
                        amrex::Real vz=0.0;

                        vx=tr_transport::compute_vel(iv_cell,0,c,state_array,dx,time,*localprobparm);
#if AMREX_SPACEDIM > 1
                        vy=tr_transport::compute_vel(iv_cell,1,c,state_array,dx,time,*localprobparm);
#if AMREX_SPACEDIM == 3
                        vz=tr_transport::compute_vel(iv_cell,2,c,state_array,dx,time,*localprobparm);
#endif
#endif
                        amrex::Real specvel=std::sqrt(vx*vx+vy*vy+vz*vz);

                        if(specvel>maxvel)
                        {
                            maxvel=specvel;
                        }
                    }
                }
                diff_array(i,j,k)=maxdcoeff;
                vel_array(i,j,k)=maxvel;
            });
        }
    }

    amrex::Real max_diff    = diff.norm0(0,0,true);
    amrex::Print()<<"max_diff:"<<max_diff<<"\n";
    if(max_diff > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
        {
            dt_diff = std::min(dt_diff, 0.5*dx[i]*dx[i]/max_diff/AMREX_SPACEDIM);
        }
    }
    if(do_adv)
    {
        amrex::Real max_vel    = vel.norm0(0,0,true);
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
        {
            dt_adv = std::min(dt_adv, dx[i]/max_vel);
        }
    }

    ParallelDescriptor::ReduceRealMin(dt_diff);
    ParallelDescriptor::ReduceRealMin(dt_adv);
}
