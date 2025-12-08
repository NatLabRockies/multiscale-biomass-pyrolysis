#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>
#include <ProbParm.H>

using namespace amrex;
namespace tr_reactions
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void production_rate(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& specsource,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time, ProbParm const& prob_parm)
        {
            IntVect iv{AMREX_D_DECL(i,j,k)};
            amrex::Real d=prob_parm.axdir;
            amrex::Real x=prob_lo[d]+(iv[d]+0.5)*dx[d];

            for(int c=0;c<NUM_SPECIES;c++)
            {
                specsource(i,j,k,c)=0.0;
            }
            specsource(i,j,k,S1_ID) = 5.0*std::pow(x,4.0)/3.0-x/6.0-2.0;
            specsource(i,j,k,PHI_ID) = -phi(i,j,k,S1_ID); 
        }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    void transform(int i, int j, int k,
                   Array4<Real> const& s_arr,
                   Array4<Real> const& state_arr,
                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                   GpuArray<int, AMREX_SPACEDIM> domlo,
                   GpuArray<int, AMREX_SPACEDIM> domhi,
                   GpuArray<Real, AMREX_SPACEDIM> dx,
                   const Real time, ProbParm const& prob_parm)
    {
        for(int dir=0;dir<AMREX_SPACEDIM;dir++)
        {
            IntVect cellid{AMREX_D_DECL(i, j, k)};

            if(cellid[dir]>domlo[dir] && cellid[dir]<domhi[dir])
            {
                IntVect ip1{AMREX_D_DECL(i,j,k)};
                IntVect im1{AMREX_D_DECL(i,j,k)};

                ip1[dir]+=1;
                im1[dir]-=1;

                state_arr(cellid,EFX_ID+dir)=-0.5*(s_arr(ip1,PHI_ID)-s_arr(im1,PHI_ID))/dx[dir];
            }
            else if(cellid[dir]==domlo[dir])
            {
                IntVect ip1{AMREX_D_DECL(i,j,k)};
                IntVect ip2{AMREX_D_DECL(i,j,k)};
                ip1[dir]+=1;
                ip2[dir]+=2;

                //fdash = (-f_{i+2}+4f_{i+1}-3f_i)/(2 dx)
                if(domhi[dir] - domlo[dir] > 2)
                {
                    state_arr(cellid,EFX_ID+dir)=-0.5*(  -s_arr(ip2,PHI_ID) 
                                     + 4.0*s_arr(ip1,PHI_ID)
                                     - 3.0*s_arr(cellid,PHI_ID))/dx[dir];
                } else 
                {
                    state_arr(cellid,EFX_ID+dir)=-(s_arr(ip1,PHI_ID) - s_arr(cellid,PHI_ID))/dx[dir];
                }
            }
            else if(cellid[dir]==domhi[dir])
            {
                IntVect im1{AMREX_D_DECL(i,j,k)};
                IntVect im2{AMREX_D_DECL(i,j,k)};
                im1[dir]-=1;
                im2[dir]-=2;

                //fdash = (f_{i-2}-4f_{i-1}+3f_i)/(2 dx)
                if(domhi[dir] - domlo[dir] > 2){
                    state_arr(cellid,EFX_ID+dir)=-0.5*(     s_arr(im2,PHI_ID)
                                     - 4.0*s_arr(im1,PHI_ID)
                                     + 3.0*s_arr(cellid,PHI_ID))/dx[dir];
                } 
                else 
                {
                    state_arr(cellid,EFX_ID+dir)=-(s_arr(cellid,PHI_ID) - s_arr(im1,PHI_ID))/dx[dir];
                }
            }
            else
            {
                //there shouldnt be another option!
            }
        }

    }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_div_sources(IntVect faceid, int dir,
                             int spec_id,
                             Array4<Real> const& sb_arr,
                             Array4<Real> const& flx_arr,
                             const GpuArray<int, AMREX_SPACEDIM> bclo,
                             const GpuArray<int, AMREX_SPACEDIM> bchi,
                             const GpuArray<int, AMREX_SPACEDIM> domlo,
                             const GpuArray<int, AMREX_SPACEDIM> domhi,
                             const Real time,
                             const Real dt,
                             GpuArray<Real, AMREX_SPACEDIM> dx,
                             ProbParm const& prob_parm,
                             int using_ib,
                             int conjsolve)

    {

    }
}
#endif
