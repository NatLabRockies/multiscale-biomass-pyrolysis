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
            specsource(i,j,k,S1_ID) = -10.0*phi(i,j,k,S1_ID);
            specsource(i,j,k,S2_ID) = -5.0*phi(i,j,k,S2_ID);
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
