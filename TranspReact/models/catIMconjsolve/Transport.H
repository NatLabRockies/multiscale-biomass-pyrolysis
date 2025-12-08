#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>
#include <ProbParm.H>

using namespace amrex;
namespace tr_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
    Real specDiff(int i, int j, int k,
                  int spec_id,
                  Array4<Real> const& phi,
                  GpuArray<Real, AMREX_SPACEDIM> dx,
                  GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                  GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                  const Real time, ProbParm const& prob_parm)
    {
        Real memfrac=phi(i,j,k,MEM_ID);
        Real imfrac=phi(i,j,k,IM_ID);
        Real liqfrac=phi(i,j,k,LE_ID);

        Real dcoeff=1.0;

        if(spec_id!=PHI_ID)
        {
            dcoeff=prob_parm.dcoeff[spec_id];
        }
        else
        {
            //1 in electrolyte and
            //catalyst conductivity in electrode
            dcoeff=(phi(i,j,k,CMASK_ID)>0.0)?prob_parm.dcoeff[spec_id]:prob_parm.cat_conductivity;
        }
        return(dcoeff);
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_vel(IntVect cellid, int dir,
                            int spec_id,
                            Array4<Real> const& phi,
                            GpuArray<Real, AMREX_SPACEDIM> dx,
                            const Real time,
                            ProbParm const& prob_parm)

    {
        Real dcoeff=prob_parm.dcoeff[spec_id];
        Real mob=dcoeff*prob_parm.Zi[spec_id]*FARADCONST/GASCONST/prob_parm.T;
        amrex::Real dirvel=mob*phi(cellid,EFX_ID+dir);
        return(dirvel);
    }
}
#endif
