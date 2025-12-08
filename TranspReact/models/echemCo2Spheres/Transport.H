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
    amrex::Real get_elyte_conductivity(amrex::Real ci[NUM_TR_SPECIES],
                                  ProbParm const& prob_parm)
    {
        amrex::Real FsqbyRT=FARADCONST*FARADCONST/(GASCONST*prob_parm.T);
        amrex::Real cond=0.0;
        amrex::Real c0=1000.0;
        for(int sp=0;sp<NUM_TR_SPECIES;sp++)
        {
            if(amrex::Math::abs(prob_parm.Zi[sp])>0)
            {
                cond += 
                prob_parm.Zi[sp]*prob_parm.Zi[sp]*FsqbyRT*
                prob_parm.dcoeff[sp]*ci[sp]/c0;
            }
        }
        //cond=0.04;
        return(cond);
    }
    
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
            amrex::Real conc[NUM_TR_SPECIES];
            for(int sp=0;sp<NUM_TR_SPECIES;sp++)
            {
                conc[sp]=phi(i,j,k,sp);
            }
            dcoeff=(phi(i,j,k,CMASK_ID)>0.0)?
            get_elyte_conductivity(conc,prob_parm):prob_parm.cat_conductivity;
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
        Real dirvel=0.0; //only ambipolar diffusion
        //Real dcoeff=prob_parm.dcoeff[spec_id];
        //Real mob=dcoeff*prob_parm.Zi[spec_id]*FARADCONST/GASCONST/prob_parm.T;

        if(spec_id==Kp_ID || spec_id==OHm_ID)
        {
            amrex::Real conc[NUM_TR_SPECIES];
            for(int sp=0;sp<NUM_TR_SPECIES;sp++)
            {
                conc[sp]=phi(cellid,sp);
            }
            amrex::Real elyte_cond=(phi(cellid,CMASK_ID)>0.0)?
            get_elyte_conductivity(conc,prob_parm):0.0;

            dirvel=elyte_cond*phi(cellid,EFX_ID+dir)/FARADCONST
            *prob_parm.tplus[spec_id]/phi(cellid,spec_id);
        }
        return(dirvel);
    }
}
#endif
