#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Chemistry.H>
#include <ProbParm.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;

AMREX_INLINE
void amrex_probinit (ProbParm & h_prob_parm, ProbParm & d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");
    pp.query("dcoeff_etrode", h_prob_parm.dc_etrode);
    pp.query("dcoeff_elyte", h_prob_parm.dc_elyte);
    pp.query("x1",h_prob_parm.x1);
    pp.query("phiL",h_prob_parm.phiL);
    pp.query("phiR",h_prob_parm.phiR);
    pp.query("left_etrode",h_prob_parm.left_etrode);
    pp.query("axial_dir",h_prob_parm.axial_dir);

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#endif

}

AMREX_GPU_DEVICE
AMREX_INLINE
void initdomaindata(Box const& bx, Array4<Real> const& phi, 
        GeometryData const& geomdata,
        ProbParm *prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();
    int axdir=prob_parm->axial_dir;
    amrex::Real axlen=prob_hi[axdir]-prob_lo[axdir];
    amrex::Real x1loc=prob_lo[axdir]+prob_parm->x1*axlen;
    
#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            Real z = prob_lo[2] + (0.5+k) * dx[2];
            Real y = prob_lo[1] + (0.5+j) * dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                IntVect iv(AMREX_D_DECL(i,j,k));
                Real x = prob_lo[0] + (0.5+i) * dx[0];

                Real xleft  = prob_lo[axdir]+iv[axdir]*dx[axdir];
                Real xright = prob_lo[axdir]+(iv[axdir]+1)*dx[axdir];
                      
                int in_etrode_l=0;
                int in_etrode_r=0;

                if(prob_parm->left_etrode)
                {
                    in_etrode_l=xleft<x1loc;
                    in_etrode_r=xright<x1loc;
                }
                else
                {
                    in_etrode_l=xleft>x1loc;
                    in_etrode_r=xright>x1loc;
                }

                int frac_cell=(in_etrode_l)*(!in_etrode_r);

                phi(i,j,k,ETRODE_ID)=(in_etrode_l && in_etrode_r)?1.0:0.0;

                //deal with fractional cells
                if(frac_cell)
                {
                    phi(i,j,k,ETRODE_ID)=(x1loc-xleft)/dx[axdir];
                }

                phi(i,j,k,CMASK_ID)=1.0-phi(i,j,k,ETRODE_ID);
                phi(i,j,k,ELYTE_ID)=1.0-phi(i,j,k,ETRODE_ID);
                phi(i,j,k,POT_ID)  = 1.0;
            }
        }
    }
}
#endif
