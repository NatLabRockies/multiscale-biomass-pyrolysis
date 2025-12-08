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
    pp.query("dcoeff_mat1", h_prob_parm.dc_mat1);
    pp.query("dcoeff_mat2", h_prob_parm.dc_mat2 );
    pp.query("dcoeff_mat3", h_prob_parm.dc_mat3 );
    pp.query("x1",h_prob_parm.x1);
    pp.query("x2",h_prob_parm.x2);
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
    amrex::Real x2loc=prob_lo[axdir]+prob_parm->x2*axlen;
    
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
                IntVect iv(i,j,k);
                Real x = prob_lo[0] + (0.5+i) * dx[0];

                Real xleft  = prob_lo[axdir]+iv[axdir]*dx[axdir];
                Real xright = prob_lo[axdir]+(iv[axdir]+1)*dx[axdir];
                        
                int in_mat1_l=xleft<x1loc;
                int in_mat1_r=xright<x1loc;
                
                int in_mat2_l=xleft>=x1loc && xleft<x2loc;
                int in_mat2_r=xright>=x1loc && xright<x2loc;
                
                int in_mat3_l=xleft>=x2loc;
                int in_mat3_r=xright>=x2loc;

                int frac_cell_mat12=(in_mat1_l)*(!in_mat1_r);
                int frac_cell_mat23=(!in_mat3_l)*(in_mat3_r);

                phi(i,j,k,MAT1_ID)=(in_mat1_l && in_mat1_r)?1.0:0.0;
                phi(i,j,k,MAT2_ID)=(in_mat2_l && in_mat2_r)?1.0:0.0;
                phi(i,j,k,MAT3_ID)=(in_mat3_l && in_mat3_r)?1.0:0.0;

                //deal with fractional cells
                if(frac_cell_mat12)
                {
                   phi(i,j,k,MAT1_ID)=(x1loc-xleft)/dx[axdir];
                   phi(i,j,k,MAT2_ID)=1.0-phi(i,j,k,MAT1_ID);
                }
                
                if(frac_cell_mat23)
                {
                   phi(i,j,k,MAT2_ID)=(x2loc-xleft)/dx[axdir];
                   phi(i,j,k,MAT3_ID)=1.0-phi(i,j,k,MAT2_ID);
                }

                phi(i,j,k,CMASK_ID)=1.0-phi(i,j,k,MAT1_ID);

                if(phi(i,j,k,CMASK_ID) < 1.0 && phi(i,j,k,CMASK_ID) > 0.0)
                {
                    //dummy value in fractional cmask cell to
                    //see if it is getting updated
                    phi(i,j,k,TEMP_ID)  = -1.0;
                }
                else
                {
                    phi(i,j,k,TEMP_ID)  = 0.5;
                }
            }
        }
    }
}
#endif
