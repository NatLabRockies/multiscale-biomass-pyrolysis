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
    //pp.query("right_temp", h_prob_parm.right_temp );
    pp.query("Kmem", h_prob_parm.Kmem );
    pp.query("Kim", h_prob_parm.Kim );
    pp.query("Kcat", h_prob_parm.Kcat );
    pp.query("Kle", h_prob_parm.Kle );
    pp.query("sourcefac", h_prob_parm.sourcefac );
    pp.query("tmap_init",h_prob_parm.tmap_init);
    pp.query("koh_init",h_prob_parm.koh_init);
    pp.query("co2_boundary",h_prob_parm.co2_boundary);
    pp.query("catelecpot",h_prob_parm.catelecpot);
    pp.query("i0_co",h_prob_parm.i0_co);
    pp.query("i0_h2",h_prob_parm.i0_h2);
    pp.query("E0_co",h_prob_parm.E0_co);
    pp.query("E0_h2",h_prob_parm.E0_h2);

    for (int c=0;c<NUM_SPECIES;c++)
    {
        h_prob_parm.dcoeff[c]=0.0;
        h_prob_parm.Zi[c]=0.0;
    }

    //all in m2/s
    h_prob_parm.dcoeff[PHI_ID]=1.0;
    h_prob_parm.dcoeff[CO2_ID]=4.3716e-9;
    h_prob_parm.dcoeff[CO_ID]=4.3716e-9;
    h_prob_parm.dcoeff[Hp_ID]=1.156e-8;
    h_prob_parm.dcoeff[OHm_ID]=9.1624e-9;
    h_prob_parm.dcoeff[HCO3m_ID]=2.3463e-9;
    h_prob_parm.dcoeff[CO32m_ID]=1.6774e-9;
    h_prob_parm.dcoeff[HCOOm_ID]=1.493e-9;
    h_prob_parm.dcoeff[HCOOH_ID]=1.493e-9;
    h_prob_parm.dcoeff[Kp_ID]=4.425e-9;
    h_prob_parm.dcoeff[TMAp_ID]=1e-20;

    //charge number
    h_prob_parm.Zi[Hp_ID]=1.0;
    h_prob_parm.Zi[OHm_ID]=-1.0;
    h_prob_parm.Zi[HCO3m_ID]=-1.0;
    h_prob_parm.Zi[CO32m_ID]=-2.0;
    h_prob_parm.Zi[HCOOm_ID]=-1.0;
    h_prob_parm.Zi[Kp_ID]=1.0;
    h_prob_parm.Zi[TMAp_ID]=1.0;

    
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
    
    const Real Lx = prob_hi[0] - prob_lo[0];
    const Real Ly = prob_hi[1] - prob_lo[1];
    const Real Lz = prob_hi[2] - prob_lo[2];

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                for(int c=0;c<NVAR;c++)
                {
                    phi(i,j,k,c)=0.0;
                }
                phi(i,j,k,CMASK_ID)=1.0;
                for(int kk=0;kk<2;kk++)
                {
                    for(int jj=0;jj<2;jj++)
                    {
                        for(int ii=0;ii<2;ii++)
                        {
                            Real x[3]={prob_lo[0]+(i+ii)*dx[0],
                                prob_lo[1]+(j+jj)*dx[1],           
                                prob_lo[2]+(k+kk)*dx[2]};


                            if(x[0]>=prob_lo[0]+0.25*Lx &&
                               x[1]<=prob_lo[1]+0.5*Ly &&
                               x[2]>=prob_lo[2]+0.5*Lz)
                            {
                                phi(i,j,k,CAT_ID)+=1.0;
                            }
                            else if(x[0]<prob_lo[0]+0.25*Lx)
                            {
                                phi(i,j,k,MEM_ID)+=1.0; 
                            }
                            else if(x[0]>=prob_lo[0]+0.25*Lx &&
                                    x[0]<prob_lo[0]+0.5*Lx)
                            {
                                phi(i,j,k,IM_ID)+=1.0; 
                            }
                            else
                            {
                                    phi(i,j,k,LE_ID)+=1.0;
                            }
                        }
                    }
                }
                phi(i,j,k,IM_ID)*=0.125;
                phi(i,j,k,MEM_ID)*=0.125;
                phi(i,j,k,CAT_ID)*=0.125;
                phi(i,j,k,LE_ID)*=0.125;
                phi(i,j,k,CMASK_ID)=(phi(i,j,k,CAT_ID)<1.0)?1.0:0.0;
                phi(i,j,k,TMAp_ID)=(phi(i,j,k,IM_ID)+phi(i,j,k,MEM_ID))*prob_parm->tmap_init;
                phi(i,j,k,OHm_ID)=(phi(i,j,k,IM_ID)+phi(i,j,k,MEM_ID))*prob_parm->tmap_init
                +phi(i,j,k,LE_ID)*prob_parm->koh_init;
                phi(i,j,k,Kp_ID)=phi(i,j,k,LE_ID)*prob_parm->koh_init;
            }
        }
    }
}
#endif
