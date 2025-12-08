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
    pp.query("tmap_conc",h_prob_parm.tmap_conc);
    pp.query("koh_init",h_prob_parm.koh_init);
    pp.query("co2_boundary",h_prob_parm.co2_boundary);
    pp.query("catelecpot",h_prob_parm.catelecpot);
    pp.query("i0_co",h_prob_parm.i0_co);
    pp.query("i0_h2",h_prob_parm.i0_h2);
    pp.query("E0_co",h_prob_parm.E0_co);
    pp.query("E0_h2",h_prob_parm.E0_h2);
    pp.query("im_thickness",h_prob_parm.im_thickness);
    pp.query("imfrac",h_prob_parm.imfrac);
    pp.query("catfrac",h_prob_parm.catfrac);
    pp.query("current_input",h_prob_parm.current_input);

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
    /*h_prob_parm.dcoeff[Hp_ID]=3e-9;
    h_prob_parm.dcoeff[OHm_ID]=3e-9;
    h_prob_parm.dcoeff[HCO3m_ID]=3e-9;
    h_prob_parm.dcoeff[CO32m_ID]=3e-9;
    h_prob_parm.dcoeff[HCOOm_ID]=3e-9;
    h_prob_parm.dcoeff[HCOOH_ID]=3e-9;
    h_prob_parm.dcoeff[Kp_ID]=3e-9;
    h_prob_parm.dcoeff[TMAp_ID]=1e-20;*/

    //charge number
    h_prob_parm.Zi[Hp_ID]=1.0;
    h_prob_parm.Zi[OHm_ID]=-1.0;
    h_prob_parm.Zi[HCO3m_ID]=-1.0;
    h_prob_parm.Zi[CO32m_ID]=-2.0;
    h_prob_parm.Zi[HCOOm_ID]=-1.0;
    h_prob_parm.Zi[Kp_ID]=1.0;
    h_prob_parm.Zi[TMAp_ID]=1.0;

    amrex::Vector<amrex::Real> locx;
    amrex::Vector<amrex::Real> locy;
    amrex::Vector<amrex::Real> locz;
    amrex::Vector<amrex::Real> imcoverage;

    pp.query("nspheres",h_prob_parm.nsph);

    locx.resize(h_prob_parm.nsph);
    locy.resize(h_prob_parm.nsph);
    locz.resize(h_prob_parm.nsph);

    pp.queryarr("loc_x",locx);
    pp.queryarr("loc_y",locy);
    pp.queryarr("loc_z",locz);
    pp.queryarr("imcoverage",imcoverage);
    pp.query("sphererad",h_prob_parm.sphererad);

    h_prob_parm.loc_x  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
    h_prob_parm.loc_y  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
    h_prob_parm.loc_z  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
    h_prob_parm.imcoverage  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                locx.begin(),
                                locx.end(),
                       h_prob_parm.loc_x);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                locy.begin(),
                                locy.end(),
                       h_prob_parm.loc_y);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                locz.begin(),
                                locz.end(),
                       h_prob_parm.loc_z);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                imcoverage.begin(),
                                imcoverage.end(),
                       h_prob_parm.imcoverage);

    
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
    
    int axdir=0;
    amrex::Real cat_interface_loc=prob_lo[axdir]+prob_parm->catfrac*
    (prob_hi[axdir]-prob_lo[axdir]);
    amrex::Real initial_currentden=0.0;

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
                IntVect iv(AMREX_D_DECL(i,j,k));
                for(int c=0;c<NVAR;c++)
                {
                    phi(i,j,k,c)=0.0;
                }
                Real xleft  = prob_lo[axdir]+iv[axdir]*dx[axdir];
                Real xright = prob_lo[axdir]+(iv[axdir]+1)*dx[axdir];
                Real x=0.5*(xleft+xright);
                
                int in_cat_l=0;
                int in_cat_r=0;
                    
                in_cat_l=xleft<cat_interface_loc;
                in_cat_r=xright<cat_interface_loc;
                
                int frac_cell=(in_cat_l)*(!in_cat_r);
                phi(i,j,k,CAT_ID)=(in_cat_l && in_cat_r)?1.0:0.0;
                
                //deal with fractional cells
                if(frac_cell)
                {
                    phi(i,j,k,CAT_ID)=(cat_interface_loc-xleft)/dx[axdir];
                }
                
                phi(i,j,k,LE_ID)=(1.0-phi(i,j,k,CAT_ID));
                phi(i,j,k,CMASK_ID)=phi(i,j,k,LE_ID);
                /*phi(i,j,k,TMAp_ID)=(phi(i,j,k,IM_ID)+phi(i,j,k,MEM_ID))*prob_parm->tmap_conc;
                phi(i,j,k,OHm_ID)=(phi(i,j,k,IM_ID)+phi(i,j,k,MEM_ID))*prob_parm->tmap_conc
                +phi(i,j,k,LE_ID)*prob_parm->koh_init;
                phi(i,j,k,Kp_ID)=phi(i,j,k,LE_ID)*prob_parm->koh_init;*/
                phi(i,j,k,CO2_ID)=phi(i,j,k,LE_ID)*prob_parm->co2_boundary;
                phi(i,j,k,OHm_ID)=phi(i,j,k,LE_ID)*prob_parm->koh_init*0.000001;
                phi(i,j,k,Kp_ID)=phi(i,j,k,LE_ID)*prob_parm->koh_init*0.0000001;

                //solid potential
                if(phi(i,j,k,CAT_ID)==1.0)
                {
                    phi(i,j,k,PHI_ID)=prob_parm->catelecpot
                    +initial_currentden/prob_parm->cat_conductivity*(x-prob_lo[axdir]);
                }
            }
        }
    }
}
#endif
