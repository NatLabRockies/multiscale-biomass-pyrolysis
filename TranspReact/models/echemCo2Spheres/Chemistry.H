#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>
#include <ProbParm.H>
#include <UnivConstants.H>
#include <Transport.H>

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
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            specsource(i,j,k,sp) = 0.0;
        }
        specsource(i,j,k,PHI_ID)=0.0;
        
        //Reactions
        amrex::Real kreac,reactrate;
        amrex::Real sourcefac=prob_parm.sourcefac;

        //H2O <-> Hp + OHm
        /*amrex::Real kwf=2.4e-2;//mol/m3/s
        amrex::Real kwb=2.4e6; //m3/mol/s
        reactrate=kwf-kwb*phi(i,j,k,Hp_ID)*phi(i,j,k,OHm_ID);
        specsource(i,j,k,Hp_ID)+=reactrate;
        specsource(i,j,k,OHm_ID)+=reactrate;*/
        
        //CO2 + H2O <-> HCO3m + Hp
        amrex::Real k2f=0.0016618;
        amrex::Real k2b=3.2872;
        reactrate=k2f*phi(i,j,k,CO2_ID)-k2b*phi(i,j,k,HCO3m_ID)*phi(i,j,k,Hp_ID);
        specsource(i,j,k,CO2_ID)-=reactrate;
        specsource(i,j,k,HCO3m_ID)+=reactrate;
        specsource(i,j,k,Hp_ID)+=reactrate;
        
        //CO32m + Hp <-> HCO3m
        amrex::Real k3f=5e4;
        amrex::Real k3b=0.065587;
        reactrate=k3f*phi(i,j,k,CO32m_ID)*phi(i,j,k,Hp_ID)-
        k3b*phi(i,j,k,HCO3m_ID);
        specsource(i,j,k,CO32m_ID)-=reactrate;
        specsource(i,j,k,Hp_ID)-=reactrate;
        specsource(i,j,k,HCO3m_ID)+=reactrate;
        
        //CO2 + OHm <-> HCO3m
        amrex::Real k4f=2.23;
        amrex::Real k4b=5.23e-5;
        reactrate=k4f*phi(i,j,k,CO2_ID)*phi(i,j,k,OHm_ID)-
        k4b*phi(i,j,k,HCO3m_ID);
        specsource(i,j,k,CO2_ID)-=reactrate;
        specsource(i,j,k,OHm_ID)-=reactrate;
        specsource(i,j,k,HCO3m_ID)+=reactrate;
        
        //HCO3m + OHm <-> CO32m + H2O
        /*amrex::Real k5f=6e6;
        amrex::Real k5b=1.07e6;
        reactrate=k5f*phi(i,j,k,HCO3m_ID)*phi(i,j,k,OHm_ID)-
        k5b*phi(i,j,k,CO32m_ID);
        specsource(i,j,k,HCO3m_ID)-=reactrate;
        specsource(i,j,k,OHm_ID)-=reactrate;
        specsource(i,j,k,CO32m_ID)+=reactrate;*/

        for(int sp=0;sp<NUM_TR_SPECIES;sp++)
        {
           specsource(i,j,k,sp)*=sourcefac;
        }
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
                
                int mask=int(s_arr(cellid,CMASK_ID));
                int mask_ip1=int(s_arr(ip1,CMASK_ID));
                int mask_im1=int(s_arr(im1,CMASK_ID));
                int normal_case=mask*mask_ip1*mask_im1;

                if(normal_case)
                {
                  state_arr(cellid,EFX_ID+dir)=-0.5*(s_arr(ip1,PHI_ID)-s_arr(im1,PHI_ID))/dx[dir];
                }
                else if(mask==1)
                {
                    //doing onesided first order
                    //can do better
                    if(mask_im1==1)
                    {
                        state_arr(cellid,EFX_ID+dir)=-(s_arr(cellid,PHI_ID)-s_arr(im1,PHI_ID))/dx[dir];
                    }
                    else if(mask_ip1==1)
                    {
                        state_arr(cellid,EFX_ID+dir)=-(s_arr(ip1,PHI_ID)-s_arr(cellid,PHI_ID))/dx[dir];
                    }
                    else
                    {
                        //not well resolved case
                        //1 valid cell between two dead cells
                        state_arr(cellid,EFX_ID+dir)=0.0;
                    }
                }
                else
                {
                    state_arr(cellid,EFX_ID+dir)=0.0;
                }
            }
            else if(cellid[dir]==domlo[dir])
            {
                int mask=int(s_arr(cellid,CMASK_ID));
                if(mask==1)
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
                    } 
                    else 
                    {
                        state_arr(cellid,EFX_ID+dir)=-(s_arr(ip1,PHI_ID) - s_arr(cellid,PHI_ID))/dx[dir];
                    }
                }
            }
            else if(cellid[dir]==domhi[dir])
            {
                int mask=int(s_arr(cellid,CMASK_ID));
                if(mask==1)
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
        IntVect lcell=faceid;
        IntVect rcell=faceid;
        lcell[dir]-=1;
        amrex::Real conc_face[NUM_TR_SPECIES],grad_conc[NUM_TR_SPECIES];
        amrex::Real efield;
        amrex::Real eps=1e-30;
        
        //flux for species is t+/F j
        //t+/F (-sig grad(phi) + kD grad(ln(c)))
        //t+/F (-sig grad(phi) + kD/c grad(c)))

        int left_phys_boundary=(faceid[dir]==domlo[dir]);
        int right_phys_boundary=(faceid[dir]==(domhi[dir]+1));
        int internal_face=(!left_phys_boundary && !right_phys_boundary);

        int mask_L=int(sb_arr(lcell,CMASK_ID));
        int mask_R=int(sb_arr(rcell,CMASK_ID));

        //1 when both mask_L and mask_R are 0
        int covered_interface=(!mask_L)*(!mask_R);
        //1 when both mask_L and mask_R are 1
        int regular_interface=(mask_L)*(mask_R);
        //1-0 or 0-1 interface
        int cov_uncov_interface=(mask_L)*(!mask_R)+(!mask_L)*(mask_R);


        if(!covered_interface && conjsolve==0 && spec_id<=PHI_ID)
        {
            if(internal_face && regular_interface)
            {
                for(int sp=0;sp<NUM_TR_SPECIES;sp++)
                {
                    conc_face[sp]=0.5*(sb_arr(rcell,sp)+sb_arr(lcell,sp));
                    grad_conc[sp]=(sb_arr(rcell,sp)-sb_arr(lcell,sp))/dx[dir];
                }
                efield=-(sb_arr(rcell,PHI_ID)-sb_arr(lcell,PHI_ID))/dx[dir];
            }
            else
            {
                IntVect intcell=faceid;
                IntVect intcell1=faceid;
                int sgn=0;
                if(cov_uncov_interface)
                {
                    sgn=(mask_L==1)?1:-1;
                    intcell=(sgn==1)?lcell:rcell;
                    intcell1=intcell;
                    intcell1[dir]-=sgn;
                }
                else
                {
                    sgn=(right_phys_boundary==1)?1:-1;
                    intcell=(sgn==1)?lcell:rcell;
                    intcell1=intcell;
                    intcell1[dir]-=sgn;
                }

                for(int sp=0;sp<NUM_TR_SPECIES;sp++)
                {
                    conc_face[sp]=1.5*sb_arr(intcell,sp)-0.5*sb_arr(intcell1,sp);
                    grad_conc[sp]=(conc_face[sp]-sb_arr(intcell,sp))/(0.5*dx[dir])*sgn;
                }
                amrex::Real phi_face=1.5*sb_arr(intcell,PHI_ID)-0.5*sb_arr(intcell1,PHI_ID);
                amrex::Real gradphi=(phi_face-sb_arr(intcell,PHI_ID))/(0.5*dx[dir])*sgn;
                efield=-gradphi;
            }

            if(spec_id<PHI_ID)
            {
                flx_arr(faceid) += tr_transport::get_elyte_conductivity(conc_face,prob_parm)
                *efield*prob_parm.tplus[spec_id]/FARADCONST;
                
                for(int sp=0;sp<NUM_TR_SPECIES;sp++)
                {
                    flx_arr(faceid) += prob_parm.kD[sp]/(conc_face[sp]+eps)*grad_conc[sp]*prob_parm.tplus[spec_id]/FARADCONST; 
                }
            }
            else
            {
                for(int sp=0;sp<NUM_TR_SPECIES;sp++)
                {
                    flx_arr(faceid) += prob_parm.kD[sp]/(conc_face[sp]+eps)*grad_conc[sp]; 
                }
            }
        }
    }
}
#endif
