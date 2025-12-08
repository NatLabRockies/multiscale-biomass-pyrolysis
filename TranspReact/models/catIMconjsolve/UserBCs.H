#ifndef _BOUNDARIES_H_
#define _BOUNDARIES_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>

using namespace amrex;
namespace tr_boundaries
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void get_ambipolar_quantities(amrex::Real ci[NUM_TR_SPECIES],
                                  amrex::Real gradci[NUM_TR_SPECIES],
                                  ProbParm const& prob_parm,
                                  amrex::Real &ntr,
                                  amrex::Real &dtr,
                                  amrex::Real &efield)
    {
        ntr=0.0;
        dtr=0.0;
        amrex::Real RTbyF=(GASCONST*prob_parm.T/FARADCONST);

        for(int sp=0;sp<NUM_TR_SPECIES;sp++)
        {
            ntr += prob_parm.dcoeff[sp]*prob_parm.Zi[sp]*gradci[sp];
            dtr += prob_parm.dcoeff[sp]*(prob_parm.Zi[sp]*prob_parm.Zi[sp])*ci[sp];
        }
        efield=(ntr/dtr)*RTbyF;
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real donnan_conc_ratio(amrex::Real dPhi, amrex::Real Zi,amrex::Real T)
    {
        amrex::Real FbyRT=FARADCONST/(GASCONST*T);
        amrex::Real donnan_ratio = std::exp(dPhi*Zi*FbyRT);
        return(donnan_ratio);
    }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real donnan_conc_ratio_der(amrex::Real dPhi, amrex::Real Zi,amrex::Real T)
    {
        amrex::Real FbyRT=FARADCONST/(GASCONST*T);
        amrex::Real donnan_ratio = std::exp(dPhi*Zi*FbyRT);
        amrex::Real donnan_ratio_der=Zi*FbyRT*donnan_ratio;
        return(donnan_ratio_der);
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real newton_solve_donnan_potential(amrex::Real Ci[NUM_TR_SPECIES],ProbParm const& prob_parm)
    {
        //equation in sigma_ zi ci exp(phid* ZiF/RT) + z_tmap c_tmap=0
        amrex::Real phid_guess=-0.01; //50 mV

        amrex::Real current_phid=phid_guess;
        int n_iters=150;
        amrex::Real tol=1e-5;
        amrex::Real func,funcder;
        int converged=0;

        for(int i=0;i<n_iters;i++)
        {
           func=0.0;
           funcder=0.0;
           for(int sp=0;sp<NUM_TR_SPECIES;sp++)
           {
                if(sp!=TMAp_ID)
                {
                   func += prob_parm.Zi[sp]*Ci[sp]
                   *donnan_conc_ratio(current_phid,prob_parm.Zi[sp],prob_parm.T);
                   funcder += prob_parm.Zi[sp]*Ci[sp]
                   *donnan_conc_ratio_der(current_phid,prob_parm.Zi[sp],prob_parm.T);
                }
           }
           //func+=prob_parm.Zi[TMAp_ID]*prob_parm.tmap_conc; 
           func+=prob_parm.Zi[TMAp_ID]*Ci[TMAp_ID]; 

           if(amrex::Math::abs(func) < tol)
           {
                //amrex::Print()<<"i,func:"<<i<<"\t"<<func<<"\tdonnan converged\n";
                converged=1;
                break;
           }
           else
           {
              current_phid += (-func/funcder);
           }
        }

        if(!converged)
        {
           amrex::AllPrint()<<"func:"<<func<<"\t"<<current_phid<<"\t"<<tol<<"\n";
           //current_phid=phid_guess;
           amrex::Abort("Donnan potential not converging");
        }

        return(current_phid);
    }
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
   amrex::Real  bvfunc_simple(amrex::Real phi_i,
                       amrex::Real phi_e,
                       amrex::Real co2conc,amrex::Real echem_temp,
                       ProbParm const& prob_parm,amrex::Real &i_co,
                       amrex::Real &i_h2)
   {
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real c0=1.0;
       amrex::Real fe=prob_parm.far_eff;
       amrex::Real eta_co=phi_e-phi_i-prob_parm.E0_co;
       i_co=prob_parm.i0_co*(co2conc/c0)*std::sinh(prob_parm.alpha_co*eta_co*FbyRT);
       i_h2=i_co/fe*(1.0-fe);
       amrex::Real bvcurrent=i_co+i_h2;
       return(bvcurrent);
   }

   AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
   amrex::Real  bvfunc_complex(amrex::Real phi_i,
                       amrex::Real phi_e,
                       amrex::Real co2conc,amrex::Real echem_temp,
                       ProbParm const& prob_parm,amrex::Real &i_co,
                       amrex::Real &i_h2)
   {
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real c0=1.0;
       amrex::Real eta_co=phi_e-phi_i-prob_parm.E0_co;
       amrex::Real eta_h2=phi_e-phi_i-prob_parm.E0_h2;
       i_co=prob_parm.i0_co*(co2conc/c0)*std::exp(-prob_parm.alpha_co*eta_co*FbyRT);
       i_h2=prob_parm.i0_h2*std::exp(-prob_parm.alpha_h2*eta_h2*FbyRT);
       amrex::Real bvcurrent=i_co+i_h2;
       return(bvcurrent);
   }

   AMREX_GPU_DEVICE AMREX_INLINE 
   amrex::Real bvfuncder_complex(amrex::Real phi_i,amrex::Real phi_e,
                         amrex::Real co2conc,
                         amrex::Real echem_temp,
                         ProbParm const& prob_parm)
   {
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real c0=1.0;
       amrex::Real eta_co=phi_e-phi_i-prob_parm.E0_co;
       amrex::Real eta_h2=phi_e-phi_i-prob_parm.E0_h2;
       amrex::Real i_coder=prob_parm.i0_co*(co2conc/c0)*std::exp(-prob_parm.alpha_co*eta_co*FbyRT)*(prob_parm.alpha_co*FbyRT);
       amrex::Real i_h2der=prob_parm.i0_h2*std::exp(-prob_parm.alpha_h2*eta_h2*FbyRT)*(prob_parm.alpha_h2*FbyRT);
       amrex::Real bvcurrentder=i_coder+i_h2der;
       return(bvcurrentder);
   }
  
   AMREX_GPU_DEVICE AMREX_INLINE amrex::Real bvsolve_simple(amrex::Real phi_i_guess,
                                       amrex::Real phi_e,
                                       amrex::Real co2conc,
                                       amrex::Real echem_temp,
                                       amrex::Real j_input,
                                       ProbParm const& prob_parm)
  {
    
      amrex::Real FbyRT=FBYR/echem_temp;
      amrex::Real fe=prob_parm.far_eff;
      amrex::Real c0=1.0;
      amrex::Real eta=std::asinh(fe*j_input/prob_parm.i0_co*(c0/co2conc))/prob_parm.alpha_co/FbyRT;
      amrex::Real phi_c=phi_e-eta-prob_parm.E0_co;
      //amrex::Print()<<"phi_c:"<<phi_c<<"\t"<<eta<<"\t"<<co2conc<<"\t"<<j_input<<"\n";
      return(phi_c);
  }

  AMREX_GPU_DEVICE AMREX_INLINE amrex::Real bvsolve_complex(amrex::Real phi_i_guess,
                                       amrex::Real phi_e,
                                       amrex::Real co2conc,
                                       amrex::Real echem_temp,
                                       amrex::Real j_input,
                                       ProbParm const& prob_parm)
  {
    //solve for ionic potential from bv current
    int n_newton_iter=30;
    amrex::Real phi_c=phi_i_guess;
    amrex::Real ih2,ico,jbv,jbvder,dphi;
    amrex::Real tol=1e-5;
    int converged=0;
    for(int it=0;it<n_newton_iter;it++)
    {
        amrex::Real jbv=bvfunc_complex(phi_c,phi_e,co2conc,echem_temp,prob_parm,ih2,ico);
        //print("jbv,j_input",jbv,j_input,phi_c)
        //amrex::Print()<<"jbv,j_input"<<jbv<<"\t"<<j_input<<"\t"<<phi_i_guess<<"\t"<<phi_e<<"\n";
        if(amrex::Math::abs(j_input-jbv)<tol)
        {
            converged=1;
            break;
        }
        else
        {
            jbvder=bvfuncder_complex(phi_c,phi_e,co2conc,echem_temp,prob_parm);
            dphi=(j_input-jbv)/jbvder;
            //print("dphi:",dphi)
            phi_c += dphi;
        }
     }
     if(!converged)
     {
        amrex::AllPrint()<<"phi_c,dphi,jbv,j_input,tol:"
        <<phi_c<<"\t"<<dphi<<"\t"<<jbv<<"\t"<<j_input<<"\t"<<tol<<"\n";
        //phi_c=phi_i_guess;
        amrex::Abort("bvsolve not converging");
     }

     return(phi_c);
  }

  AMREX_GPU_DEVICE
  AMREX_INLINE amrex::Real imsolve(amrex::Real phigrad_guess,
                            amrex::Real phi_e,
                            amrex::Real phi_donnan,
                            amrex::Real co2conc,
                            amrex::Real conc_im_donnan[NUM_TR_SPECIES],
                            amrex::Real echem_temp,
                            amrex::Real im_thickness,
                            amrex::Real specflux[NUM_TR_SPECIES],
                            ProbParm const& prob_parm)
  {
    int n_iter=20;
    amrex::Real phigrad=phigrad_guess;
    amrex::Real specgrad[NUM_TR_SPECIES]={0.0};
    amrex::Real FbyRT=FBYR/echem_temp;
    amrex::Real RTbyF=1.0/FbyRT;
    amrex::Real tol=1e-6;
    int converged=0;
    amrex::Real err1=0.0;
    amrex::Real err2=0.0;

    for(int it=0;it<n_iter;it++)
    {
        for(int i=0;i<NUM_TR_SPECIES;i++)
        {
            amrex::Real lhs=-specflux[i]/prob_parm.dcoeff[i]
            -prob_parm.Zi[i]*FbyRT*conc_im_donnan[i]*phigrad;
            amrex::Real Afac=(prob_parm.Zi[i]*FbyRT*im_thickness*phigrad+1.0);
            specgrad[i]=lhs/Afac;
        }

        amrex::Real sum_fiZi_by_Di=0.0;
        amrex::Real sum_Zi2_ci=0.0;
        for(int i=0;i<NUM_TR_SPECIES;i++)
        { 
            sum_fiZi_by_Di += specflux[i]*prob_parm.Zi[i]
            /prob_parm.dcoeff[i];
            sum_Zi2_ci += (prob_parm.Zi[i]*prob_parm.Zi[i])*(specgrad[i]*im_thickness+conc_im_donnan[i]);
        }

        phigrad=-RTbyF*(sum_fiZi_by_Di)/sum_Zi2_ci;
        err1=0.0;
        err2=0.0;
        for(int i=0;i<NUM_TR_SPECIES;i++)
        {
            err1 += prob_parm.Zi[i]*specgrad[i];
            amrex::Real flux_actual=specflux[i];
            amrex::Real flux=-prob_parm.dcoeff[i]*(specgrad[i]+prob_parm.Zi[i]*FbyRT
                                                   *(specgrad[i]*im_thickness+conc_im_donnan[i])*phigrad);
            err2 += (flux-flux_actual)*(flux-flux_actual);
        }
        err2=std::sqrt(err2);
        if(amrex::Math::abs(err1) < tol && err2 < tol)
        {
            converged=1;
            break;
        }
    }
    if(!converged)
    {
        amrex::AllPrint()<<"err1,err2,tol:"<<err1<<"\t"<<err2<<"\t"<<tol<<"\n";
        amrex::Abort("imsolve not converging");
    }

    return(phigrad);
  }

  AMREX_GPU_DEVICE AMREX_INLINE
  void bc_ib(IntVect face,
             int dir, int sgn,
             int solved_comp,
             Array4<Real> const& sb_arr,
             Array4<Real> const& acoeff,
             Array4<Real> const& rhs,
             GpuArray<int, AMREX_SPACEDIM> domlo,
             GpuArray<int, AMREX_SPACEDIM> domhi,
             GpuArray<Real, AMREX_SPACEDIM> prob_lo,
             GpuArray<Real, AMREX_SPACEDIM> prob_hi,
             GpuArray<Real, AMREX_SPACEDIM> dx,
             const Real time,
             ProbParm const& prob_parm,int conjsolve)
  {
      
      IntVect lcell=face;
      IntVect rcell=face;
      lcell[dir]-=1;
      amrex::Real posvec[2];
      posvec[0]=prob_lo[0]+(lcell[0]+0.5)*dx[0];
      posvec[1]=prob_lo[1]+(lcell[1]+0.5)*dx[1];
      posvec[dir]+=0.5*dx[dir];
      amrex::Real midx=0.5*(prob_lo[0]+prob_hi[0]);
      amrex::Real midy=0.5*(prob_lo[1]+prob_hi[1]);

      //interior cell
      IntVect intcell;
      IntVect outcell;

      if(conjsolve==0)
      {
          intcell=(sgn==1)?lcell:rcell;
          outcell=(sgn==1)?rcell:lcell;
      }
      else
      {
          intcell=(sgn==1)?rcell:lcell;
          outcell=(sgn==1)?lcell:rcell;
      }

      IntVect intcell1=intcell;
      int intcell1_adjust=(conjsolve==0)?-sgn:sgn;
      intcell1[dir]+=intcell1_adjust;

      IntVect outcell1=outcell;
      int outcell1_adjust=(conjsolve==0)?sgn:-sgn;
      outcell1[dir]+=outcell1_adjust;
      IntVect outcell2=outcell1;
      outcell2 += outcell1_adjust;

      amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
      outward_normal[dir]=(conjsolve==0)?sgn:-sgn;
      amrex::Real dx2=dx[dir]*dx[dir];
      
      amrex::Real specflux[NUM_TR_SPECIES]={0.0};
      amrex::Real phi_boundary=0.0;
      amrex::Real cat_current_den=0.0;
      
      amrex::Real ylen=prob_hi[1]-prob_lo[1];
      amrex::Real yloc=posvec[1];
      int im_coverage=(yloc < prob_lo[1]+prob_parm.imfrac*ylen)?1:0;
      //int im_coverage=0;
      amrex::Real FbyRT=FBYR/prob_parm.T;
      amrex::Real RTbyF=1.0/FbyRT;

      if(conjsolve==0)
      {

          //this is always the case no matter the sgn for outward flux
          amrex::Real dphidn=(sb_arr(outcell2,PHI_ID)-sb_arr(outcell1,PHI_ID))/dx[dir];
          amrex::Real sig=prob_parm.cat_conductivity;
          amrex::Real jbv=-sig*dphidn;

          amrex::Real ci[NUM_TR_SPECIES];
          amrex::Real ci_donnan[NUM_TR_SPECIES];

          amrex::Real phi_e=0.5*(5.0*sb_arr(outcell1,PHI_ID)-3.0*sb_arr(outcell2,PHI_ID));
          amrex::Real deltaphi_donnan;

          //current directed from catalyst to electrolyte
          amrex::Real j_input=amrex::Math::abs(jbv); 
          amrex::Real phi_i_guess=sb_arr(intcell,PHI_ID);
          amrex::Real co2conc=sb_arr(intcell,CO2_ID);

          for(int sp=0;sp<NUM_TR_SPECIES;sp++)
          {
              //add donnan potential here
              ci[sp]=sb_arr(intcell,sp);
          }
          ci[TMAp_ID]=1.5*ci[OHm_ID];

          //amrex::Print()<<"jbv,i,j,k:"<<jbv<<"\t"<<intcell<<"\n";

          //amrex::Real phi_catalyst=bvsolve(phi_i_guess,
          //                                 phi_e,co2conc,prob_parm.T,
          //                                 j_input,prob_parm);
          
          amrex::Real phi_catalyst=bvsolve_simple(phi_i_guess,
                                           phi_e,co2conc,prob_parm.T,
                                           j_input,prob_parm);

          //j_co and j_h2 are magnitudes
          amrex::Real j_co,j_h2;
          //bvfunc(phi_catalyst,phi_e,co2conc,
          //       prob_parm.T,prob_parm,j_co,j_h2);
          
          bvfunc_simple(phi_catalyst,phi_e,co2conc,
                 prob_parm.T,prob_parm,j_co,j_h2);


          //outward flux
          specflux[CO_ID]=-j_co/2.0/FARADCONST; 
          specflux[CO2_ID]=j_co/2.0/FARADCONST;
          specflux[OHm_ID]=-(j_co+j_h2)/FARADCONST;

          if(im_coverage)
          {
              deltaphi_donnan=newton_solve_donnan_potential(ci,prob_parm);
              //amrex::Print()<<"deltaphi_donnan:"<<deltaphi_donnan<<"\n";

              for(int i=0;i<NUM_TR_SPECIES;i++)
              {
                  ci_donnan[i]=ci[i]*donnan_conc_ratio(deltaphi_donnan, prob_parm.Zi[i],prob_parm.T);
              }
              amrex::Real phigrad_guess=20.0; //guess
              amrex::Real phigrad=imsolve(phigrad_guess,phi_e,deltaphi_donnan,co2conc,
                                          ci_donnan,prob_parm.T,
                                          prob_parm.im_thickness,
                                          specflux,
                                          prob_parm);

              phi_boundary=phi_catalyst+deltaphi_donnan-phigrad*prob_parm.im_thickness;
          }
          else
          {
              phi_boundary=phi_catalyst;
          }
          //amrex::Print()<<"phi_boundary:"<<phi_boundary<<"\n";
      }
      else
      {
          cat_current_den=0.0;
          for(int sp=0;sp<NUM_TR_SPECIES;sp++)
          {
              amrex::Real dphidn=(sb_arr(outcell2,PHI_ID)-sb_arr(outcell1,PHI_ID))/dx[dir];
              if(amrex::Math::abs(prob_parm.Zi[sp]) > 0.0)
              {
                  amrex::Real ci=sb_arr(outcell1,sp);
                  //amrex::Real phi_e=0.5*(5.0*sb_arr(outcell1,PHI_ID)-3.0*sb_arr(outcell2,PHI_ID));
                  amrex::Real dcidn=(sb_arr(outcell2,sp)-sb_arr(outcell1,sp))/dx[dir];
                  amrex::Real dcoeff=prob_parm.dcoeff[sp];
                  amrex::Real flux=-dcoeff*(dcidn + prob_parm.Zi[sp]*FbyRT*ci*dphidn);
                  cat_current_den += flux*prob_parm.Zi[sp]*FARADCONST;
                  //amrex::Print()<<"flux,current:"<<flux<<"\t"<<dcidn<<"\t"<<dphidn<<"\t"<<cat_current_den<<"\n";
              }
          }
          //cat_current_den*=-1;
          //amrex::Print()<<"cat_current_den:"<<cat_current_den<<"\n";
      }


      if(solved_comp==PHI_ID)
      {
          if(conjsolve==0)
          {
              //potential solve has 1 bcoeff on the electrolyte side
              acoeff(intcell) += 1.0/(0.5*dx2);
              rhs(intcell) += phi_boundary/(0.5*dx2);
          }
          else
          {
              //acoeff(intcell) += 1.0/(0.5*dx2);
              //rhs(intcell) += -0.17/(0.5*dx2);
              rhs(intcell) -= cat_current_den/dx[dir];
          }
      }
      else
      {
          rhs(intcell)   -= specflux[solved_comp]/dx[dir];
      }

  }

  AMREX_GPU_DEVICE AMREX_INLINE
  void species_bc(int i,int j,int k, int dir, int sgn, 
                  int spec_id, Array4<Real> const &phi, 
                  Array4<Real> const& bc_arr,
                  Array4<Real> const& robin_a,
                  Array4<Real> const& robin_b,
                  Array4<Real> const& robin_f,
                  GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                  GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                  GpuArray<Real, AMREX_SPACEDIM> dx,
                  const Real time,
                  ProbParm const& prob_parm)

  {
      //note: bdryLo/bdryHi grabs the face indices from bx that are at the boundary
      //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
      //so the ghost cell index at left side is i-1 while it is i on the right
      IntVect intcell{AMREX_D_DECL(i,j,k)};
      IntVect intcell1{AMREX_D_DECL(i,j,k)};
      IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
      amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
      outward_normal[dir]=sgn;
      int gcell_adjust,cell_adjust;

      //ghost_cell is one behind
      gcell_adjust= (sgn==-1)?-1:0;
      //cell is one behind
      cell_adjust = (sgn==1)?-1:0;

      ghost_cell[dir]+=gcell_adjust;
      intcell[dir]+=cell_adjust;
      //+1 for left side and -1 for right side
      intcell1[dir]=intcell[dir]-sgn;

      //default
      robin_a(ghost_cell) = 0.0;
      robin_b(ghost_cell) = 1.0;
      robin_f(ghost_cell) = 0.0;
      bc_arr(ghost_cell) = 0.0;

      if(sgn == -1) 
      { 
          // lo sides
          if(spec_id==PHI_ID)
          {
              robin_a(ghost_cell) = 1.0;
              robin_b(ghost_cell) = 0.0;
              robin_f(ghost_cell) = prob_parm.catelecpot;
              bc_arr(ghost_cell) = prob_parm.catelecpot;
          }
       }
       else
       {
           if(spec_id==CO2_ID)
           {
               bc_arr(ghost_cell) = prob_parm.co2_boundary;
           }
           else if(spec_id==OHm_ID)
           {
               //ideally a flux
               //j_ohm=m_koh*(c0-c_ohm); //c0=1000 mol/m3
               bc_arr(ghost_cell) = prob_parm.koh_init;
           }
           else if(spec_id==Kp_ID)
           {
               //ideally a flux
               //j_k+=m_koh*(c0-c_k+); //c0=1000 mol/m3
               bc_arr(ghost_cell) = prob_parm.koh_init;
           }
           else if(spec_id==PHI_ID)
           {
               bc_arr(ghost_cell) = 0.0;
           }
           else
           {
               bc_arr(ghost_cell) = 0.0;
           }
       } 
   }
}
#endif
