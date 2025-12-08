#ifndef _BOUNDARIES_H_
#define _BOUNDARIES_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>
#include <Transport.H>

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
           func+=prob_parm.Zi[TMAp_ID]*prob_parm.tmap_conc; 

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
           //Puvi's guess of 10 mV
           current_phid=-0.01;
           amrex::AllPrint()<<"Donnan potential not converging\n";
        }

        return(current_phid);
    }


   AMREX_GPU_DEVICE AMREX_FORCE_INLINE
   amrex::Real co_bv_current(amrex::Real voltdiff,amrex::Real co2conc,
                             ProbParm const& prob_parm)
   {
       amrex::Real c0=1.0;
       amrex::Real echem_temp=prob_parm.T;
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real eta=voltdiff-prob_parm.E0_co;
       amrex::Real i_co=prob_parm.i0_co*(co2conc/c0)*std::exp(-prob_parm.alpha_co*eta*FbyRT);
       return(i_co);
   
   } 
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE
   amrex::Real co_bv_current_der(amrex::Real voltdiff,amrex::Real co2conc,
                             ProbParm const& prob_parm)
   {
       amrex::Real c0=1.0;
       amrex::Real echem_temp=prob_parm.T;
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real eta=voltdiff-prob_parm.E0_co;
       amrex::Real i_coder=prob_parm.i0_co*(co2conc/c0)*std::exp(-prob_parm.alpha_co*eta*FbyRT)*(-prob_parm.alpha_co*FbyRT);
       return(i_coder);
   
   } 
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE
   amrex::Real co_bv_current_inv(amrex::Real i_co,amrex::Real co2conc,
                             ProbParm const& prob_parm)
   {
       amrex::Real c0=1.0;
       amrex::Real echem_temp=prob_parm.T;
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real eta=std::log((i_co/prob_parm.i0_co)*(c0/co2conc))/FbyRT
       /(-prob_parm.alpha_co);
       return(eta+prob_parm.E0_co);
   }
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE
   amrex::Real her_bv_current(amrex::Real voltdiff,
                             ProbParm const& prob_parm)
   {
       amrex::Real eta=voltdiff-prob_parm.E0_h2;
       amrex::Real echem_temp=prob_parm.T;
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real i_h2=prob_parm.i0_h2*std::exp(-prob_parm.alpha_h2*eta*FbyRT);
       return(i_h2);
   }
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE
   amrex::Real her_bv_current_der(amrex::Real voltdiff,
                             ProbParm const& prob_parm)
   {
       amrex::Real eta=voltdiff-prob_parm.E0_h2;
       amrex::Real echem_temp=prob_parm.T;
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real i_h2der=prob_parm.i0_h2*std::exp(-prob_parm.alpha_h2*eta*FbyRT)*
       (-prob_parm.alpha_h2*FbyRT);
       return(i_h2der);
   }
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE
   amrex::Real her_bv_current_inv(amrex::Real i_h2,
                             ProbParm const& prob_parm)
   {
       amrex::Real echem_temp=prob_parm.T;
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real eta=std::log(i_h2/prob_parm.i0_h2)/FbyRT/(-prob_parm.alpha_h2);
       return(eta+prob_parm.E0_h2);
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
       i_co=2.0*prob_parm.i0_co*(co2conc/c0)*std::sinh(prob_parm.alpha_co*eta_co*FbyRT);
       i_h2=i_co/fe*(1.0-fe);
       amrex::Real bvcurrent=i_co+i_h2;
       return(bvcurrent);
   }
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
   amrex::Real  bvfunc_simple1(amrex::Real phi_i,
                       amrex::Real phi_e,
                       amrex::Real co2conc,amrex::Real echem_temp,
                       ProbParm const& prob_parm,amrex::Real &i_co,
                       amrex::Real &i_h2)
   {
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real c0=1.0;
       amrex::Real fe=prob_parm.far_eff;
       amrex::Real eta_co=phi_e-phi_i-prob_parm.E0_co;
       amrex::Real alpha=prob_parm.alpha_co;
       i_co=prob_parm.i0_co*(co2conc/c0)*( std::exp((1.0-alpha)*eta_co*FbyRT) -
                                           std::exp(-alpha*eta_co*FbyRT) );
       i_h2=i_co/fe*(1.0-fe);
       amrex::Real bvcurrent=i_co+i_h2;
       return(bvcurrent);
   }
   
   AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
   amrex::Real bvfuncder_simple1(amrex::Real phi_i,
                       amrex::Real phi_e,
                       amrex::Real co2conc,amrex::Real echem_temp,
                       ProbParm const& prob_parm,amrex::Real &i_coder,
                       amrex::Real &i_h2der)
   {
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real c0=1.0;
       amrex::Real fe=prob_parm.far_eff;
       amrex::Real eta_co=phi_e-phi_i-prob_parm.E0_co;
       amrex::Real alpha=prob_parm.alpha_co;
       i_coder=-prob_parm.i0_co*(co2conc/c0)*( (1.0-alpha)*FbyRT*
                                             std::exp((1.0-alpha)*eta_co*FbyRT) -
                                           (-alpha*FbyRT)*std::exp(-alpha*eta_co*FbyRT) );
       i_h2der=i_coder/fe*(1.0-fe);
       amrex::Real bvcurrent_der=i_coder+i_h2der;
       return(bvcurrent_der);
   }

   AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
   amrex::Real  bvfunc_complex(amrex::Real phi_i,
                       amrex::Real phi_e,
                       amrex::Real co2conc,amrex::Real echem_temp,
                       ProbParm const& prob_parm,amrex::Real &i_co,
                       amrex::Real &i_h2)
   {
   
       i_h2=her_bv_current(phi_e-phi_i,prob_parm);
       i_co=co_bv_current(phi_e-phi_i,co2conc,prob_parm);
       amrex::Real bvcurrent=i_h2+i_co;
       return(bvcurrent);
   }
   
   AMREX_GPU_DEVICE AMREX_INLINE 
   amrex::Real bvfuncder_complex(amrex::Real phi_i,amrex::Real phi_e,
                         amrex::Real co2conc,
                         amrex::Real echem_temp,
                         ProbParm const& prob_parm)
   {
       amrex::Real bvcurrentder=co_bv_current_der(phi_e-phi_i,co2conc,prob_parm);
       bvcurrentder+=her_bv_current_der(phi_e-phi_i,prob_parm);
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
      amrex::Real eta=std::asinh(0.5*fe*j_input/prob_parm.i0_co*(c0/co2conc))/prob_parm.alpha_co/FbyRT;
      amrex::Real phi_c=phi_e-eta-prob_parm.E0_co;
      //amrex::Print()<<"phi_c:"<<phi_c<<"\t"<<eta<<"\t"<<co2conc<<"\t"<<j_input<<"\n";
      return(phi_c);
  }

   AMREX_GPU_DEVICE AMREX_INLINE amrex::Real bvsolve_asympt(amrex::Real co2conc,
                                                             amrex::Real echem_temp,
                                                             amrex::Real j_input,
                                                             ProbParm const& prob_parm)
   {
       //solve eta f*(1-alpha*eta*f)=j/i0*(c0/c)
       //ax^2+bx+c=0
       amrex::Real soln1,soln2;
       amrex::Real c0=1.0;
       amrex::Real eps=1e-50;
       amrex::Real f=FBYR/echem_temp;
       amrex::Real a=-prob_parm.alpha_co*f*f;
       amrex::Real b=f;
       amrex::Real c=-prob_parm.far_eff*j_input/prob_parm.i0_co*(c0/(co2conc+eps));
       amrex::Real D=b*b-4*a*c;
       amrex::Real eta=0.0;
       if(D>0)
       {
          soln1=(-b+std::sqrt(D))/(2.0*a);
          soln2=(-b-std::sqrt(D))/(2.0*a);
          eta=(soln1*soln2<0.0)?std::max(soln1,soln2):
          std::copysign(1.0,soln1)*std::min(amrex::Math::abs(soln1),
                                            amrex::Math::abs(soln2));
          amrex::Print()<<"eta:"<<eta<<"\n";
       }
       else
       {
           amrex::Print()<<"a,b,c,D:"<<a<<"\t"<<b<<"\t"<<c<<"\t"<<D<<"\n";
           amrex::Abort();
       }
       return(eta);
   }

   AMREX_GPU_DEVICE AMREX_INLINE amrex::Real bvsolve_simple1(amrex::Real phi_i_guess,
                                                             amrex::Real phi_e,
                                                             amrex::Real co2conc,
                                                             amrex::Real echem_temp,
                                                             amrex::Real j_input,
                                                             ProbParm const& prob_parm)
   {
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real fe=prob_parm.far_eff;
       amrex::Real c0=1.0;
       amrex::Real alpha=prob_parm.alpha_co;
       amrex::Real eta_guess=0.90*std::asinh(fe*0.5*j_input/prob_parm.i0_co*(c0/co2conc))/alpha/FbyRT;
       amrex::Real phi_c_guess=phi_e-eta_guess-prob_parm.E0_co;
       
       //amrex::Real eta_guess=bvsolve_asympt(co2conc,echem_temp,j_input,prob_parm);
       //amrex::Real phi_c_guess=phi_i_guess;
       
       amrex::Real ih2,ico,ih2der,icoder;
       
       //solve for ionic potential from bv current
       int n_newton_iter=30;
       amrex::Real tol=1e-7;
       int converged=0;
       int it;
       amrex::Real phi_c=phi_c_guess;
       amrex::Real jbv,jbvder,dphi;

       for(it=0;it<n_newton_iter;it++)
       {
           jbv=bvfunc_simple1(phi_c,phi_e,co2conc,echem_temp,prob_parm,ih2,ico);
           if(amrex::Math::abs(j_input-jbv)<tol)
           {
               converged=1;
               break;
           }
           else
           {
               jbvder=bvfuncder_simple1(phi_c,phi_e,co2conc,echem_temp,prob_parm,ih2der,icoder);
               dphi=(j_input-jbv)/jbvder;
               phi_c += dphi;
           }
       }

       if(!converged)
       {
           amrex::AllPrint()<<"it,phi_c,dphi,jbv,j_input,tol:"
           <<it<<"\t"<<phi_c<<"\t"<<dphi<<"\t"<<jbv<<"\t"<<j_input<<"\t"<<tol<<"\n";
           //matching open circuit potential (this is wrong!)
           phi_c=phi_e-prob_parm.E0_co;
           amrex::AllPrint()<<"bvsolve not converging\n";
       }

       /*else
         {
         amrex::Real j_bv=bvfunc_complex(phi_c,phi_e,co2conc,echem_temp,prob_parm,ih2,ico);
         amrex::Print()<<"=====================hari\n";
         amrex::Print()<<"converged, j_input, j_bv, phi_c, phi_r:"<<j_input<<"\t"<<j_bv<<"\t"<<phi_c<<"\t"<<phi_e<<"\n";
         amrex::Print()<<"=====================\n";
         }*/
       //amrex::Print()<<"phi_c:"<<phi_c<<"\t"<<eta<<"\t"<<co2conc<<"\t"<<j_input<<"\n";
       /*amrex::Real phic_simple=bvsolve_simple(phi_i_guess,
                                              phi_e,co2conc,echem_temp,
                                              j_input,prob_parm);

       if(amrex::Math::abs(phi_c-phic_simple)>tol*1e-3)
       {
           amrex::AllPrint()<<"it,phi_c,phic_simple,dphi,jbv,j_input,tol:"
           <<it<<"\t"<<phi_c<<"\t"<<phic_simple<<"\t"<<dphi<<"\t"<<jbv<<"\t"<<j_input<<"\t"<<tol<<"\n";
           //matching open circuit potential (this is wrong!)
           amrex::AllPrint()<<"bvsolve not matching up\n";
       }*/

       return(phi_c);
   }

   AMREX_GPU_DEVICE AMREX_INLINE amrex::Real bvsolve_complex(amrex::Real dphi_guess,
                                                             amrex::Real phi_e,
                                                             amrex::Real co2conc,
                                                             amrex::Real echem_temp,
                                                             amrex::Real j_in,
                                                             ProbParm const& prob_parm)
   {
       //cut off to avoid 0
       amrex::Real j_cutoff=1e-5;
       amrex::Real j_max=prob_parm.max_current_density;
       amrex::Real j_inp1=(amrex::Math::abs(j_in)<j_cutoff)?j_cutoff:j_in;
       amrex::Real j_input=(amrex::Math::abs(j_inp1)>j_max)?std::copysign(1.0,j_inp1)*j_max:j_inp1;
       //solve for ionic potential from bv current
       int n_newton_iter=40;
       amrex::Real vdiff_co_guess=co_bv_current_inv(j_input,co2conc,prob_parm);
       amrex::Real vdiff_h2_guess=her_bv_current_inv(j_input,prob_parm);
       
       amrex::Real tol=1e-4;
       int converged=0;
       amrex::Real vdiff=0.5*(vdiff_co_guess+vdiff_h2_guess);
       amrex::Real ih2,ico;
       amrex::Real jbv;
       //amrex::Print()<<"\n=======================\n";
       for(int it=0;it<n_newton_iter;it++)
       {
           amrex::Real phi_i=phi_e-vdiff;
           jbv=bvfunc_complex(phi_i,phi_e,co2conc,echem_temp,prob_parm,ih2,ico);
           //amrex::Print()<<"vdiff:"<<vdiff<<"\t"<<vdiff_co_guess<<"\t"<<vdiff_h2_guess
           //<<"\t"<<amrex::Math::abs(jbv-j_input)<<"\n";
           if(amrex::Math::abs(jbv-j_input)<tol)
           {
               converged=1;
               break;
           }
           else
           {
               amrex::Real jbvder=bvfuncder_complex(phi_i,phi_e,co2conc,echem_temp,prob_parm);
               amrex::Real newton_change=-(jbv-j_input)/jbvder;
               vdiff += newton_change;
           }
       }
       if(!converged)
       {
           amrex::AllPrint()<<"bvsolve not converging\n";
           amrex::AllPrint()<<"j_input,co2conc:"<<j_input<<"\t"<<co2conc<<"\t"<<echem_temp<<"\t"<<
           vdiff_co_guess<<"\t"<<vdiff_h2_guess<<"\t"<<vdiff<<"\t"<<amrex::Math::abs(jbv-j_input)<<"\t"<<phi_e<<"\n";
           amrex::Abort();
          amrex::Print()<<"=====================hari\n";
       }
       return(phi_e-vdiff);
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
       int n_iter=40;
       amrex::Real phigrad=phigrad_guess;
       amrex::Real specgrad[NUM_TR_SPECIES]={0.0};
       amrex::Real FbyRT=FBYR/echem_temp;
       amrex::Real RTbyF=1.0/FbyRT;
       amrex::Real tol=1e-3;
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
           amrex::AllPrint()<<"err1,err2,tol:"<<err1<<"\t"<<err2<<"\t"<<tol<<"\t"<<phigrad<<"\n";
           amrex::AllPrint()<<"imsolve not converging, using phigrad as 0\n";
           phigrad=0.0;
       }

       return(phigrad);
   }
   AMREX_GPU_DEVICE AMREX_INLINE
   int check_cut_cell(amrex::Real val)
   {
       amrex::Real tol=1e-15;
       int cutcell=(val*(1.0-val)>tol)?1:0;
       return(cutcell);
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
      int activate_surf_chem=prob_parm.activate_surf_chem;
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

      /*IntVect outcell1=outcell;
      int outcell1_adjust=(conjsolve==0)?sgn:-sgn;
      outcell1[dir]+=outcell1_adjust;
      IntVect outcell2=outcell1;
      outcell2[dir] += outcell1_adjust;*/
     
      //assuming boundaries match cartesian boundaries 
      IntVect outcell1=outcell;
      int outcell1_adjust=(conjsolve==0)?sgn:-sgn;
      IntVect outcell2=outcell1;
      outcell2[dir] += outcell1_adjust;

      amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
      outward_normal[dir]=(conjsolve==0)?sgn:-sgn;
      amrex::Real dx2=dx[dir]*dx[dir];
      
      amrex::Real ylen=prob_hi[1]-prob_lo[1];
      amrex::Real yloc=posvec[1];
      //int im_coverage=(yloc < prob_lo[1]+prob_parm.imfrac*ylen)?1:0;
      int im_coverage=prob_parm.activate_im_coverage;
      amrex::Real FbyRT=FBYR/prob_parm.T;
      amrex::Real RTbyF=1.0/FbyRT;
      
      amrex::Real specflux[NUM_TR_SPECIES]={0.0};
      amrex::Real phi_boundary=0.0;
      amrex::Real cat_current_den=0.0;
      amrex::Real ci[NUM_TR_SPECIES]={0.0};

      //under-resolved check
      int c_out1=int(sb_arr(outcell1,CMASK_ID));
      int c_out2=int(sb_arr(outcell2,CMASK_ID));

      int well_resolved=(conjsolve==0)?(!c_out1*!c_out2):(c_out1*c_out2);
      int under_resolved=!well_resolved;
      
     /*if(under_resolved)
     {
        amrex::AllPrint()<<"under resolved boundaries..\n";
     }*/

      if(conjsolve==0)
      {
          amrex::Real phi_catalyst=0.0;
          amrex::Real ci_donnan[NUM_TR_SPECIES];
          amrex::Real phi_e=(well_resolved)?
          0.5*(3.0*sb_arr(outcell1,PHI_ID)-sb_arr(outcell2,PHI_ID))
          :sb_arr(outcell1,PHI_ID);
          
          //dummy value so that solve doesnt nan
          amrex::Real co2conc=1.0; 
          amrex::Real phi_i_guess=0.0;

          for(int sp=0;sp<NUM_TR_SPECIES;sp++)
          {
              //add donnan potential here
              ci[sp]=sb_arr(intcell,sp);
          }

          if(activate_surf_chem)
          {
              //this is always the case no matter the sgn for outward flux
              amrex::Real dphidn=(well_resolved)?
              (sb_arr(outcell2,PHI_ID)-sb_arr(outcell1,PHI_ID))/dx[dir]:0.0;
              amrex::Real sig=prob_parm.cat_conductivity;
              amrex::Real j_solid=-sig*dphidn;

              //current directed from electrolyte to electrode
              //for electrochemical reduction
              amrex::Real jbv=amrex::Math::abs(j_solid); 
              phi_i_guess=sb_arr(intcell,PHI_ID);
              co2conc=sb_arr(intcell,CO2_ID);

              if(prob_parm.bvsimple)
              {
                  phi_catalyst=bvsolve_simple(phi_i_guess,
                                              phi_e,co2conc,prob_parm.T,
                                              jbv,prob_parm);
              }
              else
              {
                  //amrex::Print()<<"solving complex, co2conc, jbv:"<<co2conc<<"\t"<<jbv<<"\n";
                  phi_catalyst=bvsolve_complex(-1.0,
                                               phi_e,co2conc,prob_parm.T,
                                               jbv,prob_parm);
              }

              //j_co and j_h2 are magnitudes
              amrex::Real j_co,j_h2;

              if(prob_parm.bvsimple)
              {
                  bvfunc_simple(phi_catalyst,phi_e,co2conc,
                               prob_parm.T,prob_parm,j_co,j_h2);
              }
              else
              {
                  bvfunc_complex(phi_catalyst,phi_e,co2conc,
                                 prob_parm.T,prob_parm,j_co,j_h2);
              }

              //activate reverse reaction
              if((co2conc < prob_parm.co2star) && prob_parm.bvsimple)
              {
                  j_co*=-1;
                  j_h2*=-1;
              }

              //outward flux
              specflux[CO_ID]=-j_co/2.0/FARADCONST; 
              specflux[CO2_ID]=j_co/2.0/FARADCONST;
              specflux[OHm_ID]=-(j_co+j_h2)/FARADCONST;
              specflux[H2_ID]=-j_h2/FARADCONST;

              if(im_coverage)
              {
                  amrex::Real deltaphi_donnan=newton_solve_donnan_potential(ci,prob_parm);
                  //amrex::Print()<<"deltaphi_donnan:"<<deltaphi_donnan<<"\n";

                  for(int i=0;i<NUM_TR_SPECIES;i++)
                  {
                      ci_donnan[i]=ci[i]*donnan_conc_ratio(deltaphi_donnan, prob_parm.Zi[i],prob_parm.T);
                  }
                  amrex::Real phigrad_guess=amrex::Math::abs(prob_parm.catelecpot/(prob_hi[0]-prob_lo[0])); //guess
                  //amrex::Real phigrad=imsolve(phigrad_guess,phi_e,deltaphi_donnan,co2conc,
                  //                            ci_donnan,prob_parm.T,
                  //                            prob_parm.im_thickness,
                  //                            specflux,
                  //                            prob_parm);
                  amrex::Real phigrad=0.0;

                  phi_boundary=phi_catalyst+deltaphi_donnan-phigrad*prob_parm.im_thickness;
              }
              else
              {
                  phi_boundary=phi_catalyst;
              }
          }
          else
          {
              //use BV zero current
              amrex::Real jbv=0.0; 
              if(prob_parm.bvsimple)
              {
                  phi_catalyst=bvsolve_simple(phi_i_guess,
                                              phi_e,co2conc,prob_parm.T,
                                              jbv,prob_parm);
              }
              else
              {
                  phi_catalyst=bvsolve_complex(-1.0,
                                               phi_e,co2conc,prob_parm.T,
                                               jbv,prob_parm);
              }

              if(im_coverage)
              {
                  //Puvi's estimate
                  amrex::Real deltaphi_donnan=-0.01;
                  phi_boundary=phi_catalyst+deltaphi_donnan;
              }
              else
              {
                  phi_boundary=phi_catalyst;
              }
          }
          //amrex::Print()<<"phi_boundary:"<<phi_boundary<<"\n";
      }
      else //conjsolve=1
      {
          cat_current_den=0.0;
          amrex::Real diff_current_den=0.0;
          amrex::Real eps=1e-30;
          amrex::Real dphidn=(well_resolved)?
          (sb_arr(outcell2,PHI_ID)-sb_arr(outcell1,PHI_ID))/dx[dir]:0.0;
          amrex::Real conc_face[NUM_TR_SPECIES]={0.0};
          for(int sp=0;sp<NUM_TR_SPECIES;sp++)
          {
              conc_face[sp]=(well_resolved)?1.5*sb_arr(outcell1,sp)
              -0.5*sb_arr(outcell2,sp):sb_arr(outcell1,sp);

              amrex::Real grad_conc=well_resolved?
              (sb_arr(outcell2,sp)-sb_arr(outcell1,sp))/dx[dir]:0.0;

              diff_current_den += prob_parm.kD[sp]/(conc_face[sp]+eps)*grad_conc;
          }
          amrex::Real sig=tr_transport::get_elyte_conductivity(conc_face,prob_parm);
          cat_current_den=-sig*dphidn + diff_current_den;
      }

      if(solved_comp==PHI_ID)
      {
          if(conjsolve==0)
          {
              //potential solve has 1 bcoeff on the electrolyte side
              amrex::Real sig=tr_transport::get_elyte_conductivity(ci,prob_parm);
              acoeff(intcell) += sig*1.0/(0.5*dx2);
              rhs(intcell) += sig*phi_boundary/(0.5*dx2);
          }
          else
          {
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
              if(phi(intcell,CMASK_ID)<1.0)
              {
                  robin_a(ghost_cell) = 1.0;
                  robin_b(ghost_cell) = 0.0;
                  robin_f(ghost_cell) = 0.0;
              }
              else
              {
                  robin_a(ghost_cell) = 0.0;
                  robin_b(ghost_cell) = 1.0;
                  robin_f(ghost_cell) = 0.0;
              }
          }
      }
      else
      {
          // hi sides
          if(spec_id==PHI_ID)
          {
              if(phi(intcell,CMASK_ID)<1.0)
              {
                  robin_a(ghost_cell) = 1.0;
                  robin_b(ghost_cell) = 0.0;
                  robin_f(ghost_cell) = prob_parm.catelecpot;
              }
              else
              {
                  robin_a(ghost_cell) = 0.0;
                  robin_b(ghost_cell) = 1.0;
                  robin_f(ghost_cell) = 0.0;
              }
          }
          
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
