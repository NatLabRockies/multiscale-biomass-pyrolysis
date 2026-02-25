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
    
    /*AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real bv(amrex::Real phi_jump)
    {
        amrex::Real ocp=0.2;
        amrex::Real j0=3.0;
        amrex::Real phi0=1.0;

       Real jbv=j0*std::sinh((phi_jump-ocp)/phi0);
       return(jbv);
    }
    
    AMREX_GPU_DEVICE AMREX_INLINE 
    amrex::Real bvinv(amrex::Real jbv)
    {
       Real ocp=0.2;
       Real j0=3.0;
       Real phi0=1.0;
    
      Real phi_jump=std::asinh(jbv/j0)*phi0+ocp;
      return(phi_jump);
    }*/

    AMREX_GPU_DEVICE AMREX_INLINE
    void bc_ib(IntVect face,
               int dir, int sgn,
               int solved_comp,
               Array4<Real> const& phi,
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
        //sgn is outward normal with the usual conjsolve=0
        //solve
        IntVect lcell=face;
        IntVect rcell=face;
        lcell[dir]-=1;
        amrex::Real posvec[2];
        posvec[0]=prob_lo[0]+(lcell[0]+0.5)*dx[0];
        posvec[1]=prob_lo[1]+(lcell[1]+0.5)*dx[1];
        posvec[dir]+=0.5*dx[dir];

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

        //this is always the case no matter the sgn
        amrex::Real dphidn=(phi(outcell2,POT_ID)-phi(outcell1,POT_ID))/dx[dir];
        amrex::Real dcoeff_inv=phi(outcell1,ETRODE_ID)/prob_parm.cat_conductivity;
        dcoeff_inv += phi(outcell1,ELYTE_ID)/prob_parm.elyte_conductivity;
        amrex::Real dcoeff=1.0/dcoeff_inv;
        amrex::Real jbv=-dcoeff*dphidn;


        if(solved_comp==POT_ID)
        {
            if(conjsolve==0)
            {
                amrex::Real phi_etrode=0.5*(5.0*phi(outcell1,POT_ID)-3.0*phi(outcell2,POT_ID));
                
                //BV current always directed from electrode to electrolyte
                amrex::Real jbv_etrode_to_elyte=-jbv;
                amrex::Real phi_elyte=bvsolve_complex(-1.0,
                                    phi_etrode,prob_parm.co2conc,prob_parm.T,
                                    jbv_etrode_to_elyte,prob_parm);
               
                amrex::Print()<<"jbv etrde to elyte;"<<jbv_etrode_to_elyte<<"\n"; 
                amrex::Real dc_inv=phi(intcell,ETRODE_ID)/prob_parm.cat_conductivity;
                dc_inv += phi(intcell,ELYTE_ID)/prob_parm.elyte_conductivity;
                amrex::Real dc=1.0/dc_inv;

                //remember we are always solving dphidt+ del.flux=0
                //flux=-dc*dphidx
                
                // For sgn=-1 case (left cell is 0 and right cell is 1)
                // =======================================================
                // (out flux_right (handled by mlabec) + outward flux_left (handled here))/dx=0
                // flux_left= -dc (phi(intcell)-phi_dirichlet)/(0.5*dx) \hat{x}
                // outward flux left = flux_left . sgn*\hat{x} = -flux_left
                // outward flux left = dc*phi(intcell)/(0.5*dx) - dc*phi_dirichlet/(0.5*dx)
                // so add dc/(0.5*dx)/dx to acoeff and add dc*phi_dirichlet/(0.5*dx)/dx to bcoeff
                
                // For sgn=1 case (left cell is 1 and right cell is 0)
                // =======================================================
                // (out flux_right (handled here) + outward flux_left (handled by mlabec))/dx=0
                // flux_right= -dc (phi_dirichlet-phi(intcell))/(0.5*dx) \hat{x}
                // outward flux right = flux_right . sgn*\hat{x} = flux_right
                // outward flux right = dc*phi(intcell)/(0.5*dx) - dc*phi_dirichlet/(0.5*dx)
                // so add dc/(0.5*dx)/dx to acoeff and add dc*phi_dirichlet/(0.5*dx)/dx to bcoeff

                acoeff(intcell) += dc/(0.5*dx2);
                rhs(intcell) += dc*phi_elyte/(0.5*dx2);
            }
            else
            {
                // =======================================================
                // (out flux_right (handled by mlabec) + outward flux_left (handled here))/dx=0
                //remember we are always solving dphidt+ del.flux=0
                //jbv here is always outward facing flux
                //finite volume scheme says:
                //(outward flux inside boundary (mlabec handles this) + 
                //outward flux at this boundary(handled here))/dx=0
                //so rhs always gets the negative of the outward flux contribution
                // =======================================================
                rhs(intcell)   -= jbv/dx[dir];
            }
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
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn ==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        //default
        robin_a(ghost_cell) = 0.0;
        robin_b(ghost_cell) = 1.0;
        robin_f(ghost_cell) = 0.0;
        bc_arr(ghost_cell) = 0.0;

        if(sgn == -1) 
        { 
            // lo sides
            if(spec_id==POT_ID)
            {
                robin_a(ghost_cell) = 1.0;
                robin_b(ghost_cell) = 0.0;
                robin_f(ghost_cell) = prob_parm.phiL;
                bc_arr(ghost_cell) = prob_parm.phiL;
            }
        }
        else
        {
            if(spec_id==POT_ID)
            {
                robin_a(ghost_cell) = 1.0;
                robin_b(ghost_cell) = 0.0;
                robin_f(ghost_cell) = prob_parm.phiR;
                bc_arr(ghost_cell) = prob_parm.phiR;
            }
        } 

    }
}
#endif
