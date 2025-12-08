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
    }

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
        amrex::Real dcoeff_inv=phi(outcell1,ETRODE_ID)/prob_parm.dc_etrode;
        dcoeff_inv += phi(outcell1,ELYTE_ID)/prob_parm.dc_elyte;
        amrex::Real dcoeff=1.0/dcoeff_inv;
        amrex::Real jbv=-dcoeff*dphidn;


        if(solved_comp==POT_ID)
        {
            if(conjsolve==0)
            {
                amrex::Real phi_etrode=0.5*(5.0*phi(outcell1,POT_ID)-3.0*phi(outcell2,POT_ID));
                amrex::Real jbv_etrode_to_elyte=-jbv; //BV current always directed from electrode to electrolyte
                amrex::Real phi_elyte=phi_etrode-bvinv(jbv_etrode_to_elyte);
                amrex::Real dc_inv=phi(intcell,ETRODE_ID)/prob_parm.dc_etrode;
                dc_inv += phi(intcell,ELYTE_ID)/prob_parm.dc_elyte;
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
