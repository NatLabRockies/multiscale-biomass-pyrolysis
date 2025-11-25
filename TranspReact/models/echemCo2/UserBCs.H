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

        amrex::Real outflx;

        if(sgn == -1) 
        { 
            // lo sides
            if(spec_id==PHI_ID)
            {
                bc_arr(ghost_cell) = prob_parm.catpot;
            }
            else if(spec_id==CO2_ID)
            {
                outflx=1.0e-3*phi(cell_int,CO2_ID)*outward_normal[dir]; //mol/m2/s
                bc_arr(ghost_cell)=-outflx/prob_parm.dcoeff[CO2_ID];
            }
            else if(spec_id==CO_ID)
            {
                outflx=-0.5e-3*phi(cell_int,CO2_ID)*outward_normal[dir]; //mol/m2/s
                bc_arr(ghost_cell)=-outflx/prob_parm.dcoeff[CO_ID];
            }
            else if(spec_id==HCOOm_ID)
            {
                outflx=-0.5e-3*phi(cell_int,CO2_ID)*outward_normal[dir]*0.0; //mol/m2/s
                bc_arr(ghost_cell)=-outflx/prob_parm.dcoeff[HCOOm_ID];
            }
            else if(spec_id==OHm_ID)
            {
                outflx=-(1.0e-3*2*phi(cell_int,CO2_ID)+0.5e-3*phi(cell_int,CO2_ID))*outward_normal[dir];
                bc_arr(ghost_cell)=-outflx/prob_parm.dcoeff[OHm_ID];
            }
            else
            {
                //Real mob=prob_parm.dcoeff[spec_id]*prob_parm.Zi[spec_id]*FARADCONST/GASCONST/prob_parm.T;
                //amrex::Real flx=mob*phi(cell_int,spec_id)*phi(cell_int,EFX_ID+dir);
                //amrex::Real flx_n=mob*phi(cell_int,spec_id)*phi(cell_int,EFX_ID+dir)*outward_normal[dir];
                //bc_arr(ghost_cell)=(flx_n>0.0)?(-flx/prob_parm.dcoeff[spec_id]):0.0;
                bc_arr(ghost_cell)=0.0;
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
                bc_arr(ghost_cell) = prob_parm.koh_init;
            }
            else if(spec_id==Kp_ID)
            {
                bc_arr(ghost_cell) = prob_parm.koh_init;
            }
            // lo sides
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
