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

        for(int sp=0;sp<NUM_TR_SPECIES;sp++)
        {
            ntr += prob_parm.dcoeff[sp]*prob_parm.Zi[sp]*gradci[sp];
            dtr += prob_parm.dcoeff[sp]*(prob_parm.Zi[sp]*prob_parm.Zi[sp])*ci[sp];
        }
        efield=(ntr/dtr)*(GASCONST*prob_parm.T/FARADCONST);
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
               ProbParm const& prob_parm)
    {
        IntVect lcell=face;
        IntVect rcell=face;
        lcell[dir]-=1;
        amrex::Real xloc[2];
        xloc[0]=prob_lo[0]+(lcell[0]+0.5)*dx[0];
        xloc[1]=prob_lo[1]+(lcell[1]+0.5)*dx[1];
        xloc[dir]+=0.5*dx[dir];

        amrex::Real midy=0.5*(prob_lo[1]+prob_hi[1]);

        //interior cell
        IntVect intcell=(sgn==1)?lcell:rcell;
        IntVect intcell1=intcell;
        intcell1[dir]-=sgn;

        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        amrex::Real dx2=dx[dir]*dx[dir];

        amrex::Real gradci[NUM_TR_SPECIES];
        amrex::Real ci[NUM_TR_SPECIES];
        amrex::Real ntr_amb,dtr_amb,efield_amb;

        amrex::Real FbyRT=FARADCONST/GASCONST/prob_parm.T;
        amrex::Real RTbyF=1.0/FbyRT;
        amrex::Real phi_e=prob_parm.catelecpot;
        amrex::Real eta_co=phi_e-sb_arr(intcell,PHI_ID)-prob_parm.E0_co;
        amrex::Real eta_h2=phi_e-sb_arr(intcell,PHI_ID)-prob_parm.E0_h2;

        //amrex::Real eta_co=phi_e-prob_parm.E0_co;
        //amrex::Real eta_h2=phi_e-prob_parm.E0_h2;
        amrex::Real c0=1.0; 
        amrex::Real i_co=prob_parm.i0_co*(sb_arr(intcell,CO2_ID)/c0)*
        std::exp(-prob_parm.alpha_co*eta_co*FbyRT);
        amrex::Real i_h2=prob_parm.i0_h2*
        std::exp(-prob_parm.alpha_h2*eta_h2*FbyRT);
        amrex::Real i_tot=i_co+i_h2;

        //inward is -ve
        amrex::Real j_co=-i_co/2.0/FARADCONST;
        //outward is +ve
        amrex::Real j_co2=-j_co;
        //inward is -ve
        amrex::Real j_h2=-i_h2/FARADCONST;
        //inward is -ve
        amrex::Real j_ohm=-(i_co+i_h2)/FARADCONST;

        for(int sp=0;sp<NUM_TR_SPECIES;sp++)
        {
            ci[sp]=sb_arr(intcell,sp);

            //will backward difference when intcell is left
            //and forward difference when intcell is right
            gradci[sp]=sgn*(sb_arr(intcell,sp)-sb_arr(intcell1,sp))/dx[dir];

            //gradci from liquid to catalyst is
            gradci[sp]*=sgn;
        }

        get_ambipolar_quantities(ci,gradci,prob_parm,ntr_amb,dtr_amb,efield_amb);

        if(solved_comp==PHI_ID)
        {
            //We are solving -del^2 (phi)=source
            //in transpreact
            //so (-del(phi)|i+1/2 + del(phi)|i-1/2)/dx = source
            //(E_{i+1/2}-E_{i-1/2})/dx=Source
            //sum (Efield.outward_normal)/dx=source

            //amrex::Real phi_dirc=prob_parm.catelecpot;
            //acoeff(intcell)+= 1.0/(0.5*dx2);
            //rhs(intcell) += phi_dirc/(0.5*dx2);
            //

            //outward here means liquid to catalyst
            amrex::Real gradphi_outward=RTbyF*(-i_tot/FARADCONST - ntr_amb)/dtr_amb;
            amrex::Real efield_outward=-gradphi_outward;
            rhs(intcell) -= efield_outward/dx[dir];
            amrex::Print()<<"efield:"<<efield_outward<<"\n";
        }
        else if(solved_comp==CO2_ID) 
        {
            //Buttler Vohlmer flux
            rhs(intcell)   -= j_co2/dx[dir];
        }
        else if(solved_comp==OHm_ID)
        {
            //Buttler Vohlmer flux
            //rhs(intcell)   -= j_ohm/dx[dir];
        }
        else if(solved_comp==CO_ID)
        {
            //Buttler Vohlmer flux
            rhs(intcell)   -= j_co/dx[dir];
        }
        else
        {
            //zero flux
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

        amrex::Real outflx;

        if(sgn == -1) 
        { 
            // lo sides
            if(spec_id==PHI_ID)
            {
                bc_arr(ghost_cell) = 0.0;
            }
            else if(spec_id==CO2_ID)
            {
                bc_arr(ghost_cell)=0.0;
            }
            else if(spec_id==CO_ID)
            {
                bc_arr(ghost_cell)=0.0;
            }
            else if(spec_id==TMAp_ID)
            {
                bc_arr(ghost_cell)=0.0;
            }
            else //all ions: H+,OH-,HCO3-,CO32-,HCOO-,K+
            {
                Real mob=prob_parm.dcoeff[spec_id]*prob_parm.Zi[spec_id]*FARADCONST/GASCONST/prob_parm.T;
                amrex::Real flx=mob*phi(cell_int,spec_id)*phi(cell_int,EFX_ID+dir);
                //amrex::Real flx_n=mob*phi(cell_int,spec_id)*phi(cell_int,EFX_ID+dir)*outward_normal[dir];
                //bc_arr(ghost_cell)=(flx_n>0.0)?(-flx/prob_parm.dcoeff[spec_id]):0.0;
                bc_arr(ghost_cell)=-flx/prob_parm.dcoeff[spec_id];
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
                //grad phi = sum_i (zi Di grad(ci))/sum_i (zi^2 Di ci)
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
