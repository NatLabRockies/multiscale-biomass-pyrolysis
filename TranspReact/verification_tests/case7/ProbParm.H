#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real phiL=-1.0;
    amrex::Real phiR=0.0;
    amrex::Real  x1=0.2;
    
    amrex::Real i0_co=10.0;
    amrex::Real i0_h2=10.0;
    amrex::Real E0_co=-0.889;
    amrex::Real E0_h2=-0.828;
    amrex::Real alpha_co=0.4;
    amrex::Real alpha_h2=0.35;
    amrex::Real co2conc=100.0;
    amrex::Real T=333.15;
    amrex::Real current_input=1000; //A/m2
    amrex::Real cat_conductivity=10.0; //S/m
    amrex::Real elyte_conductivity=1.0; //S/m
    amrex::Real max_current_density=1e4;

    int axial_dir=0;
    int left_etrode=1;
};
#endif
