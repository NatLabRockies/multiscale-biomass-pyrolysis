#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real dc_etrode = 3.0;
    amrex::Real dc_elyte = 1.0;
    amrex::Real  x1=0.2;
    amrex::Real phiL=1.0;
    amrex::Real phiR=0.0;
    int axial_dir=0;
    int left_etrode=1;
};
#endif
