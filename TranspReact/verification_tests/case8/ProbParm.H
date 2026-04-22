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
    amrex::Real phiL=0.0;
    amrex::Real J0=2.0;
    int axial_dir=0;
    amrex::Real E0=0.2;
    amrex::Real jco=0.0003;
    amrex::Real alpha=50;
};
#endif
