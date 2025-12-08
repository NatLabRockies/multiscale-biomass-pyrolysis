#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real dc_mat1 = 0.1;
    amrex::Real dc_mat2 = 1.0;
    amrex::Real dc_mat3 = 0.5;
    amrex::Real  x1=0.2;
    amrex::Real x2=0.6;
    int axial_dir=0;
};
#endif
