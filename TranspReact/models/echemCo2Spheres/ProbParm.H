#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <Species.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real dcoeff[NUM_SPECIES]={0.0};
    amrex::Real Zi[NUM_SPECIES]={0.0};
    amrex::Real tplus[NUM_SPECIES]={0.0};
    amrex::Real kD[NUM_SPECIES]={0.0};
    amrex::Real T=333.15;
    amrex::Real Kmem=35.0;
    amrex::Real Kim=35.0;
    amrex::Real Kle=35.0;
    amrex::Real Kcat=1e6;
    amrex::Real sourcefac=1e-10;
    amrex::Real catelecpot=-1.0;
    int axial_dir=0;
    int flat_interface=0;
    amrex::Real tmap_conc=1002.0; //mol/m3
    amrex::Real koh_init=1000.0; //mol/m3
    amrex::Real co2_boundary=100.0; //mass transfer condition 
    amrex::Real i0_co=0.1352*1e5;
    amrex::Real i0_h2=0.09;
    amrex::Real E0_co=-0.889;
    amrex::Real E0_h2=-0.828;
    amrex::Real alpha_co=0.5;
    amrex::Real alpha_h2=0.5;
    amrex::Real imfrac=0.2;
    amrex::Real far_eff=0.95;
    amrex::Real co2star=0.06;

    amrex::Real current_input=1000; //A/m2
    amrex::Real cat_conductivity=10.0; //S/m
    amrex::Real elyte_const_conductivity=1.0; //S/m

    int activate_surf_chem=1;
    int activate_im_coverage=0;
    
    int nsph_x=0;
    int nsph_y=0;
    int nsph_z=0;
    int nsph=0;
    amrex::Real* loc_x=nullptr;
    amrex::Real* loc_y=nullptr;
    amrex::Real* loc_z=nullptr;
    amrex::Real *imcoverage=nullptr;
    amrex::Real sphererad=1e-6;
    amrex::Real im_thickness=0.1e-7;
    amrex::Real catfrac=0.1;
    amrex::Real sph_twiddle_fac=0.0;

    int bvsimple=1;
    // New fields for TIFF-based domain mapping
    int nx_xct = 1;  // X dimension of domain map
    int ny_xct = 1;  // Y dimension of domain map
    int nz_xct = 1;  // Z dimension of domain map
    int* domain_map = nullptr;  // pointer to flattened 3D domain map

    amrex::Real current_test=3000.0;
    amrex::Real co2conc_test=1.0;
    amrex::Real max_current_density=1e4;
};
#endif
