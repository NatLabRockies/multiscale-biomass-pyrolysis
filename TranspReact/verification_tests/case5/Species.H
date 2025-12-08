#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#define MAT1_ID 0
#define MAT2_ID 1
#define MAT3_ID 2
#define TEMP_ID 3
#define NUM_SPECIES 4

namespace tr_species
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
}

#endif
