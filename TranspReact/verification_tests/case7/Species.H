#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#define ETRODE_ID 0
#define ELYTE_ID 1
#define POT_ID 2
#define NUM_SPECIES 3

namespace tr_species
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
}

#endif
