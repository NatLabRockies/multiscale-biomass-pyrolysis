#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#define CO2_ID 0
#define Hp_ID 1
#define OHm_ID 2
#define CO_ID 3
#define HCO3m_ID 4
#define CO32m_ID 5
#define Kp_ID 6
#define H2_ID 7
#define TMAp_ID 8

#define PHI_ID 9
#define EFX_ID 10
#define EFY_ID 11
#define EFZ_ID 12
#define MEM_ID 13
#define CAT_ID 14
#define IM_ID 15
#define LE_ID 16
#define NUM_SPECIES 17
#define NUM_TR_SPECIES 9

namespace tr_species
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
}

#endif
