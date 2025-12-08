#ifndef _Utils_H_
#define _Utils_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void update_interface_cells(int i,int j,int k,int specid,
                              Array4<Real> const &sb_arr,
                              Array4<Real> const &phi_arr)
{
    IntVect iv_cell{AMREX_D_DECL(i, j, k)};

    if(sb_arr(iv_cell,CMASK_ID)>0.0 && sb_arr(iv_cell,CMASK_ID)<1.0)
    {
        Real nadds=0.0;
        Real neighborsum=0.0;

#if AMREX_SPACEDIM == 3
        for(int kk=-1;kk<=1;kk++)
#endif
        {
#if AMREX_SPACEDIM > 1
            for(int jj=-1;jj<=1;jj++)
#endif
            {
                for(int ii=-1;ii<=1;ii++)
                {
                    IntVect iv_neigh{AMREX_D_DECL(i+ii, j+jj, k+kk)};
                    if(iv_neigh != iv_cell)
                    {
                        nadds=nadds+1.0;
                        neighborsum += sb_arr(iv_neigh,specid);
                    }
                }
            }
        }

        //nadds should be 26
        phi_arr(iv_cell,specid)=neighborsum/nadds;
    } 
}
#endif
