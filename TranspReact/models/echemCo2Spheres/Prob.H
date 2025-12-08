#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Chemistry.H>
#include <ProbParm.H>
#include <UnivConstants.H>
#include <VarDefines.H>
#include <UserBCs.H>

using namespace amrex;

    AMREX_INLINE
void amrex_probinit (ProbParm & h_prob_parm, ProbParm & d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");
    //pp.query("right_temp", h_prob_parm.right_temp );
    pp.query("Kmem", h_prob_parm.Kmem );
    pp.query("Kim", h_prob_parm.Kim );
    pp.query("Kcat", h_prob_parm.Kcat );
    pp.query("Kle", h_prob_parm.Kle );
    pp.query("sourcefac", h_prob_parm.sourcefac );
    pp.query("tmap_conc",h_prob_parm.tmap_conc);
    pp.query("koh_init",h_prob_parm.koh_init);
    pp.query("co2_boundary",h_prob_parm.co2_boundary);
    pp.query("catelecpot",h_prob_parm.catelecpot);
    pp.query("i0_co",h_prob_parm.i0_co);
    pp.query("i0_h2",h_prob_parm.i0_h2);
    pp.query("E0_co",h_prob_parm.E0_co);
    pp.query("E0_h2",h_prob_parm.E0_h2);
    pp.query("im_thickness",h_prob_parm.im_thickness);
    pp.query("imfrac",h_prob_parm.imfrac);
    pp.query("catfrac",h_prob_parm.catfrac);
    pp.query("current_input",h_prob_parm.current_input);
    pp.query("activate_surf_chem",h_prob_parm.activate_surf_chem);
    pp.query("activate_im_coverage",h_prob_parm.activate_im_coverage);
    pp.query("cat_conductivity",h_prob_parm.cat_conductivity);
    pp.query("elyte_const_conductivity",h_prob_parm.elyte_const_conductivity);
    pp.query("co2star",h_prob_parm.co2star);
    amrex::Real kD=1e-9;
    pp.query("kD",kD);
    pp.query("bvsimple",h_prob_parm.bvsimple);
    pp.query("alpha_co",h_prob_parm.alpha_co);
    pp.query("alpha_h2",h_prob_parm.alpha_h2);

    pp.query("current_test",h_prob_parm.current_test);
    pp.query("co2conc_test",h_prob_parm.co2conc_test);

    for (int c=0;c<NUM_SPECIES;c++)
    {
        h_prob_parm.dcoeff[c]=1.0;
        h_prob_parm.Zi[c]=0.0;
        h_prob_parm.tplus[c]=0.0;
        h_prob_parm.kD[c]=0.0;
    }

    //all in m2/s
    /*h_prob_parm.dcoeff[CO2_ID]=4.3716e-9;
    h_prob_parm.dcoeff[CO_ID]=4.3716e-9;
    h_prob_parm.dcoeff[Hp_ID]=1.156e-8;
    h_prob_parm.dcoeff[OHm_ID]=9.1624e-9;
    h_prob_parm.dcoeff[HCO3m_ID]=2.3463e-9;
    h_prob_parm.dcoeff[CO32m_ID]=1.6774e-9;
    h_prob_parm.dcoeff[Kp_ID]=4.425e-9;
    h_prob_parm.dcoeff[TMAp_ID]=1e-20;*/

    h_prob_parm.dcoeff[CO2_ID]=4.3716e-9;
    h_prob_parm.dcoeff[CO_ID]=4.3716e-9;
    h_prob_parm.dcoeff[Hp_ID]=1.156e-8;
    h_prob_parm.dcoeff[OHm_ID]=9.1624e-9;
    h_prob_parm.dcoeff[HCO3m_ID]=2.3463e-9;
    h_prob_parm.dcoeff[CO32m_ID]=1.6774e-9;
    h_prob_parm.dcoeff[Kp_ID]=4.425e-9;
    h_prob_parm.dcoeff[H2_ID]=4.425e-9;
    h_prob_parm.dcoeff[TMAp_ID]=1e-20;

    //charge number
    h_prob_parm.Zi[Hp_ID]=1.0;
    h_prob_parm.Zi[OHm_ID]=-1.0;
    h_prob_parm.Zi[HCO3m_ID]=-1.0;
    h_prob_parm.Zi[CO32m_ID]=-2.0;
    h_prob_parm.Zi[Kp_ID]=1.0;
    h_prob_parm.Zi[TMAp_ID]=1.0;
    
    h_prob_parm.tplus[OHm_ID]=-0.5;
    h_prob_parm.tplus[Kp_ID]=0.5;
    
    h_prob_parm.kD[OHm_ID]=kD;
    h_prob_parm.kD[Kp_ID]=kD;
    
    amrex::Vector<Real> plo;
    amrex::Vector<Real> phi;
    amrex::ParmParse ppgeom("geometry");
    ppgeom.queryarr("prob_lo",plo);
    ppgeom.queryarr("prob_hi",phi);


    pp.query("nsph_x",h_prob_parm.nsph_x);
    pp.query("nsph_y",h_prob_parm.nsph_y);
    pp.query("nsph_z",h_prob_parm.nsph_z);
    pp.query("sph_twiddle_fac",h_prob_parm.sph_twiddle_fac);
    h_prob_parm.nsph=h_prob_parm.nsph_x*h_prob_parm.nsph_y*h_prob_parm.nsph_z;
    if(h_prob_parm.nsph > 0)
    {
        amrex::Vector<amrex::Real> locx;
        amrex::Vector<amrex::Real> locy;
        amrex::Vector<amrex::Real> locz;
        amrex::Vector<amrex::Real> imcoverage;

        locx.resize(h_prob_parm.nsph);
        locy.resize(h_prob_parm.nsph);
        locz.resize(h_prob_parm.nsph);

        amrex::Real dx_sph=(phi[0]-plo[0])/h_prob_parm.nsph_x;
        amrex::Real dy_sph=(phi[1]-plo[1])/h_prob_parm.nsph_y;
        amrex::Real dz_sph=(phi[2]-plo[2])/h_prob_parm.nsph_z;

        pp.queryarr("imcoverage",imcoverage);
        pp.query("sphererad",h_prob_parm.sphererad);

        int nsph_x=h_prob_parm.nsph_x;
        int nsph_y=h_prob_parm.nsph_y;
        int nsph_z=h_prob_parm.nsph_z;

        for(int k=0;k<nsph_z;k++)
        {
           for(int j=0;j<nsph_y;j++)
           {
                for(int i=0;i<nsph_x;i++)
                {
                    int index=k*(nsph_y*nsph_x)+j*nsph_x+i;
                    if(h_prob_parm.sph_twiddle_fac>0.0)
                    {
                        locx[index]=amrex::RandomNormal(plo[0]+(i+0.5)*dx_sph,h_prob_parm.sph_twiddle_fac*dx_sph);
                        locy[index]=amrex::RandomNormal(plo[1]+(j+0.5)*dy_sph,h_prob_parm.sph_twiddle_fac*dy_sph);
                        locz[index]=amrex::RandomNormal(plo[2]+(k+0.5)*dz_sph,h_prob_parm.sph_twiddle_fac*dz_sph);
                    }
                    else
                    {
                        locx[index]=plo[0]+(i+0.5)*dx_sph;
                        locy[index]=plo[1]+(j+0.5)*dy_sph;
                        locz[index]=plo[2]+(k+0.5)*dz_sph;
                    }
                }
           }
        }

        h_prob_parm.loc_x  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
        h_prob_parm.loc_y  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
        h_prob_parm.loc_z  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));
        h_prob_parm.imcoverage  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm.nsph*sizeof(amrex::Real));

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         locx.begin(),
                         locx.end(),
                         h_prob_parm.loc_x);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         locy.begin(),
                         locy.end(),
                         h_prob_parm.loc_y);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         locz.begin(),
                         locz.end(),
                         h_prob_parm.loc_z);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         imcoverage.begin(),
                         imcoverage.end(),
                         h_prob_parm.imcoverage);

    }
    
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#endif

}

AMREX_GPU_DEVICE
AMREX_INLINE
void initdomaindata(Box const& bx, Array4<Real> const& phi, 
                    GeometryData const& geomdata,
                    ProbParm *prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    amrex::Real phi_e=-2.0;
    amrex::Real phi_i=tr_boundaries::bvsolve_complex(-1.0,phi_e,prob_parm->co2conc_test,
                                   prob_parm->T,prob_parm->current_test,
                                   *prob_parm);

    amrex::Print()<<"vdiff:"<<phi_e-phi_i<<"\n";
    amrex::Real i_co,i_h2;
    amrex::Real bvcur=tr_boundaries::bvfunc_complex(phi_i,phi_e,prob_parm->co2conc_test,
                                     prob_parm->T,*prob_parm,i_co,i_h2);

    amrex::Print()<<"bvcur,i_co,i_h2:"<<bvcur<<"\t"<<i_co<<"\t"<<i_h2<<"\n";

    int axdir=0;
    amrex::Real cat_interface_loc=prob_lo[axdir]+prob_parm->catfrac*
    (prob_hi[axdir]-prob_lo[axdir]);
    amrex::Real initial_currentden=0.0;

    int nx = prob_parm->nx_xct;
    int ny = prob_parm->ny_xct;
    int nz = prob_parm->nz_xct;
    amrex::Real dx_xct[AMREX_SPACEDIM];
    dx_xct[0]=(prob_hi[0]-prob_lo[0])/prob_parm->nx_xct;
    dx_xct[1]=(prob_hi[1]-prob_lo[1])/prob_parm->ny_xct;
    dx_xct[2]=(prob_hi[2]-prob_lo[2])/prob_parm->nz_xct;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                IntVect iv(AMREX_D_DECL(i,j,k));
                for(int c=0;c<NVAR;c++)
                {
                    phi(i,j,k,c)=0.0;
                }
                phi(i,j,k,CMASK_ID) = 1.0;
                
                for(int kk=0;kk<2;kk++)
                {
                    for(int jj=0;jj<2;jj++)
                    {
                        for(int ii=0;ii<2;ii++)
                        {
                            Real x[3]={prob_lo[0]+(i+ii)*dx[0],
                                prob_lo[1]+(j+jj)*dx[1],prob_lo[0]+(k+kk)*dx[2]};

                            for(int s=0;s<prob_parm->nsph;s++)
                            {
                                amrex::Real xcent[3]={prob_parm->loc_x[s],
                                    prob_parm->loc_y[s],prob_parm->loc_z[s]};

                                amrex::Real dist2=(x[0]-xcent[0])*(x[0]-xcent[0]);
                                dist2+=(x[1]-xcent[1])*(x[1]-xcent[1]);
                                dist2+=(x[2]-xcent[2])*(x[2]-xcent[2]);

                                amrex::Real dist2xy=(x[0]-xcent[0])*(x[0]-xcent[0]);
                                dist2xy+=(x[1]-xcent[1])*(x[1]-xcent[1]);
                                
                                if(dist2 < prob_parm->sphererad*prob_parm->sphererad)
                                //if(dist2xy < prob_parm->sphererad*prob_parm->sphererad)
                                {
                                    phi(i,j,k,CAT_ID)+=1.0;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                /*for(int kk=0;kk<2;kk++)
                {
                    for(int jj=0;jj<2;jj++)
                    {
                        for(int ii=0;ii<2;ii++)
                        {
                            Real x[3]={prob_lo[0]+(i+ii)*dx[0],
                                prob_lo[1]+(j+jj)*dx[1],prob_lo[0]+(k+kk)*dx[2]};

                            if(x[2] < prob_lo[2]+0.333*(prob_hi[2]-prob_lo[2]))
                            {
                                phi(i,j,k,CAT_ID)+=1.0;
                            }
                        }
                    }
                }*/

                phi(i,j,k,CAT_ID)*=0.25;
                if(phi(i,j,k,CAT_ID)>0.0)
                {
                    phi(i,j,k,CAT_ID)=1.0;
                }
                phi(i,j,k,LE_ID)=(1.0-phi(i,j,k,CAT_ID));
                phi(i,j,k,CMASK_ID)=phi(i,j,k,LE_ID);

                phi(i,j,k,CO2_ID)=phi(i,j,k,LE_ID)*prob_parm->co2_boundary*0.01;
                phi(i,j,k,OHm_ID)=phi(i,j,k,LE_ID)*prob_parm->koh_init;
                phi(i,j,k,Kp_ID)=phi(i,j,k,LE_ID)*prob_parm->koh_init;

                //solid potential
                if(phi(i,j,k,CAT_ID)==1.0)
                {
                    phi(i,j,k,PHI_ID)=prob_parm->catelecpot;
                }
            }
        }
    }
}
#endif
