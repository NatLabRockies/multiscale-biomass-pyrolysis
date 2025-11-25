#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_TimeIntegrator.H>
#include <AMReX_MLTensorOp.H>
#include <ProbParm.H>
#include <TranspReact.H>
#include <Species.H>
#include <Transport.H>
#include <UserBCs.H>
#include <Chemistry.H>
#include <AMReX_MLABecLaplacian.H>
#include <compute_adv_flux.H>

/*Advances the chemisty state from time -> time + dt_lev
This is done using TimeIntegrator from AMReX, which can
be either implicit or explicit depending on the parameters
in the input file. Implicit currently requires an existing
SUNDIALS installation, with USE_SUNDIALS=TRUE at compile time,
along with setting the following in the input file

Explicit (no sundials):
   integration.type = RungeKutta
   integration.rk.type = 3 (for 3rd order)
Implicit (with sundials):
   integration.type = SUNDIALS
   integration.sundials.type = ERK
 */
void TranspReact::chemistry_advance(int lev, Real time, Real dt_lev,
       MultiFab &adsrc_lev, 
                                             MultiFab &phi_old_lev, MultiFab &phi_new_lev)
{
    MultiFab& S_new = phi_new_lev; // old value
    MultiFab& S_old = phi_old_lev; // current value

    auto rhs_function = [&] ( Vector<MultiFab> & dSdt_vec, 
                             const Vector<MultiFab>& S_vec, const Real time) {
        auto & dSdt = dSdt_vec[0];
        MultiFab S(S_vec[0], amrex::make_alias, 0, S_vec[0].nComp());
        update_rxnsrc_at_level(lev, S, dSdt, time);
        amrex::MultiFab::Saxpy(dSdt, 1.0, adsrc_lev, 0, 0, NUM_SPECIES, 0);
                
    };
    
    auto rhs_null_function = [&] ( Vector<MultiFab> & dSdt_vec, 
                             const Vector<MultiFab>& S_vec, const Real time) {
        auto & dSdt = dSdt_vec[0];
        dSdt.setVal(0.0);
                
    };
    Vector<MultiFab> state_old, state_new;

    // This term has the current state
    state_old.push_back(MultiFab(S_old, amrex::make_alias, 0, S_new.nComp()));
    // This is where the integrator puts the new state, hence aliased to S_new
    state_new.push_back(MultiFab(S_new, amrex::make_alias, 0, S_new.nComp()));
    // Define the integrator
    TimeIntegrator<Vector<MultiFab>> integrator(state_old);
    if(integration_type=="SUNDIALS")
    {
        if(integration_sd_type=="ERK" || integration_sd_type=="DIRK")
        {
                integrator.set_rhs(rhs_function);
        }
        else if(integration_sd_type=="IMEX-RK")
        {
           //only implicit
           integrator.set_imex_rhs(rhs_function,rhs_null_function);
        }
        else if(integration_sd_type=="EX-MRI" || 
                integration_sd_type=="IM-MRI" ||
                integration_sd_type=="IMEX-MRI")
        {
           //always assume reaction is the fast rhs
           integrator.set_rhs(rhs_null_function);
           integrator.set_fast_rhs(rhs_function);
        }
        else
        {
          amrex::Abort("Wrong sundials integration type in transpreact\n");
        }
    }
    // Advance from time to time + dt_lev
    //S_new/phi_new should have the new state
    integrator.advance(state_old, state_new, time, dt_lev); 
}

void TranspReact::update_advsrc_at_all_levels(int specid,Vector<MultiFab>& Sborder,
                                              Vector<MultiFab>& adv_src, 
                                              amrex::Real cur_time)
{
    int time=cur_time;
    ProbParm const* localprobparm = d_prob_parm;

    Vector< Array<MultiFab,AMREX_SPACEDIM> > flux(finest_level+1);
    
    for(int lev=0;lev<=finest_level;lev++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);

            flux[lev][idim].define(ba, dmap[lev], 1, 0);
            flux[lev][idim].setVal(0.0);
        }
    }

    for(int lev=0; lev <= finest_level; lev++)
    {
        adv_src[lev].setVal(0.0);
        flux[lev][0].setVal(0.0);
#if AMREX_SPACEDIM > 1
        flux[lev][1].setVal(0.0);
#if AMREX_SPACEDIM == 3
        flux[lev][2].setVal(0.0);
#endif
#endif
    }

    for(int lev=0; lev <= finest_level; lev++)
    {
        compute_scalar_advection_flux(specid, lev, Sborder[lev], 
                                      flux[lev], all_bcs_lo[specid], 
                                      all_bcs_hi[specid], 
                                      cur_time);
    }

    // =======================================================
    // Average down the fluxes before using them
    // =======================================================
    for (int lev = finest_level; lev > 0; lev--)
    {
        average_down_faces(amrex::GetArrOfConstPtrs(flux[lev  ]),
                           amrex::GetArrOfPtrs(flux[lev-1]),
                           refRatio(lev-1), Geom(lev-1));
    }

    for(int lev=0;lev<=finest_level;lev++)
    {
        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();

        for (MFIter mfi(adv_src[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            Array4<Real> advsrc_arr = adv_src[lev].array(mfi);

            Array4<Real> sborder_arr = Sborder[lev].array(mfi);
            GpuArray<Array4<Real>, AMREX_SPACEDIM> flux_arr{
                AMREX_D_DECL(flux[lev][0].array(mfi), 
                             flux[lev][1].array(mfi), flux[lev][2].array(mfi))};

            // update residual
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                advsrc_arr(i, j, k) = (flux_arr[0](i, j, k) - flux_arr[0](i + 1, j, k)) / dx[0];
#if AMREX_SPACEDIM > 1
                advsrc_arr(i,j,k) += (flux_arr[1](i, j, k) - flux_arr[1](i, j + 1, k)) / dx[1];
#if AMREX_SPACEDIM == 3
                advsrc_arr(i,j,k) += (flux_arr[2](i, j, k) - flux_arr[2](i, j, k + 1)) / dx[2]; 
#endif
#endif
            });
        }
    }

    // Additional source terms for axisymmetric geometry
    if(geom[0].IsRZ())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dx = geom[lev].CellSizeArray();
            auto prob_lo = geom[lev].ProbLoArray();

            for (MFIter mfi(adv_src[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

                Array4<Real> s_arr = Sborder[lev].array(mfi);
                Array4<Real> advsrc_arr = adv_src[lev].array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    // calculate r which is always x
                    // ideally x will be positive for axisymmetric cases
                    amrex::Real rval = amrex::Math::abs(prob_lo[0]+(i+0.5)*dx[0]);

                    // Calculate the advective source term component
                    IntVect iv_cell{AMREX_D_DECL(i, j, k)};
                    amrex::Real velr=tr_transport::compute_vel(iv_cell,0,specid,s_arr,dx,time,*localprobparm);
                    advsrc_arr(i,j,k) -= velr / rval;
                });
            }
        }
    }

}

void TranspReact::compute_scalar_advection_flux(int specid,int lev, MultiFab& Sborder, 
                                                Array<MultiFab,AMREX_SPACEDIM>& flux, 
                                                Vector<int>& bc_lo, Vector<int>& bc_hi,
                                                Real current_time)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int captured_specid = specid;
    //class member variable
    int captured_hyporder = hyp_order; 

    amrex::Real lev_dt=dt[lev];

    // Get the boundary ids
    const int* domlo_arr = geom[lev].Domain().loVect();
    const int* domhi_arr = geom[lev].Domain().hiVect();

    GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
    GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

    GpuArray<int,AMREX_SPACEDIM> bclo={AMREX_D_DECL(bc_lo[0], bc_lo[1], bc_lo[2])};
    GpuArray<int,AMREX_SPACEDIM> bchi={AMREX_D_DECL(bc_hi[0], bc_hi[1], bc_hi[2])};


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Box bx_x = mfi.nodaltilebox(0);
#if AMREX_SPACEDIM > 1
            Box bx_y = mfi.nodaltilebox(1);
#if AMREX_SPACEDIM == 3
            Box bx_z = mfi.nodaltilebox(2);
#endif
#endif

            Real time = current_time; // for GPU capture

            Array4<Real> sborder_arr = Sborder.array(mfi);

            GpuArray<Array4<Real>, AMREX_SPACEDIM> 
            flux_arr{AMREX_D_DECL(flux[0].array(mfi), 
                                  flux[1].array(mfi), flux[2].array(mfi))};

            amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 0, captured_specid, sborder_arr, 
                             bclo, bchi, domlo, domhi, flux_arr[0], 
                             time, dx, lev_dt, *localprobparm, captured_hyporder);
            });

#if AMREX_SPACEDIM > 1
            amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 1, captured_specid, sborder_arr, 
                             bclo, bchi, domlo, domhi, flux_arr[1], 
                             time, dx, lev_dt, *localprobparm, captured_hyporder);
            });

#if AMREX_SPACEDIM == 3
            amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 2, captured_specid, sborder_arr, 
                             bclo, bchi, domlo, domhi, flux_arr[2], 
                             time, dx, lev_dt, *localprobparm, captured_hyporder);
            });
#endif
#endif
        }
    }
}

void TranspReact::update_rxnsrc_at_level(int lev, MultiFab &S, MultiFab &dSdt, amrex::Real time)
{
    ProbParm const* localprobparm = d_prob_parm;
    dSdt.setVal(0.0);
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();


    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real> s_arr = S.array(mfi);
        Array4<Real> dsdt_arr = dSdt.array(mfi);

        // update residual
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                tr_reactions::production_rate(i, j, k, s_arr, dsdt_arr,
                        prob_lo, prob_hi, dx, time, *localprobparm);

                });
    }

}

void TranspReact::update_rxnsrc_at_all_levels(Vector<MultiFab>& Sborder,
        Vector<MultiFab>& rxn_src, 
        amrex::Real cur_time)
{
    amrex::Real time = cur_time;
    ProbParm const* localprobparm = d_prob_parm;

    // Zero out reactive source MFs
    for(int lev=0; lev <= finest_level; lev++)
    {
        rxn_src[lev].setVal(0.0);
    }

    for(int lev=0;lev<=finest_level;lev++)
    {
        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();

        for (MFIter mfi(rxn_src[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> sborder_arr = Sborder[lev].array(mfi);
            Array4<Real> rxn_arr = rxn_src[lev].array(mfi);

            // update residual
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    tr_reactions::production_rate(i, j, k, sborder_arr, rxn_arr,
                            prob_lo, prob_hi, dx, time, *localprobparm);

                    });
        }
    }
}

void TranspReact::implicit_solve_scalar(Real current_time, Real dt, int spec_id, 
        Vector<MultiFab>& Sborder, 
        Vector<MultiFab>& Sborder_old,
        Vector<MultiFab>& rxn_src, 
        Vector<MultiFab>& adv_src) 
{
    BL_PROFILE("TranspReact::implicit_solve_species(" + std::to_string( spec_id ) + ")");


    // FIXME: add these as inputs
    int max_coarsening_level = linsolve_max_coarsening_level;
    int linsolve_verbose = 1;
    int captured_spec_id=spec_id;
    int steady_solve=steadyspec[spec_id];

    //==================================================
    // amrex solves
    // read small a as alpha, b as beta

    //(A a - B del.(b del)) phi = f
    //
    // A and B are scalar constants
    // a and b are scalar fields
    // f is rhs
    // in this case: A=0,a=0,B=1,b=conductivity
    // note also the negative sign
    //====================================================
    ProbParm const* localprobparm = d_prob_parm;

    const Real tol_rel = linsolve_reltol;
    const Real tol_abs = linsolve_abstol;

    // set A and B, A=1/dt, B=1
    Real ascalar = 1.0;
    Real bscalar = 1.0;

#ifdef AMREX_USE_HYPRE
    if(use_hypre)
    {
        amrex::Print()<<"using hypre\n";
    }
#endif

    // default to inhomogNeumann since it is defaulted to flux = 0.0 anyways
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_lo 
        = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_hi 
        = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    int mixedbc=0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        //lower side bcs
        if (all_bcs_lo[spec_id][idim] == PERBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Periodic;
        }
        if (all_bcs_lo[spec_id][idim] == DIRCBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Dirichlet;
        }
        if (all_bcs_lo[spec_id][idim] == HNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Neumann;
        }
        if (all_bcs_lo[spec_id][idim] == IHNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::inhomogNeumann;
        }
        if (all_bcs_lo[spec_id][idim] == ROBINBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
        if (all_bcs_lo[spec_id][idim] == AXISBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Neumann;
        }

        //higher side bcs
        if (all_bcs_hi[spec_id][idim] == PERBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Periodic;
        }
        if (all_bcs_hi[spec_id][idim] == DIRCBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Dirichlet;
        }
        if (all_bcs_hi[spec_id][idim] == HNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Neumann;
        }
        if (all_bcs_hi[spec_id][idim] == IHNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::inhomogNeumann;
        }
        if (all_bcs_hi[spec_id][idim] == ROBINBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
        if (all_bcs_hi[spec_id][idim] == AXISBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Neumann;
        }
    }

    Vector<MultiFab> specdata(finest_level+1);
    Vector<MultiFab> acoeff(finest_level+1);
    Vector<MultiFab> bcoeff(finest_level+1);
    Vector<MultiFab> solution(finest_level+1);
    Vector<MultiFab> rhs(finest_level+1);

    Vector<MultiFab> robin_a(finest_level+1);
    Vector<MultiFab> robin_b(finest_level+1);
    Vector<MultiFab> robin_f(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        specdata[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        bcoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    linsolve_ptr.reset(new MLABecLaplacian(Geom(0,finest_level), 
                                           boxArray(0,finest_level), 
                                           DistributionMap(0,finest_level), info));
    MLMG mlmg(*linsolve_ptr);
    mlmg.setMaxIter(linsolve_maxiter);
    mlmg.setVerbose(linsolve_verbose);

#ifdef AMREX_USE_HYPRE
        if (use_hypre)
        {
            mlmg.setHypreOptionsNamespace("tr.hypre");
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        }
#endif
    linsolve_ptr->setDomainBC(bc_linsolve_lo, bc_linsolve_hi);
    linsolve_ptr->setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        // Copy args (FabArray<FAB>& dst, FabArray<FAB> const& src, 
        // int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
        specdata[ilev].setVal(0.0);
        amrex::Copy(specdata[ilev], Sborder_old[ilev], captured_spec_id, 
                    0, 1, num_grow);

        acoeff[ilev].setVal(0.0);
        if(!steady_solve)
        {
           acoeff[ilev].setVal(1.0/dt);
        }

        bcoeff[ilev].setVal(1.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        rhs[ilev].setVal(0.0);
        amrex::MultiFab::Saxpy(rhs[ilev], 1.0, rxn_src[ilev], spec_id, 0, 1, 0);

        if(!steady_solve)
        {
            amrex::MultiFab::Saxpy(rhs[ilev], 1.0/dt, specdata[ilev], 0, 0, 1, 0);
        }
        if(do_advection)
        {
            amrex::MultiFab::Saxpy(rhs[ilev], 1.0, adv_src[ilev], 0, 0, 1, 0);
        }

        amrex::Copy(specdata[ilev], Sborder[ilev], captured_spec_id, 
                0, 1, num_grow);

        solution[ilev].setVal(0.0);
        
        //for some reason, the previous solution initialization
        //fails MLMG, must be a tolerance thing
        if(!steady_solve)
        {
            amrex::MultiFab::Copy(solution[ilev], specdata[ilev], 0, 0, 1, 0);
        }

        // fill cell centered diffusion coefficients and rhs
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Real time = current_time; // for GPU capture

            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Array4<Real> bcoeff_arr = bcoeff[ilev].array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                bcoeff_arr(i,j,k)=tr_transport::specDiff(i,j,k,captured_spec_id,sb_arr,
                                                         dx,prob_lo,prob_hi,time,*localprobparm); 

            });
        }



        // average cell coefficients to faces, this includes boundary faces
        Array<MultiFab, AMREX_SPACEDIM> face_bcoeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoeff[ilev].boxArray(), 
                                                IntVect::TheDimensionVector(idim));
            face_bcoeff[idim].define(ba, bcoeff[ilev].DistributionMap(), 1, 0);
        }
        // true argument for harmonic averaging
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoeff), 
                                          bcoeff[ilev], geom[ilev], true);


        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> bc_arr = specdata[ilev].array(mfi);
            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Real time = current_time; // for GPU capture

            Array4<Real> robin_a_arr = robin_a[ilev].array(mfi);
            Array4<Real> robin_b_arr = robin_b[ilev].array(mfi);
            Array4<Real> robin_f_arr = robin_f[ilev].array(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (!geom[ilev].isPeriodic(idim))
                {
                    //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                    //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                    //so the ghost cell index at left side is i-1 while it is i on the right
                    if (bx.smallEnd(idim) == domain.smallEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryLo(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            tr_boundaries::species_bc(i, j, k, idim, -1, 
                                                      captured_spec_id, sb_arr, bc_arr, robin_a_arr,
                                                      robin_b_arr, robin_f_arr, 
                                                      prob_lo, prob_hi, dx, time, *localprobparm);
                        });
                    }
                    if (bx.bigEnd(idim) == domain.bigEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            tr_boundaries::species_bc(i, j, k, idim, +1, 
                                                      captured_spec_id, sb_arr, bc_arr, robin_a_arr, 
                                                      robin_b_arr, robin_f_arr,
                                                      prob_lo, prob_hi, dx, time, *localprobparm);
                        });
                    }
                }
            }
        }

        linsolve_ptr->setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        linsolve_ptr->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));

        // bc's are stored in the ghost cells
        if(mixedbc)
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]), &(robin_a[ilev]), 
                                     &(robin_b[ilev]), &(robin_f[ilev]));
        }
        else
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]));
        }
    }

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

    //bound species density
    if(bound_specden)
    { 
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            amrex::Real minconc=min_species_conc; 
            for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> soln_arr = solution[ilev].array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                      if(soln_arr(i,j,k) < minconc)
                      {
                        soln_arr(i,j,k)=minconc;
                      } 
                });
            }
        }
    }

    // copy solution back to phi_new
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, spec_id, 1, 0);
    }
    
    Print()<<"Solved species:"<<allvarnames[spec_id]<<"\n";

    //clean-up
    specdata.clear();
    acoeff.clear();
    bcoeff.clear();
    solution.clear();
    rhs.clear();

    robin_a.clear();
    robin_b.clear();
    robin_f.clear();
}

void TranspReact::transform_variables(Vector<MultiFab>& Sborder,amrex::Real cur_time)
{
    amrex::Real time = cur_time;
    ProbParm const* localprobparm = d_prob_parm;
    
    for(int lev=0;lev<=finest_level;lev++)
    {
        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();
        const int* domlo_arr = geom[lev].Domain().loVect();
        const int* domhi_arr = geom[lev].Domain().hiVect();

        GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
        GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

        for (MFIter mfi(phi_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> sborder_arr = Sborder[lev].array(mfi);
            Array4<Real> phi_arr = phi_new[lev].array(mfi);

            // update residual
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                tr_reactions::transform(i, j, k, sborder_arr, phi_arr,
                                        prob_lo, prob_hi, domlo, domhi, dx, time, *localprobparm);

            });
        }
    }
}


