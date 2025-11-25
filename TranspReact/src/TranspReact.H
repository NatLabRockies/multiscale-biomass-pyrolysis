#ifndef TranspReact_H_
#define TranspReact_H_

#include <AMReX_MLMG.H>

#include <string>
#include <limits>
#include <memory>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <ProbParm.H>

#include <AMReX_MLABecLaplacian.H>
#include <UnivConstants.H>
#include <VarDefines.H>


using namespace amrex;

class TranspReact : public amrex::AmrCore
{
public:
    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    TranspReact();
    virtual ~TranspReact();

    // advance solution to final time
    void Evolve_coupled();
    void Evolve_split();

    void Evolve()
    {
        if(!split_chemistry)
        {
           Evolve_coupled();
        }
        else
        {
           Evolve_split();
        }
    }

    // initializes multilevel data
    void InitData();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse(int lev, amrex::Real time, 
                                        const amrex::BoxArray& ba, 
                                        const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel(int lev, amrex::Real time, 
                             const amrex::BoxArray& ba, 
                             const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel(int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch(int lev, amrex::Real time, 
                                         const amrex::BoxArray& ba, 
                                         const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, 
                          amrex::Real time, int ngrow) override;

   void compute_scalar_advection_flux(int specid, int lev, MultiFab& Sborder, 
                                             Array<MultiFab,AMREX_SPACEDIM>& flux, 
                                             Vector<int>& bc_lo, Vector<int>& bc_hi,
                                             Real time);
  
   void update_advsrc_at_all_levels(int specid,Vector<MultiFab>& Sborder,
                                              Vector<MultiFab>& adv_src, 
                                              amrex::Real cur_time);

    void update_rxnsrc_at_all_levels(Vector<MultiFab>& Sborder,
                                            Vector<MultiFab>& rxn_src, 
                                            amrex::Real cur_time);
    
    void implicit_solve_scalar(Real current_time, Real dt, int spec_id,
            Vector<MultiFab>& Sborder, Vector<MultiFab>& Sborder_old,
                                  Vector<MultiFab>& rxn_src, 
                                  Vector<MultiFab>& adv_src); 
    
    void transform_variables(Vector<MultiFab>& Sborder,amrex::Real cur_time);

    // compute dt from CFL considerations
    void find_transp_timescales(int lev,amrex::Real cur_time,
            amrex::Real& dt_adv,amrex::Real &dt_diff);

    void chemistry_advance(int lev, Real time, 
            Real dt_lev, MultiFab &adsrc_lev,
            MultiFab &phi_old_lev, 
            MultiFab &phi_new_lev);

    void update_rxnsrc_at_level(int lev, MultiFab &S, 
            MultiFab &dSdt, amrex::Real time);

    static ProbParm* h_prob_parm;
    static ProbParm* d_prob_parm;
   
    std::unique_ptr<amrex::MLABecLaplacian> linsolve_ptr;
    
    ////////////////
    amrex::Gpu::ManagedVector<Real> refine_phi;
    amrex::Gpu::ManagedVector<Real> refine_phigrad;
    amrex::Gpu::ManagedVector<int> refine_phi_comps;

    amrex::Vector<int> istep;     // which step?
    amrex::Vector<int> nsubsteps; // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;
    
    int hyp_order=1;
    int do_reactions=1;
    int do_transport=1;
    int do_advection=1;
    int transform_vars=0;
    int adaptive_dt=0;
    int split_chemistry=0;
    amrex::Real fixed_dt=1e-12;
    amrex::Real advective_cfl=0.5;
    amrex::Real diffusive_cfl=0.5;

    amrex::Real dt_max=1e50;
    amrex::Real dt_min=1e-50;
    std::string integration_type="ForwardEuler";
    std::string integration_sd_type="ERK";

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::MultiFab> phi_new;
    amrex::Vector<amrex::MultiFab> phi_old;

    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<BCRec> bcspec; // 1-component

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();


    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";


    amrex::Array<amrex::Vector<int>,NUM_SPECIES> all_bcs_lo;
    amrex::Array<amrex::Vector<int>,NUM_SPECIES> all_bcs_hi;

    amrex::Vector<int> pot_bc_lo{ROBINBC, ROBINBC, ROBINBC};
    amrex::Vector<int> pot_bc_hi{ROBINBC, ROBINBC, ROBINBC};
    
    amrex::Real linsolve_reltol = 1e-10;
    amrex::Real linsolve_abstol = 1e-10;
    amrex::Real linsolve_bot_reltol = 1e-10;
    amrex::Real linsolve_bot_abstol = -1.0;
    amrex::Real min_species_conc=1e-20;
    amrex::Vector<std::string> allvarnames;
    int linsolve_num_pre_smooth = 2;
    int linsolve_num_post_smooth = 2;
    int linsolve_num_final_smooth = 8;
    int linsolve_num_bottom_smooth = 0;

    int linsolve_maxiter=100;
    int linsolve_max_coarsening_level=10;
    int bound_specden=0;
    int ngrow_for_fillpatch=2;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;
    int num_timestep_correctors=1;
    int num_split_correctors=2;

#ifdef AMREX_USE_HYPRE
    int use_hypre=0;
#endif

    // plotfile prefix and frequency
    std::string plot_file{"plt"};
    int plot_int = -1;
    int plot_int_old = -1;
    amrex::Real plot_time = -1.0;

    // checkpoint prefix and frequency
    std::string chk_file{"chk"};
    int chk_int = -1;
    int chk_int_old = -1;
    amrex::Real chk_time = -1.0;

private:
    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo(int crse_lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch(int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch(int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // utility to copy in data from phi_old and/or phi_new into another multifab
    void GetData(int lev, amrex::Real time, 
                 amrex::Vector<amrex::MultiFab*>& data, amrex::Vector<amrex::Real>& datatime);

    // a wrapper for EstTimeStep(0
    void ComputeDt(amrex::Real cur_time,amrex::Real dt_diff,amrex::Real dt_adv);

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF() const;

    // write plotfile to disk
    void WritePlotFile(int plotfilenum) const;

    // write checkpoint file to disk
    void WriteCheckpointFile(int chkfilenum) const;

    // read checkpoint file from disk
    void ReadCheckpointFile();

    // utility to skip to next line in Header
    static void GotoNextLine(std::istream& is);
    
    int steadyspec[NUM_SPECIES]={0};
    int unsolvedspec[NUM_SPECIES]={0};

};

#endif
