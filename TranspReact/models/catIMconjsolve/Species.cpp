#include<Species.H>

namespace tr_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);

    void init()
    {
        specnames[CO2_ID]="CO2";
        specnames[Hp_ID]="Hp";
        specnames[OHm_ID]="OHm";
        specnames[CO_ID]="CO";
        specnames[HCO3m_ID]="HCO3m";
        specnames[CO32m_ID]="CO32m";
        specnames[HCOOm_ID]="HCOOm";
        specnames[Kp_ID]="Kp";
        specnames[HCOOH_ID]="HCOOH";
        specnames[TMAp_ID]="TMAp";
        specnames[PHI_ID]="Potential";
        specnames[EFX_ID]="Efieldx";
        specnames[EFY_ID]="Efieldy";
        specnames[EFZ_ID]="Efieldz";
        specnames[MEM_ID]="Membrane";
        specnames[CAT_ID]="Catalyst";
        specnames[IM_ID]="Ionomer";
        specnames[LE_ID]="LiqElectrolyte";
    }    
    void close()
    {
        specnames.clear();
    }
    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(specnames.begin(),specnames.end(),specname);
        if(it != specnames.end())
        {
            loc=it-specnames.begin();
        }
        return(loc);
    }
}
