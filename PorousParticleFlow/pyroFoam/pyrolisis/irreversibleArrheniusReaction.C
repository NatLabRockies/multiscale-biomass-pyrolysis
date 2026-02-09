#include "irreversibleArrheniusReaction.h"

using namespace Foam;
using namespace pyrolisis;

irreversibleArrheniusReaction::irreversibleArrheniusReaction(const word& name, const dictionary& dict, const List<word>& species_name)
:
    m_name(name),
    m_reactant_stoch(dict.lookup("reactant_stochiometric_coefficients")),
    m_product_stoch(dict.lookup("product_stochiometric_coefficients")),
    m_A(readScalar(dict.lookup("A"))),
    m_Ta(readScalar(dict.lookup("Ta"))),
    m_beta(dict.lookupOrDefault<scalar>("beta",0.)),
    m_Q(dict.lookupOrDefault<scalar>("heatOfReaction",0.)),
    m_reactant_index(m_reactant_stoch.size(), -1),
    m_product_index(m_product_stoch.size(), -1),
    m_custom_reaction(dict.lookupOrDefault<bool>("customReaction", false))
{
    Info << "Reading reaction " << m_name << endl;

    List<word> reactant_list(dict.lookup("reactants"));
    List<word> product_list(dict.lookup("products"));

    if( reactant_list.size() != m_reactant_stoch.size() )
    {
        FatalErrorInFunction << "Number of reactants and stochiometric coefficients do not match" << endl;
    }

    if( product_list.size() != m_product_stoch.size() )
    {
        FatalErrorInFunction << "Number of products and stochiometric coefficients do not match" << endl;
    }

    // Populate index arrays
    forAll(species_name, specieI)
    {
        forAll(reactant_list, reactI)
        {
            if ( species_name[specieI] == reactant_list[reactI] )
            {
                m_reactant_index[reactI] = specieI;
            }
        }

        forAll(product_list, prodI)
        {
            if ( species_name[specieI] == product_list[prodI] )
            {
                m_product_index[prodI] = specieI;
            }
        }
    }

    forAll(m_reactant_index, id)
    {
        if (m_reactant_index[id] == -1)
        {
            FatalErrorInFunction << "Error: unknown reactant " << reactant_list[id] << "\n" << abort(FatalError);
        }
    }

    forAll(m_product_index, id)
    {
        if (m_product_index[id] == -1)
        {
            FatalErrorInFunction << "Error: unknown product " << product_list[id] << "\n" << abort(FatalError);
        }
    }
    
    if (m_custom_reaction) 
    {
	    m_mf = scalar (dict.get<scalar>("mf"));
    }

}

irreversibleArrheniusReaction::~irreversibleArrheniusReaction()
{}


scalar irreversibleArrheniusReaction::computeReactionRate(const scalar& T, const scalarField& species) const
{
    /* Compute Arrhenius rate */
    scalar k = pow(T, m_beta) * m_A * exp(-m_Ta / T);

    /* Compute and return reaction rate */
    scalar prodCnu(1.0);

    forAll(m_reactant_stoch, id)
    {
        prodCnu *= pow(species[m_reactant_index[id]], m_reactant_stoch[id]);
    }

    return k * prodCnu;
}

scalar irreversibleArrheniusReaction::computeCustomReactionRate(const scalar& T, const scalarField& species, const scalarField& molarWeight, const scalar& m0) const
{
    
    //Pout<<"computing custom reaction rate"<<endl;
	/* Compute Arrhenius rate */
    scalar k = pow(T, m_beta) * m_A * exp(-m_Ta / T);

    /* Compute and return reaction rate */
    scalar prodCnu(1.0);

    scalar mtot = 0.0;
    scalar sum_mw = 0.0;
    scalar sum_m0 = m0;
    //Pout<<"computing reactant mass rate"<<endl;
    forAll(m_reactant_stoch, id)
    {
	//Pout<<"computing concentration for species "<<id<<endl;
	prodCnu *= pow(species[m_reactant_index[id]], m_reactant_stoch[id]);
	scalar conc = species[m_reactant_index[id]];
	scalar m = conc*molarWeight[m_reactant_index[id]];
	sum_mw += molarWeight[m_reactant_index[id]]; 
	//sum_m0 += m0[m_reactant_index[id]]; 
	mtot+=m;
	//Pout << "LDPE concentration is" << conc << endl;
	//Pout << "m0 is" << m0 << endl;
	//Pout << "sum_m0 is" << sum_m0 << endl;
    }
    //Pout << "LDPE mass is" << mtot << endl;
    scalar mfinal = sum_m0*m_mf;
    //Pout << "final mass is" <<mfinal<<endl;

    scalar alpha = 0;
    if (mfinal > 0)
    {
	    alpha = max(mag((sum_m0-mtot)/(sum_m0-mfinal)), 1e-6);
    }

    //Pout << "alpha is" << alpha << endl;

    scalar falpha = 2*(1-alpha)*std::sqrt(-std::log(max(mag(1-alpha), 1e-6)));
    //scalar falpha = std::log(max(mag(1-alpha), 1e-6));
    //Pout << "falpha is" << falpha << endl;


    scalar mscale = sum_m0-mfinal;

    return k * falpha * mscale / sum_mw;
}

scalar irreversibleArrheniusReaction::computeSources(
    const scalar& T,
    const scalarField& molarFraction,
    const scalarField& molarWeight,
    const scalar& m0,
    scalarField& ndot
) const
{
    scalar R=0.0;
    if(m_custom_reaction){
	    R = computeCustomReactionRate(T, molarFraction, molarWeight, m0);
    }
    else
    {
	    R = computeReactionRate(T, molarFraction);
    }
    
    //Pout<<"R is: "<<R<<endl;
    
    //Pout<<"computing reactant ndot"<<endl;
    forAll(m_reactant_index, id)
    {
        ndot[m_reactant_index[id]] += -m_reactant_stoch[id] * R;
    }
    //Pout<<"computing product ndot"<<endl;

    forAll(m_product_index, id)
    {
        ndot[m_product_index[id]] += m_product_stoch[id] * R;
    }
    
    //Pout<<"computing mdot"<<endl;

    scalar mdot = 0.; // Mass reacted

    // Compute Kg/m3 reacted
    // Notice that the magnitude of the time derivative of the molar density
    // of the reactant is taken. This means that mdot represents the total
    // mass that underwent this reaction per unit time.
    forAll(m_reactant_index, id)
    {
        mdot += mag(ndot[m_reactant_index[id]]) * molarWeight[m_reactant_index[id]];
    }

    //Pout<<"mdot is: "<<mdot<<endl;

    return m_Q * mdot;
}
