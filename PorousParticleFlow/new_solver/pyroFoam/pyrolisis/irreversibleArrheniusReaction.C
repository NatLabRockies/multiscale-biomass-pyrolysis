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
    m_reactant_index(m_reactant_stoch.size(), -1),
    m_product_index(m_product_stoch.size(), -1)
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

scalarField irreversibleArrheniusReaction::computeMolarSources(const scalar& T, const scalarField& species) const
{
    scalar R = computeReactionRate(T, species);

    scalarField molarSources(species.size(), 0.);

    forAll(m_reactant_index, id)
    {
        molarSources[m_reactant_index[id]] = -m_reactant_stoch[id] * R;
    }

    forAll(m_product_index, id)
    {
        molarSources[m_product_index[id]] = m_product_stoch[id] * R;
    }

    return molarSources;
}
