#pragma once

#include "fvCFD.H"
#include "dictionary.H"

namespace Foam
{
namespace pyrolisis
{

class irreversibleArrheniusReaction
{
    word m_name;

    List<scalar> m_reactant_stoch;
    List<scalar> m_product_stoch;

    scalar m_A;
    scalar m_Ta;
    scalar m_beta;
    scalar m_Q; // Heat of reaction in J/Kg
    scalar m_mf=0; // final mass multiplier from experiment
    List<int> m_reactant_index;
    List<int> m_product_index;
    bool m_custom_reaction;
    scalar computeReactionRate (const scalar& T, const scalarField& species) const;
    scalar computeCustomReactionRate (const scalar& T, const scalarField& species, const scalarField& molarWeight, const scalar& m0) const;

public:

    irreversibleArrheniusReaction(const word& name, const dictionary& dict, const List<word>& species_name);
    ~irreversibleArrheniusReaction();

    scalar computeSources(
        const scalar& T,
        const scalarField& molarFraction,
        const scalarField& molarWeight,
        const scalar& m0,
        scalarField& ndot
    ) const;

    bool customReaction() {return m_custom_reaction;};

};

}
}
