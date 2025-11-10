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

    List<int> m_reactant_index;
    List<int> m_product_index;

    scalar computeReactionRate (const scalar& T, const scalarField& species) const;

public:

    irreversibleArrheniusReaction(const word& name, const dictionary& dict, const List<word>& species_name);
    ~irreversibleArrheniusReaction();

    scalarField computeMolarSources(const scalar& T, const scalarField& species) const;
};

}
}
