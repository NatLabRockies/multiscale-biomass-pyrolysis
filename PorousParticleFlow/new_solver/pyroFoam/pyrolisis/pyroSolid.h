#pragma once

#include "fvCFD.H"
#include <utility>

namespace Foam {

namespace pyrolisis {

class irreversibleArrheniusReaction;

class pyroSolid {

private:
  const fvMesh &m_mesh;
  IOdictionary m_dict;
  volScalarField m_porosity;
  volScalarField m_rhoField;
  volScalarField m_T;
  volScalarField m_rhoCp;
  volScalarField m_Qdot; //Heat of reaction
 // volTensorField m_Kp;
  volScalarField m_htc;
  volScalarField m_kappa_tot;
  label m_nSubTimeSteps;
  scalar m_poreSize;
  List<word> m_speciesName;
  PtrList<volScalarField> m_species;
  PtrList<volScalarField> m_reactionRates;
  scalarField m_rho;
  scalarField m_cp;
  scalarField m_molWeight;
  scalarField m_kappa;
  List<std::pair<bool, label>> m_is_gas;
  PtrList<irreversibleArrheniusReaction> m_reactions;
  List<bool> m_addToPoro;
  List<scalar> m_formationEnthalpy;

  label getSpecieId(const word &name);

  void solveEnergy();

 // void updateKp();
  void updateHTC();

public:
  pyroSolid(const fvMesh &mesh);
  ~pyroSolid();

  void evolve();

  const volScalarField &porosity();

  bool isGasTransferSpecie(const word &name);

  const volScalarField &RR(const word &name);

  const scalar &poreSize();

  const volTensorField& Kp();
  const volScalarField& HTC();
  const volScalarField& T();

  volScalarField mdot();
};

} // namespace pyrolisis
} // namespace Foam
