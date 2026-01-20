ps.solveEnergy();

volScalarField& he = thermo.he();

fvScalarMatrix EEqn
(
    fvm::ddt(alpha, rho, he) + mvConvection->fvmDiv(phi, he)
  + fvc::ddt(alpha, rho, K) + fvc::div(phi, K)
  + (
        he.name() == "e"
      ? fvc::div
        (
            fvc::absolute(phi/fvc::interpolate(alpha*rho), U),
            p,
            "div(phiv,p)"
        )
      : -dpdt
    )
  - fvm::laplacian(alpha*turbulence->alphaEff(), he)
//      + fvm::Sp(ps.HTC()/thermo.Cpv(), he)
 ==
    alpha * Qdot
//      +  ( ps.HTC() * ( ps.T() - T + he/thermo.Cpv() ) )
  + fvOptions(alpha, rho, he)
);

if (MRF.active())
{
    EEqn += fvc::div(MRF.phi(), p);
}

EEqn.relax();

fvOptions.constrain(EEqn);

EEqn.solve();

fvOptions.correct(he);

thermo.correct();

rhoCv = rho*thermo.Cpv();

fvScalarMatrix TEqn
(
    fvm::ddt(alpha,rhoCv, T)
  - fvc::ddt(alpha,rhoCv, T)
  + fvm::Sp(ps.HTC(), T)
  ==
    ps.HTC()*ps.T()
);

TEqn.relax();
TEqn.solve();

// Update enthalpy and thermo
he = thermo.he(p, T);
thermo.correct();

Info<< "min/max(T) = "
    << min(T).value() << ", " << max(T).value() << endl;
