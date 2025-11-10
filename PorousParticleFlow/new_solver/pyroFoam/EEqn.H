{
    //!NOTE: This should be changed!

    volScalarField& he = thermo.he();


    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + fvm::Sp(ps.HTC()/thermo.Cpv(), he)
      + thermophysicalTransport->divq(he)
      + fvm::SuSp(RRQdot/he,he)
     ==
        reaction->Qdot()
      +  ( ps.HTC() * ( ps.T() - T + he/thermo.Cpv() ) )
      + fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
