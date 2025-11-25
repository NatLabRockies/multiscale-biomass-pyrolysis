tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

RRQdot *= 0.;

reaction->correct();

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport->divj(Yi)
         ==
            reaction->R(Yi)
          + fvModels.source(rho, Yi)
        );

        if( ps.isGasTransferSpecie(Yi.name()) )
        {
            YiEqn -= ps.RR(Yi.name());
            RRQdot -= ps.RR(Yi.name())*dimensionedScalar("Hf",dimEnergy/dimMass, composition.Hf(i) );
        }

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();
