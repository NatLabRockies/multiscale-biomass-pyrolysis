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

mdot *= 0.;

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(alpha, rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(alpha * turbulence->muEff(), Yi)
             ==
                alpha*reaction->R(Yi)
              + fvOptions(rho, Yi)
            );


            if( ps.isGasTransferSpecie(Yi.name()) )
            {
                YiEqn -= ps.RR(Yi.name());
                mdot -= ps.RR(Yi.name());
                RRQdot -= ps.RR(Yi.name())*dimensionedScalar("Hf",dimEnergy/dimMass, composition.Hc(i) );
            }

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);
}
