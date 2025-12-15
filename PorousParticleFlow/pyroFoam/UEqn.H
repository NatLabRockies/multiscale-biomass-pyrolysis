    MRF.correctBoundaryVelocity(U);

    Kp = tensor::I
        * turbulence->muEff() * ( 1.0 - alpha )
        * (
            (150. * (alpha) )
            /
            ( sqr(alpha) * alpha * sqr(pore_size) )
        );

    /*
       forAll(Kp, cellI)
       {
            Kp[cellI] = tensor::I;
            Kp[cellI].xx() = mu[cellI]/   sqr(pore_size);
            Kp[cellI].yy() = ...
            Kp[cellI].zz() = 0;
       }
    */

    Kp.correctBoundaryConditions();


    fvVectorMatrix UEqn
    (
        fvm::ddt(alpha,rho, U) + fvm::div(phi, U,"div(phi,U)")
      + turbulence->divDevRhoReff(U)
      + fvm::Sp(tr(Kp), U) + ( (Kp - (tensor::I*tr(Kp)))&U )
     ==
        fvOptions(alpha, rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
