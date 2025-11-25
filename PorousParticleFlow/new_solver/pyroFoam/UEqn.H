// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

Kp = tensor::I 
    * mu * ( 1.0 - poro )
    * (
        (150. * (poro) )
        /
        ( sqr(poro) * poro * sqr(pore_size) ) 
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

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
  + fvm::Sp(tr(Kp), U) + ( (Kp - (tensor::I*tr(Kp)))&U )
 ==
    fvModels.source(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
