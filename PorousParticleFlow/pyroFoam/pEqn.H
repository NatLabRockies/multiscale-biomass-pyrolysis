rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(alpha*rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(alpha*rho*HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
         fvc::ddt(alpha, rho) + fvc::div(phiHbyA) - mdot
       - fvm::laplacian(rhorAUf, p_rgh)
    );

    p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        alphaPhi = phiHbyA + p_rghEqn.flux();
        phi = alphaPhi/alphaf;

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
}

p = p_rgh + rho*gh;

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"
