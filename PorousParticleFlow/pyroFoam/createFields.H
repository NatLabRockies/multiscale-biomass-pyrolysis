#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

Info << "Initialize pyrolysis\n";
pyrolisis::pyroSolid ps(mesh);

const volScalarField& alpha = ps.porosity();
surfaceScalarField alphaf = fvc::interpolate(alpha);

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();


#include "compressibleCreatePhi.H"

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<rhoReactionThermo>> reaction
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);
Info<< "Reading g\n" << endl;
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"


volScalarField RRQdot
(
    IOobject
    (
        "RRQdot",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("RRQdot",dimEnergy/dimTime/dimVolume,0.)
);

volScalarField mdot
(
    IOobject
    (
        "mdot",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mdot",dimDensity/dimTime,0.)
);
dimensionedScalar pore_size("ps",dimLength,ps.poreSize());

volTensorField Kp
(
    IOobject
    (
        "Kp",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("K",dimViscosity*dimDensity/dimArea,tensor::zero)
);

scalar globalContErr = 0;

volScalarField rhoCv(rho*thermo.Cpv());
