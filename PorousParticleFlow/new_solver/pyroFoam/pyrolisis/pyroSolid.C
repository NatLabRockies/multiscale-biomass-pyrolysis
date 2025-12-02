#include "pyroSolid.h"
#include "irreversibleArrheniusReaction.h"

using namespace Foam;
using namespace pyrolisis;

pyroSolid::pyroSolid(const fvMesh& mesh)
:
    m_mesh(mesh),
    m_dict (
        IOobject
        (
            "pyrolisisDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    m_porosity
    (
        IOobject
        (
            "porosity",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh
    ),
    m_rhoField
    (
        IOobject
        (
            "rho.solid",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh,
        dimensionedScalar("rho",dimDensity,0.)
    ),
    m_T
    (
        IOobject
        (
            "T.solid",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh
    ),
    m_rhoCp
    (
        IOobject
        (
            "rhoCp.solid",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        m_mesh,
        dimensionedScalar("rhoCp",dimEnergy/dimTemperature/dimVolume,0.)
    ),
    m_Qdot
    (
        IOobject
        (
            "Qdot.solid",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh,
        dimensionedScalar("Qdot",dimEnergy/dimTime/dimVolume,0.)
    ),
    // m_Kp
    // (
    //     IOobject
    //     (
    //         "Kp",
    //         m_mesh.time().timeName(),
    //         m_mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     m_mesh,
    //     dimensionedTensor("Kp",dimViscosity*dimDensity/dimArea,tensor::zero)
    // ),
    m_htc
    (
        IOobject
        (
            "htc",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh,
        dimensionedScalar("htc",dimPower/dimVolume/dimTemperature,0.)
    ),
    m_kappa_tot
    (
        IOobject
        (
            "kappa.solid",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh,
        dimensionedScalar("kappa",dimPower/dimLength/dimTemperature,0.)
    ),
    m_nSubTimeSteps(m_dict.lookupOrDefault<label>("nSubTimeSteps",1)),
    m_poreSize(readScalar(m_dict.lookup("poreSize"))),
    m_speciesName(m_dict.lookup("species")),
    m_species(m_speciesName.size()),
    m_rho(m_speciesName.size()),
    m_cp(m_speciesName.size()),
    m_molWeight(m_speciesName.size()),
    m_kappa(m_speciesName.size()),
    m_is_gas(m_speciesName.size()),
    m_addToPoro(m_speciesName.size(),true),
    m_formationEnthalpy(m_speciesName.size())
{
    Info << "Reading pyrolisis model" << endl;
    Info << "Found species: " << m_speciesName << endl;

    forAll(m_speciesName, specieI)
    {

        dictionary& speciesDict(m_dict.subDict("speciesCoeffs").subDict(m_speciesName[specieI]));
        m_species.set (
            specieI,
            new volScalarField
            (
                IOobject
                (
                    m_speciesName[specieI] + ".solid",
                    m_mesh.time().timeName(),
                    m_mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                m_mesh,
                dimensionedScalar("Y",dimless,0.)
            )
        );

        m_species[specieI].storeOldTime();

        m_rho[specieI] = scalar(readScalar(speciesDict.lookup("rho")));
        m_cp[specieI] = scalar(readScalar(speciesDict.lookup("cp")));
        m_kappa[specieI] = scalar(readScalar(speciesDict.lookup("kappa")));
        m_molWeight[specieI] = scalar(readScalar(speciesDict.lookup("molWeight")));
        m_molWeight[specieI] *= 1e-3;
        m_is_gas[specieI].first = speciesDict.lookupOrDefault<bool>("gas",false);

        m_addToPoro[specieI] = speciesDict.lookupOrDefault<bool>("addToPorosity",true);
        m_formationEnthalpy[specieI] = speciesDict.lookupOrDefault<scalar>("hf",0.); // In energy per mass

        if (m_is_gas[specieI].first)
        {
            m_is_gas[specieI].second = m_reactionRates.size();

            m_reactionRates.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        m_speciesName[specieI]  + ".RR",
                        m_mesh.time().timeName(),
                        m_mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    m_mesh,
                    dimensionedScalar("Rdot",dimDensity/dimTime,0.)
                )
            );
        }

        forAll(m_rhoCp,cellI)
        {
            m_rhoCp[cellI] += m_species[specieI][cellI]*m_rho[specieI]*m_cp[specieI];
        }

    }

    m_rhoCp.correctBoundaryConditions();
    m_rhoCp.storeOldTime();


    List<word> reactionList(m_dict.lookup("reactions"));

    Info << "Found reactions: " << reactionList << endl;

    m_reactions.resize(reactionList.size());

    forAll(reactionList, reactI)
    {
        if( !m_dict.subDict("reactionCoeffs").found(reactionList[reactI]) )
        {
            FatalErrorInFunction << "Cannot find reaction " << reactionList[reactI] << "\n" << abort(FatalError);
        }

        m_reactions.set (
            reactI,
            new irreversibleArrheniusReaction
            (
                reactionList[reactI],
                m_dict.subDict("reactionCoeffs").subDict(reactionList[reactI]),
                m_speciesName
            )
        );
    }
}

pyroSolid::~pyroSolid()
{}

void pyroSolid::evolve()
{
//    updateKp();

    solveEnergy();

    Info << "Updating solid composition" << endl;

    scalarField Y(m_species.size(),0.);
    //scalarField Y0(m_species.size(),0.);
    /* Reset reaction rate */
    forAll(m_species, specieI)
    {
        if( m_is_gas[specieI].first )
        {
            m_reactionRates[m_is_gas[specieI].second] *= 0.;
            m_reactionRates[m_is_gas[specieI].second].correctBoundaryConditions();
        }
    }


    /* Go cell-by-cell */
    forAll(m_mesh.C(), cellI)
    {
        // Skip if outside the solid
        if ( m_porosity[cellI] > 0.999 )
        {
            continue;
        }

        m_Qdot[cellI] = 0.; // Reset heat of reaction

        scalar dt = m_mesh.time().deltaT().value();

        scalar sub_dt = dt / m_nSubTimeSteps;

        forAll(m_species, specieI)
        {
            //- Convert from [vol_specie/m3] to [mol_specie/m3] by multiplying by density and dividing
            //  by molar weight.
            //  Always use oldTime to make it work with pimple outer iterations.
            Y[specieI] =  m_species[specieI][cellI] * ( m_rho[specieI] / m_molWeight[specieI] );
            //Y0[specieI] = Y[specieI];
        }
        scalarField ntot(m_species.size(),0.);

        /* Use Euler marching */
        for (label ts = 0; ts < m_nSubTimeSteps; ts++)
        {
            scalarField ndot(m_species.size(),0.);

            forAll(m_reactions, reactI)
            {
                // This updates ndot and adds the heat of reaction for m_reactions[reactI] to
                // m_Qdot. This is added every sub time step.
                m_Qdot[cellI] += m_reactions[reactI].computeSources(m_T[cellI], Y, m_molWeight, ndot) / scalar(m_nSubTimeSteps);
            }

            Y += sub_dt * ndot;

            Y = max(Y,0.);
            ntot += ndot / scalar(m_nSubTimeSteps);
        }

        /* Update reaction rate and concentration*/
        forAll(m_species, specieI)
        {
            //scalar Y_0 = m_species[specieI].oldTime()[cellI] * m_rho[specieI];
            m_species[specieI][cellI] = Y[specieI] * m_molWeight[specieI] / m_rho[specieI];
            //Info << "is gas " << m_species[specieI].name() << ": " << m_is_gas[specieI].first << endl;
            if( m_is_gas[specieI].first )
            {
                m_reactionRates[m_is_gas[specieI].second][cellI] = ntot[specieI] *  m_molWeight[specieI];
                //Info << "ndot " << m_species[specieI].name() << ": " << ntot[specieI] << endl;
                m_species[specieI][cellI] = 0.;
            }
        }

        /* Update local porosity and density */
        scalar vol_species(0.);
        scalar mass_species(0.);
        forAll(m_species, specieI)
        {
            if(!m_addToPoro[specieI]) continue;

            vol_species += m_species[specieI][cellI];
            mass_species += m_species[specieI][cellI] * m_rho[specieI];
        }

        m_porosity[cellI] = 1.0 - vol_species; //1.0 - (1.0 - m_porosity[cellI])*mass_species/m_rhoField[cellI];
        m_rhoField[cellI] = mass_species;
    }

    m_porosity.correctBoundaryConditions();
    m_rhoField.correctBoundaryConditions();

}

const volScalarField& pyroSolid::porosity()
{
    return m_porosity;
}

label pyroSolid::getSpecieId(const word& name)
{
    label id = -1;
    forAll(m_speciesName, specieI)
    {
        if( m_speciesName[specieI] == name )
        {
            id = specieI;
        }
    }

    return id;
}


bool pyroSolid::isGasTransferSpecie(const word& name)
{
    label id = getSpecieId(name);

    if (id == -1)
    {
        FatalErrorInFunction << "Error: unknown gas specie " << name << "\n" << abort(FatalError);
    }

    return m_is_gas[id].first;
}

const volScalarField& pyroSolid::RR(const word& name)
{
    label id = getSpecieId(name);

    if (id == -1)
    {
        FatalErrorInFunction << "Error: unknown gas specie " << name << "\n" << abort(FatalError);
    }

    Info << "Gas specie " << name  << " id: " << m_is_gas[id].second << " transfer mass: " << fvc::domainIntegrate(m_reactionRates[m_is_gas[id].second]) << endl;
    return m_reactionRates[m_is_gas[id].second];
}

const scalar& pyroSolid::poreSize()
{
    return m_poreSize;
}

// void pyroSolid::updateKp()
// {
//     //- Basic Darcy for now

//     const auto& mu = m_mesh.lookupObject<volScalarField>("thermo:mu");

//     forAll(m_Kp, cellI)
//     {
//         m_Kp[cellI] =  tensor::I
//             * mu[cellI]
//             * (
//                 (150. * (m_porosity[cellI]) )
//                 /
//                 sqr(m_porosity[cellI]) * m_porosity[cellI] * sqr(m_poreSize)
//             );
//     }

//     m_Kp.correctBoundaryConditions();

// }

void pyroSolid::updateHTC()
{
    const auto& mu = m_mesh.lookupObject<volScalarField>("thermo:mu");
    dimensionedScalar dimPS("dimPS",dimLength,m_poreSize);
    const auto& rho = m_mesh.lookupObject<volScalarField>("rho");
    const auto& U = m_mesh.lookupObject<volVectorField>("U");
    //const volTensorField Kp = m_mesh.lookupObject<volTensorField>("Kp");

    const volScalarField Re = mag(U) * rho * dimPS / mu;
    const scalar Pr = 0.7;
    //volScalarField Nu = sqrt( mag(Kp) / mu ) / dimPS;

    // Withaker correlation
    const volScalarField Nu = (1. - m_porosity) * ( 2. + 1.1 * pow(Re, 0.6) * pow(Pr,1.0/3.0));


    forAll(m_htc, cellI)
    {
        scalar kappa(0.);

        forAll(m_kappa, specieI)
        {
            kappa += m_kappa[specieI]*m_species[specieI][cellI];
        }

        // Add surface area per unit volume to be consistent with equations
        m_htc[cellI] = ( Nu[cellI]/m_poreSize * kappa ) * ( 6.0 * (1.0 - m_porosity[cellI] ) / m_poreSize);
    }

    m_htc.correctBoundaryConditions();

}

// const volTensorField& pyroSolid::Kp()
// {
//     return m_Kp;
// }

const volScalarField& pyroSolid::HTC()
{
    return m_htc;
}

const volScalarField& pyroSolid::T()
{
    return m_T;
}

void pyroSolid::solveEnergy()
{
    //- Compute effective conductivity
    dimensionedScalar kappaDim("kappaDim",dimPower/dimLength/dimTemperature, 1.);
    surfaceScalarField porosityf(fvc::interpolate(m_porosity));
    const auto& T_fluid = m_mesh.lookupObject<volScalarField>("T");

    m_kappa_tot *= 0.;

    forAll(m_kappa, specieI)
    {
        m_kappa_tot += m_species[specieI]*kappaDim*m_kappa[specieI];
    }

    surfaceScalarField kappaf = fvc::interpolate(m_kappa_tot);

    forAll(m_rhoCp,cellI)
    {
        scalar rhocp = 1e-5;
        forAll(m_kappa, specieI)
        {
            rhocp += m_species[specieI][cellI]*m_rho[specieI]*m_cp[specieI];

        }

        m_rhoCp[cellI] = rhocp;
    }
    m_rhoCp.correctBoundaryConditions();

    //- Set to zero at solid boundary faces to avoid diffusion outside solid
    // forAll(kappaf, faceI)
    // {
    //     if ( porosityf[faceI] > 0.5 )
    //     {
    //         kappaf[faceI] *= 0.;
    //     }
    // }

    updateHTC();

    //const volScalarField& RRQdot = m_mesh.lookupObject<volScalarField>("RRQdot");

    // Formation enthaly term
    // volScalarField formH(fvc::ddt(m_rhoCp, m_T)*0.); // Just initialize
    // forAll(m_speciesName, specieI)
    // {
    //     if (m_formationEnthalpy[specieI] < 1e-16) continue;
    //
    //     const volScalarField ddtspecies = fvc::ddt(m_species[specieI]);
    //     const scalar coeff = m_rho[specieI]*m_formationEnthalpy[specieI];
    //
    //     forAll(formH,cellI)
    //     {
    //         formH[cellI] += coeff*ddtspecies[cellI];
    //     }
    // }


    if( !(m_dict.lookupOrDefault<bool>("solveSolidEnergy",true)) ) return;

    fvScalarMatrix TEqn
    (
        fvm::ddt(m_rhoCp, m_T)
      - fvm::laplacian(kappaf, m_T)
      + fvm::Sp(m_htc, m_T)
      ==
        ( m_htc * T_fluid )
      + m_Qdot
//      - RRQdot
    );

    TEqn.relax();

    TEqn.solve();
    m_T.correctBoundaryConditions();

    Info<< "min/max(T) solid = "
        << min(m_T).value() << ", " << max(m_T).value() << endl;

}

volScalarField pyroSolid::mdot()
{

    volScalarField mdot
    (
        IOobject
        (
            "mdot",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        m_mesh,
        dimensionedScalar("mdot",dimDensity/dimTime,0.)
    );

    forAll(m_reactionRates, rateI)
    {
        mdot += m_reactionRates[rateI];
    }

    mdot.correctBoundaryConditions();

    return mdot;
}
