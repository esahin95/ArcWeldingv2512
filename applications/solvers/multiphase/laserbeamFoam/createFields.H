#include "createRDeltaT.H"

const dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);

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

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Creating multiphaseSystem\n" << endl;
autoPtr<multiphaseInter::multiphaseSystem> fluidPtr =
    multiphaseInter::multiphaseSystem::New(mesh);

multiphaseInter::multiphaseSystem& fluid = fluidPtr();

if (!fluid.incompressible())
{
    FatalError << "One or more phases are not incompressible. " << nl
        << "This is a incompressible solver." << abort(FatalError);
}

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    fluid.rho()
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

const word interfaceTrackingScheme
(
    fluid.lookup("interfaceTrackingScheme")
);

if
(
    interfaceTrackingScheme != "MULES"
 && interfaceTrackingScheme != "isoAdvector"
)
{
    FatalErrorIn("laserbeamFoam.C")
        << "Unknown interface tracking scheme "
        << interfaceTrackingScheme
        << ". Available options are 'MULES' and 'isoAdvector'"
        << exit(FatalError);
}

IOdictionary phaseChangeProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading phase change properties for metal" << endl;
dictionary& dict = phaseChangeProperties.subDict("metal");

const dimensionedScalar Tsol
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    dict
);

const dimensionedScalar Tliq
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    dict
);

const dimensionedScalar Lmelt
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    dict
);

Info<< "Reading phase change properties for metal" << endl;
dict = phaseChangeProperties.subDict("gas");

const dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    dict
);

const dimensionedScalar Tvap
(
    "Tvap",
    dimensionSet(0, 0, 0, 1, 0),
    dict
);

const dimensionedScalar Mm
(
    "Mm",
    dimensionSet(1, 0, 0, 0, -1),
    dict
);

const dimensionedScalar Lvap
(
    "Lvap",
    dimensionSet(0, 2, -2, 0, 0),
    dict
);

const dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    scalar(8.314)
);

const dimensionedScalar SB
(
    "SB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    scalar(5.67e-8)
);

const dimensionedScalar T0 
(
    "T0",
    dimensionSet(dimTemperature),
    dict
);

volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DC", dimensionSet(1,-3,-1,0,0), 1.0e14)
);

volScalarField alphaLiquid
(
    IOobject
    (
        "alphaLiquid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(min((T - TSolidus)/(TLiquidus - TSolidus), 1.0), 0.0)
);

volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

volScalarField Tcorr
(
    IOobject
    (
        "Tcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tcorr", dimensionSet(0, 0, 0, 1, 0), 0.0)
);

// Recoil Pressure
volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

// Evaporation Cooling
volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);

surfaceScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho) * (2.0 / (rho1 + rho2))
);

surfaceScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    damper
);


volScalarField e1temp
(
    IOobject
    (
        "e1temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("e1temp", dimensionSet(0, 0, 0, -0, 0), 0.0)
);
e1temp = fvc::average(alphaLiquid);

/*
dimensionedScalar elec_resistivity_metal
(
    "elec_resistivity",
    dimensionSet(1, 3, -3, 0, 0),
    transportPropertiesMetal//.lookup("elec_resistivity")
);
*/

// Electrical resistivity of condensate only
/*
volScalarField electrical_resistivity
(
    IOobject
    (
        "electrical_resistivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    elec_resistivity_metal
);
*/

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho * fluid.Cp()
);
rhoCp.oldTime();

surfaceScalarField rhoPhiCp
(
    IOobject
    (
        "rhoPhiCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(fluid.Cp())*rhoPhi
);

const dimensionedScalar DClarge
(
    "DClarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DCsmall
(
    "DCsmall", dimless, 1.0e-12
);

// Construct incompressible turbulence model
autoPtr
<
    CompressibleTurbulenceModel<multiphaseInter::multiphaseSystem>
> turbulence
(
    CompressibleTurbulenceModel<multiphaseInter::multiphaseSystem>::New
    (
        rho,
        U,
        rhoPhi,
        fluid
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// Create laser heat source
//laserHeatSource laser(mesh);

// depositionTemp=laser.deposition();
// depositionTemp.correctBoundaryConditions();

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"


#include "isoAdvector/createPorosity.H"
isoAdvection advector(alpha1, phi, U);