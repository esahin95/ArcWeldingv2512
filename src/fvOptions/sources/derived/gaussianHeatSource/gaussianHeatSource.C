/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gaussianHeatSource.H"
#include "fvMatrices.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcGrad.H"
#include "constants.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(gaussianHeatSource, 0);
    addToRunTimeSelectionTable(option, gaussianHeatSource, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::gaussianHeatSource::alpha() const
{
    auto talpha = volScalarField::New
    (
        "talpha",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, 0.0)
    );
    talpha.ref() = mesh_.lookupObject<volScalarField>(alphaName_);

    return talpha;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::gaussianHeatSource::gaussianHeatSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    alphaName_(coeffs_.getOrDefault<word>("phase", "alpha.gas")),
    Q_(coeffs_.get<scalar>("Q")),
    r_(coeffs_.get<scalar>("r")),
    pos_(coeffs_.get<vector>("pos")),
    d_(normalised(coeffs_.get<vector>("d"))),
    Qv_
    (
        IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimPower/dimVolume, Zero),
        "zeroGradient"
    )
{
    DebugInfo << "Generating fOptions model" << endl;

    fieldNames_.resize(1, "T");

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::gaussianHeatSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo << "Applying fvOptions to : " << eqn.psi().name() << endl;
    
    // Surface energy flux
    const scalar q(Q_ / constant::mathematical::pi / sqr(r_));
    
    // Correction bounds
    const scalar rhoAvg(0.5 * (gMin(rho.primitiveField()) + gMax(rho.primitiveField())));

    // Volumetric energy source
    const volScalarField& alpha = mesh_.lookupObject<volScalarField>(alphaName_);
    const volScalarField magGradAlpha(mag(fvc::grad(alpha)));
    forAll(Qv_, cellI)
    {
        // Distance to source
        const vector x(mesh_.C()[cellI]);
        const vector d(x - (x&d_)*d_);
        const scalar r(mag(d));

        // Volumetric heat source
        Qv_[cellI] = r <= r_? q * magGradAlpha[cellI] * rho[cellI] / rhoAvg : 0.0;
    }
    Qv_.correctBoundaryConditions();

    scalar totalQ = gWeightedSum(mesh_.V(), Qv_.primitiveField());
    Info << "Total energy absorbed [W]: " << totalQ << endl;

    // Add to energy source
    eqn += Qv_.internalField();
}


// ************************************************************************* //
