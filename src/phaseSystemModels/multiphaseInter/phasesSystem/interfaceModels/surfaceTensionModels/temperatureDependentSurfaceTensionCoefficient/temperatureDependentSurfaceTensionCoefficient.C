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

#include "temperatureDependentSurfaceTensionCoefficient.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseInter
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(temperatureDependentSurfaceTensionCoefficient, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        temperatureDependentSurfaceTensionCoefficient,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInter::surfaceTensionModels::temperatureDependentSurfaceTensionCoefficient::
temperatureDependentSurfaceTensionCoefficient
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    surfaceTensionModel(dict, pair, registerObject),
    TName_(dict.getOrDefault<word>("T", "T")),
    sigma_(Function1<scalar>::New("sigma", dict))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::multiphaseInter::surfaceTensionModels::temperatureDependentSurfaceTensionCoefficient::
sigma() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();

    auto tsigma = volScalarField::New
    (
        "sigma",
        IOobject::NO_REGISTER,
        mesh,
        dimensionSet(1,0,-2,0,0)
    );
    auto& sigma = tsigma.ref();

    const volScalarField& T = mesh.lookupObject<volScalarField>(TName_);

    sigma.field() = sigma_->value(T.field());

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaBf, patchi)
    {
        sigmaBf[patchi] = sigma_->value(TBf[patchi]);
    }

    return tsigma;
}


// ************************************************************************* //
