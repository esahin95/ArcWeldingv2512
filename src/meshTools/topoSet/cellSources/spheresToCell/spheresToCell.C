/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "spheresToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(spheresToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, spheresToCell, word);
    addToRunTimeSelectionTable(topoSetSource, spheresToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, spheresToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, spheresToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        spheresToCell,
        word,
        spheres
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        spheresToCell,
        istream,
        spheres
    );
}


Foam::topoSetSource::addToUsageTable Foam::spheresToCell::usage_
(
    spheresToCell::typeName,
    "\n    Usage: spheresToCell List of (centreX centreY centreZ) radius\n\n"
    "    Select all cells with cellCentre within bounding spheres\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::spheresToCell::combine(topoSet& set, const bool add) const
{
    const tmp<pointField> tctrs(this->transform(mesh_.cellCentres()));
    const pointField& ctrs = tctrs();

    /*
    const scalar orad2 = sqr(radius_);
    const scalar irad2 = innerRadius_ > 0 ? sqr(innerRadius_) : -1;

    // Treat innerRadius == 0 like unspecified innerRadius (always accept)

    forAll(ctrs, elemi)
    {
        const scalar d2 = magSqr(ctrs[elemi] - origin_);

        if ((d2 < orad2) && (d2 > irad2))
        {
            addOrDelete(set, elemi, add);
        }
    }
    */

    forAll(ctrs, cellI)
    {
        bool centreInSphere = false;

        forAll(centres_, sphereI)
        {
            scalar offset = magSqr(centres_[sphereI] - ctrs[cellI]);
            if (offset <= radii2_[sphereI])
            {
                centreInSphere = true;
                break;
            }
        }

        addOrDelete(set, cellI, centreInSphere);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spheresToCell::spheresToCell
(
    const polyMesh& mesh,
    const DynamicField<point>& centres,
    const DynamicField<scalar>& radii2,
    const scalar offset
)
:
    topoSetCellSource(mesh),
    centres_(centres),
    radii2_(radii2),
    offset_(offset)
{}


Foam::spheresToCell::spheresToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    centres_(),
    radii2_(),
    offset_(dict.lookupOrDefault("offset", 0.0))
{
    // Raw list of data
    List<scalar> data(dict.lookup("data"));
    
    // Read dynamic fields
    for(int i = 0; i < data.size(); i += 5)
    {
        centres_.append(point(data[i], data[i+1], data[i+2]));
        radii2_.append(sqr(scalar(data[i+3]) + offset_));
    }
}

/*
Foam::spheresToCell::spheresToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    origin_(checkIs(is)),
    radius_(readScalar(checkIs(is))),
    innerRadius_(0)
{}
*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spheresToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells with centre within any of " 
                << centres_.size() << " spheres" << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Given action = " << action << " not supported" << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
