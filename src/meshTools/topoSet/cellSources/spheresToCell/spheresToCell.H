/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::spheresToCell

Description
    A \c topoSetCellSource to select cells based
    on cell centres inside a given bounding sphere.

    Operands:
    \table
      Operand   | Type    | Location
      output    | cellSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        cellSet;
        action      <action>;

        // Mandatory entries
        source      spheresToCell;
        origin      (0 0 0);
        radius      0.5;

        // Optional entries
        innerRadius 0;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: spheresToCell           | word |  yes  | -
      origin     | The origin (centre) of the sphere   | vector | yes | -
      radius     | The (outer) radius of sphere        | scalar | yes | -
      innerRadius | The inner radius of sphere         | scalar | no  | 0
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    spheresToCell.C

\*---------------------------------------------------------------------------*/

#ifndef spheresToCell_H
#define spheresToCell_H

#include "topoSetCellSource.H"

#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class spheresToCell Declaration
\*---------------------------------------------------------------------------*/

class spheresToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Centres
        DynamicField<point> centres_;

        //- Squared Radii
        DynamicField<scalar> radii2_;

        //- Offset
        scalar offset_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("spheresToCell");


    // Constructors

        //- Construct from components
        spheresToCell
        (
            const polyMesh& mesh,
            const DynamicField<point>& centres,
            const DynamicField<scalar>& radii2,
            const scalar offset
        );

        //- Construct from dictionary
        spheresToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        //spheresToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~spheresToCell() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
