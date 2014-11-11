/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::surfaceChargeCoupledFvPatchScalarField

Description
    Mixed boundary condition for electric potential, to be used for electrostatics
    across surfaces with charges.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
    - refGradient = sigma / myKValue()
    - refValue = neighbour value
    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    where KDelta is relative permittivity coefficient K * deltaCoeffs

    Example usage:
        myInterfacePatchName
        {
            type        surfaceChargeCoupledFvPatchScalarField;
            neighbourFieldName  phiV;
            myEpr                3.9;
            nbrEpr              80.0;
        }

    Needs to be on underlying mapped(Wall)FvPatch.

    /--------------------------------------------------------------------------\

    Here we modify compressible::turbulentTemperatureCoupledBaffleMixed boundary
    condition, therefore the features and notations inherits. 

    Note: K : heat conduction at patch. Gets supplied how to lookup/calculate K:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'basicThermo' : use basicThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use basicSolidThermo K()
    - 'directionalSolidThermo' directionalK()

    \--------------------------------------------------------------------------/

    Note: runs in parallel with arbitrary decomposition. Uses mapped
    functionality to calculate exchange.

SourceFiles
    surfaceChargeCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceChargeCoupledFvPatchScalarField_H
#define surfaceChargeCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class surfaceChargeCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceChargeCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{

    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;
        const scalar surfCharge_;
        const scalar myEpr_;
	const scalar nbrEpr_;

public:

    //- Runtime type information
    TypeName("surfaceChargeCoupled");


    // Constructors

        //- Construct from patch and internal field
        surfaceChargeCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceChargeCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceChargeCoupledFvPatchScalarField onto a
        //  new patch
        surfaceChargeCoupledFvPatchScalarField
        (
            const surfaceChargeCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceChargeCoupledFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        surfaceChargeCoupledFvPatchScalarField
        (
            const surfaceChargeCoupledFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceChargeCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
