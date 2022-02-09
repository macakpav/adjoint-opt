/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "outletAngle.H"
#include "createZeroField.H"
#include "coupledFvPatch.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(outletAngle, 1);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    outletAngle,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

outletAngle::outletAngle
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    outletPatch_(
        mesh.boundaryMesh().findPatchID
        (
            dict.get<word>("outletPatch")
        )
    ),
    avgUy_(0),
    UyTarget_(dict.get<scalar>("UyTarget"))
{
    if (outletPatch_==-1) {
        FatalErrorInFunction
            << "No valid patch name on which to minimize " << type() << endl
            << exit(FatalError);
    }

    // Allocate boundary field pointers
    bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    bdJdvtPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar outletAngle::J()
{
    J_ = Zero;

    // References
    const volVectorField& U = vars_.UInst();
    const vector j(0, 1, 0);
    const scalarField& phiPatch = vars_.phiInst().boundaryField()[outletPatch_];
    const scalarField& magSf = mesh_.boundary()[outletPatch_].magSf(); //face area
    const vectorField& Ub = U.boundaryField()[outletPatch_];

    J_ = gSum
    (
        sqr( (Ub & j) - UyTarget_ ) * magSf
    );
    avgUy_ = gSum((Ub & j) * phiPatch) / gSum(phiPatch) ;
    return J_;
}

void outletAngle::update_boundarydJdp() // probably unnecessary if bdJdpPtr_ not initialized in constructor
{
    bdJdpPtr_()[outletPatch_] = vector(Zero);
}

void outletAngle::update_boundarydJdv()
{
    const volVectorField& U = vars_.U();
    const vector j(0, 1, 0);
    const fvPatchVectorField& Ub = U.boundaryField()[outletPatch_];

    bdJdvPtr_()[outletPatch_] = 2*j*( (Ub&j) - UyTarget_);

}

void outletAngle::update_boundarydJdvn() // probably unnecessary if bdJdvnPtr_ not initialized in constructor
{
    bdJdvnPtr_()[outletPatch_] = 0.0;
}

void outletAngle::update_boundarydJdvt()
{
    const volVectorField& U = vars_.U();
    const vector j(0, 1, 0);
    const fvPatchVectorField& Ub = U.boundaryField()[outletPatch_];

    bdJdvtPtr_()[outletPatch_] = 2*j*( (Ub&j) - UyTarget_);
}


void outletAngle::addHeaderColumns() const
{
    objFunctionFilePtr_()
        << setw(width_) << "avgUy " << "UyTarget ";

}


void outletAngle::addColumnValues() const
{
    objFunctionFilePtr_()
        << setw(width_) << avgUy_ << " " << UyTarget_ << " ";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
