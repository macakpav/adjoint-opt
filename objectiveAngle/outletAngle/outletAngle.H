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


Class
    Foam::objectives::outletAngle

Description
    Optimizes to reach target Uy at outletPatch. 
    Works only for planar patch perpendicular to x direction.
    angle
    {
        weight          1.;
        type            OutletAngle;
        outletPatch     "outlet";
        UyTarget        50.5;
    }

SourceFiles
    outletAngle.C

\*---------------------------------------------------------------------------*/

#ifndef outletAngle_H
#define outletAngle_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                      Class outletAngle Declaration
\*---------------------------------------------------------------------------*/

class outletAngle
:
    public objectiveIncompressible
{
    // Private data

        //- Outlet patch
        label outletPatch_;

        //- Mean Uy velocity
        scalar avgUy_;

        //- Target velocity
        scalar UyTarget_;


public:

    //- Runtime type information
    TypeName("OutletAngle");


    // Constructors

        //- from components
        outletAngle
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~outletAngle() = default;


    // Member Functions

        //- Return the objectiveReportPatches
        void initialize();

        //- Return the objective function value
        scalar J();

        //- Update values to be added to the adjoint inlet velocity
        void update_boundarydJdp();

        //- Update values to be added to the adjoint outlet velocity
        void update_boundarydJdv();

        //- Update values to be added to the adjoint outlet pressure
        void update_boundarydJdvn();

        //- Update values to be added to the adjoint outlet tangential velocity
        void update_boundarydJdvt();

        // Helper write functions

            //- Write headers for additional columns
            virtual void addHeaderColumns() const;

            //- Write information to additional columns
            virtual void addColumnValues() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
