/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::objectives::objectiveTargetForce

Description
    Optimizes for target force. Works for only one patch. 
    Original file objectiveForce.C

    targetForce
    {
        weight          1.;
        type            TargetForce;
        patches         (PROF);
        targetForce     (0 6.5 0);
        rhoInf          1.225;
    }

SourceFiles
    objectiveTargetForce.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveTargetForce_H
#define objectiveTargetForce_H

#include "objectiveIncompressible.H"
#include "wallFvPatch.H"
#include "createZeroField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                        Class objectiveTargetForce Declaration
\*---------------------------------------------------------------------------*/

class objectiveTargetForce
:
    public objectiveIncompressible
{
protected:

    // Protected Data

        labelHashSet forcePatches_;
        vector forceDirection_;
        // scalar Aref_; //PM remove
        scalar rhoInf_;
        // scalar UInf_; //PM remove
        scalar target_; //PM add target_
        scalar force_; //PM add force_

        autoPtr<volVectorField> stressXPtr_;
        autoPtr<volVectorField> stressYPtr_;
        autoPtr<volVectorField> stressZPtr_;


public:

    //- Runtime type information
    TypeName("TargetForce");


    // Constructors

        //- Construct from components
        objectiveTargetForce
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveTargetForce() = default;


    // Member Functions

        //- Return the objective function value
        scalar J();

        //- Update values to be added to the adjoint wall velocity
        void update_boundarydJdp();

        //- Update delta(n dS)/delta b multiplier
        void update_dSdbMultiplier();

        //- Update delta(x)/delta b multiplier
        void update_dxdbMultiplier();

        //- Update dJ/dStress multiplier
        void update_dJdStressMultiplier();

        //- PM added calcForce funcs
        scalar calcForce() const;
        scalar calcForce(label patchI) const;

        //- Return force direction
        const vector& forceDirection() const;

        //PM Helper write functions

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
