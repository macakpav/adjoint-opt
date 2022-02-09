/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "objectiveTargetForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveTargetForce, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveTargetForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveTargetForce::objectiveTargetForce
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    forcePatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    ),
    forceDirection_(Zero), //PM set later
    rhoInf_(dict.get<scalar>("rhoInf")),
    target_(0), //PM read target
    stressXPtr_
    (
        Foam::createZeroFieldPtr<vector>
        (
            mesh_, "stressX", dimLength/sqr(dimTime)
        )
    ),
    stressYPtr_
    (
        Foam::createZeroFieldPtr<vector>
        (
            mesh_, "stressY", dimLength/sqr(dimTime)
        )
    ),
    stressZPtr_
    (
        Foam::createZeroFieldPtr<vector>
        (
            mesh_, "stressZ", dimLength/sqr(dimTime)
        )
    )
{
    // Sanity check and print info
    if (forcePatches_.empty())
    {
        FatalErrorInFunction
            << "No valid patch name on which to minimize " << type() << endl
            << exit(FatalError);
    }
    if (forcePatches_.size()>1)
    {
        FatalErrorInFunction
            << type() << " works for only one patch." << endl
            << exit(FatalError);
    }
    if (debug)
    {
        Info<< "Minimizing " << type() << " in patches:" << endl;
        for (const label patchI : forcePatches_)
        {
            Info<< "\t " << mesh_.boundary()[patchI].name() << endl;
        }
    }

    // Allocate boundary field pointers
    bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdxdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdStressPtr_.reset(createZeroBoundaryPtr<tensor>(mesh_));
    
    // PM interpret target force
    vector targetForce = dict.get<vector>("targetForce");
    forceDirection_ = targetForce / mag(targetForce);
    target_ = targetForce & forceDirection_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveTargetForce::J()
{

    force_ = calcForce(); //PM wrapped force calc into func

    J_ = sqr(force_ - target_); //PM reformulate J

    return J_;
}


void objectiveTargetForce::update_boundarydJdp()
{
    for (const label patchI : forcePatches_)
    {
        bdJdpPtr_()[patchI] = 2*rhoInf_ *(calcForce()-target_) //PM sqr() derivative
           * forceDirection_; //PM remove denom
    }
}


void objectiveTargetForce::update_dSdbMultiplier()
{
    // Compute contributions with mean fields, if present
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();
    const autoPtr<incompressible::RASModelVariables>&
        turbVars = vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

    tmp<volSymmTensorField> tdevReff = turbVars->devReff(lamTransp, U);
    const volSymmTensorField& devReff = tdevReff();

    for (const label patchI : forcePatches_)
    {
        bdSdbMultPtr_()[patchI] =  2*rhoInf_ *(calcForce()-target_) //PM sqr() derivative
        *(
            (
                forceDirection_& devReff.boundaryField()[patchI]
            )
          + (forceDirection_)*p.boundaryField()[patchI]
        )
       ;//PM remove denom
    }
}


void objectiveTargetForce::update_dxdbMultiplier()
{
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();

    const autoPtr<incompressible::RASModelVariables>&
        turbVars = vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

    //tmp<volSymmTensorField> tdevReff = turbVars->devReff(lamTransp, U);
    //const volSymmTensorField& devReff = tdevReff();

    volScalarField nuEff(lamTransp.nu() + turbVars->nutRef());
    volTensorField gradU(fvc::grad(U));
    volTensorField::Boundary& gradUbf = gradU.boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> nf = patch.nf();
            gradUbf[patchI] = nf*U.boundaryField()[patchI].snGrad();
        }
    }

    volTensorField stress(nuEff*(gradU + T(gradU)));

    stressXPtr_().replace(0, stress.component(0));
    stressXPtr_().replace(1, stress.component(1));
    stressXPtr_().replace(2, stress.component(2));

    stressYPtr_().replace(0, stress.component(3));
    stressYPtr_().replace(1, stress.component(4));
    stressYPtr_().replace(2, stress.component(5));

    stressZPtr_().replace(0, stress.component(6));
    stressZPtr_().replace(1, stress.component(7));
    stressZPtr_().replace(2, stress.component(8));

    volTensorField gradStressX(fvc::grad(stressXPtr_()));
    volTensorField gradStressY(fvc::grad(stressYPtr_()));
    volTensorField gradStressZ(fvc::grad(stressZPtr_()));

    // the notorious second-order derivative at the wall. Use with caution!
    volVectorField gradp(fvc::grad(p));

    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        bdxdbMultPtr_()[patchI] =  2*rhoInf_ *(calcForce()-target_) //PM sqr() derivative
        *(
            (
                (
                   -(forceDirection_.x() * gradStressX.boundaryField()[patchI])
                   -(forceDirection_.y() * gradStressY.boundaryField()[patchI])
                   -(forceDirection_.z() * gradStressZ.boundaryField()[patchI])
                ) & nf
            )
            + (forceDirection_ & nf)*gradp.boundaryField()[patchI]
        ); //PM remove denom
    }
}


void objectiveTargetForce::update_dJdStressMultiplier()
{
    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        bdJdStressPtr_()[patchI] =  2*(calcForce()-target_) //PM sqr() derivative
        *(forceDirection_ * nf); //PM remove denom
    }
}

scalar objectiveTargetForce::calcForce() const
{
    
    vector pressureForce(Zero);
    vector viscousForce(Zero);
    vector cumulativeForce(Zero);


    const volScalarField& p = vars_.pInst();
    const autoPtr<incompressible::turbulenceModel>&
       turbulence = vars_.turbulence();

    volSymmTensorField devReff(turbulence->devReff());

    for (const label patchI : forcePatches_)
    {
        pressureForce += gSum
        (
            mesh_.Sf().boundaryField()[patchI] * p.boundaryField()[patchI]
        );
        // Viscous term calculated using the full tensor derivative
        viscousForce +=  gSum
        (
            devReff.boundaryField()[patchI]
          & mesh_.Sf().boundaryField()[patchI]
        );
    }

    cumulativeForce = pressureForce + viscousForce;

    return rhoInf_ * cumulativeForce & forceDirection_; //PM added rhoInf_
}


const vector& objectiveTargetForce::forceDirection() const
{
    return forceDirection_;
}


void objectiveTargetForce::addHeaderColumns() const
{
    objFunctionFilePtr_()
        << " force " << " targetForce ";
        //<< setw(width_) << "force" << "targetForce ";
        //PM setw not in namespace?

}


void objectiveTargetForce::addColumnValues() const
{
    objFunctionFilePtr_()
        << " " << force_ << " " << target_ << " ";
        //<< setw(width_) << force_ << " " << target_ << " ";
        //PM setw not in namespace?
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
