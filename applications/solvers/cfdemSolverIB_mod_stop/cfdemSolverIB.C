/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverIB

Description
    Transient solver for incompressible flow.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6, 
    where additional functionality for CFD-DEM coupling using immersed body
    (fictitious domain) method is added.
Contributions
    Alice Hager
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "cfdemCloudIB.H"
#include "implicitCouple.H"

#include "averagingModel.H"
#include "regionModel.H"
#include "voidFractionModel.H"

#include "dynamicFvMesh.H" //dyM

#include "cellSet.H"

#if defined(version22)
    #include "meshToMeshNew.H"
    #include "fvIOoptionList.H"
#endif


// AYCOCK EDITS
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"

    #include "createDynamicFvMesh.H"

    #include "createFields.H"

    #include "initContinuityErrs.H"

    #if defined(version23)
    #include "createFvOptions.H"
    #endif

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloudIB particleCloud(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Initialize a IOdict for communicating the pressure gradient:
    IOdictionary myGradPDict
    (
        IOobject
        (
            "myGradP",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    myGradPDict.add("gradient", scalar(0));


    // ****************************************************************************************************** //
    // OUTLINE OF FD-DEM COUPLING PROCESS (UPDATED 03/2016):
    // 1) Initialize the DEM solver, but do not execute any 6DOF (get particile initial positions and velocities from LIGGGHTS input file)
    // 2) Calculate the fluid forces & moments on the particles using the current velocity field (used in 3)
    // 3) Perform the DEM calculation (includes 6DOF and collisions)
    // 4) While the DEM solver is running, solve for the fluid motion using the PISO algorithm
    //     ^ Are the particles present in this step? Hager says no, but there is an extra term
    //     'particleCloud.ddtVoidfraction' in the pressure - Poisson equation...
    // 5) Update the velocity field from (4) using the particle velocities obtained in (3)
    // 6) Perform a velocity and pressure correction to satisfy continuity (not working?)
    // 7) Repeat the process, starting at (2). 
    // ****************************************************************************************************** //

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //=== dyM ===================
        interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        mesh.update(); //dyM

        #include "readPISOControls.H"
        #include "CourantNo.H"

        // do particle stuff

        // ****************************************************************************************************** //
        // Outline of the steps in function 'particleCloud.evolve()': (UPDATED 03/2016)
        // (look in derived/cfdemCloudIB.C)
        // ****************************************************************************************************** //
        // 1) get DEM data (position, velocity, angular velocity)
        // 2) find cellIDs contained in particle domain (see locate model 'engineSearchIB')
        // 3) set the void fraction field (see void fraction model 'IBVoidFraction')
        // 4) set the forces and moments (via force model(s) chosen in forceDict;
        //    ^ see 'subModels' directory; 'ShirgaonkarIB' for viscous & pressure drag
        //                                 'ShirgaonkarIBTorques' for viscous-induced moment or "torque" on particle
        //                                 'ArchimedesIB' for buoyancy force 
        // 5) Call LIGGGHTS and run for the specified number of sub-timesteps (usually dt_DEM = 0.1 * dt_CFD)
        //    ^ this occurs during the call to 'dataExchangeM().couple()' 
        // 6) write DEM data (in preparation for the next timestep)
        // ****************************************************************************************************** //
        Info << "- evolve()" << endl;
        particleCloud.evolve();
        // ****************************************************************************************************** //

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(voidfraction,U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
                #if defined(version23)
                ==
                fvOptions(U)
                #endif
            );

            UEqn.relax();

            #if defined(version23)
            fvOptions.constrain(UEqn);
            #endif

            if (momentumPredictor)
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rUA = 1.0/UEqn.A();
                surfaceScalarField rUAf(fvc::interpolate(rUA));

                U = rUA*UEqn.H();
                #ifdef version23
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + rUAf*fvc::ddtCorr(U, phi);
                #else
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, U, phi);
                #endif
                adjustPhi(phi, U, p);

                #if defined(version23)
                fvOptions.makeRelative(phi);
                //fvOptions.relativeFlux(phi);
                #endif

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        // ADDITION OF A PRESSURE TERM FROM 'particleCloud.ddtVoidfraction();
                        // What is this extra term?
                        // Outline 'particleCloud.ddtVoidfraction()':
                        // (look in cfdemCloud/cfdemCloud.C; 
                        // ^- how does cfdemCloudIB reference cfdemCloud? forceModel?)
                        // ddtVoidfraction_ = fvc::ddt(voidfraction);
                        // dimensions: (0,0,-1,0,0)
                        fvm::laplacian(rUA, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

        turbulence->correct();


        // Perform velocity correction; I assume this step retrieves the velocity information from 
        // the DEM solver? Check the code for calcVelocityCorrection()
        Info << "particleCloud.calcVelocityCorrection() " << endl;
        volScalarField voidfractionNext=mesh.lookupObject<volScalarField>("voidfractionNext");

        // ****************************************************************************************************** //
        // Outline the steps for calcVelocityCorrection():
        // (look in derived/cfdemCloudIB.C)
        // 1) Calculate particle velocities:
        //     a) calculate 'r' vector
        //     b) get angular velocities from 'angularVelocities()'
        //     c) calculate the cross product of angular velocities and r vector (omega -> U)
        //
        // 2) --> IMPOSE THE PARTICLE VELOCITIES ON THE FLUID VELOCITY FIELD <--
        //
        // 3) Make the velocity field divergence free (also correct pressure field) <-- Where does this occur?? Because it's not working! :)
        // DONE
        particleCloud.calcVelocityCorrection(p,U,phiIB,voidfractionNext);

        // ************************************** //
        // Correct phi:
        phi = linearInterpolate(U) & mesh.Sf();
        // ************************************** //

        // ****************************************************************************************************** //

        // Output some additional information for debugging:
        Info << endl << endl;
        Info << "Current particle position: " << endl;
        Info << particleCloud.position(0) << endl;
        Info << "Current particle velocity: " << endl;
        Info << mag(particleCloud.velocity(0)) << endl;

        scalar xMax  = 0.01;
        scalar vMin = 1e-6;
        if (particleCloud.position(0)[0] > xMax )
        {
             Info<< "Particle has escaped" << endl;
             FatalError
                << "Particle has escaped" << nl
                << exit(FatalError);
        }

        if (mag(particleCloud.velocity(0)) < vMin )
        {
             Info<< "Particle has been captured" << endl;
             FatalError
                << "Particle has been captured" << nl
                << exit(FatalError);
        }

        if (particleCloud.numberOfParticles() == 0 )
        {
             Info<< "Particle has been lost" << endl;
             FatalError
                << "Particle has been lost" << nl
                << exit(FatalError);
        }

        // *******************************************************************************************************//
        // ADDED 03/02/2016
        // ****************************************************************************************** //
        // Volume check (Aycock)
        scalar voidFractionV = fvc::domainIntegrate(1-voidfractionNext).value();
        scalar totalNominalV = 0;
        for(int particleI = 0; particleI < particleCloud.numberOfParticles(); particleI++)
        {
            scalar radius=particleCloud.radius(particleI);
            scalar nominalV = 4./3.*(M_PI)*radius*radius*radius;
            totalNominalV += nominalV;
        }
        scalar volumeRatio = voidFractionV/totalNominalV;
        Info << "##################################" << endl;
        Info << "Total nominal volume: " << totalNominalV << endl;
        Info << "Total particle volume: " << voidFractionV << endl;
        Info << "Volume ratio (current / nominal): " << volumeRatio << endl;
        Info << "##################################" << endl << endl;

        // You may need to relax the volume ratio 'tolerance' below 0.99 here depending on the refinement of your mesh:
        if(volumeRatio < 0.9)
        { 
            Info << "WARNING: Particle volume has decreased more than 10\% below the nominal value" << endl;

            //FatalError << "Particle volume has decreased more than 1\% below the nominal value" << endl
            //           << "(possible locate model or void fraction error)" << exit(FatalError);
        }
        // ****************************************************************************************** //

        #if defined(version23)
        fvOptions.correct(U);
        #endif

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
