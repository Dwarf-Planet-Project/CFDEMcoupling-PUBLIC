/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "fileName.H"
#include "cfdemCloudIB.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "mpi.h"
#include "IOmanip.H"
#include <lammps.h>         // these are LAMMPS include files 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudIB::cfdemCloudIB
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    angularVelocities_(NULL),
    pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
    pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))),
    haveEvolvedOnce_(false),
    skipLagrangeToEulerMapping_(false),
    checkPeriodicCells_(false),
    rotationCoupling_(true) // Keep on by default to avoid messing up previous cases
{

    if(this->couplingProperties().found("skipLagrangeToEulerMapping"))
    {
        Info << "Will skip lagrange-to-Euler mapping..." << endl;
        skipLagrangeToEulerMapping_=true;
    }

    // ***************************** //
    // Additions from AYCOCK:
    // ***************************** //
    if(this->couplingProperties().found("checkPeriodicCells"))
    {
        Info << "Considering periodic BCs... (! only works for periodic x currently !)" << endl;
        checkPeriodicCells_=true;
    }

    if(this->couplingProperties().found("rotationCouplingOff"))
    {
        Info << "Turning off rotation coupling (DEM to CFD; for two way, add ShirgaonkarIBTorques) - AYCOCK" << endl;
        rotationCoupling_=false;
    } else
    {
    Info << "Turning on rotation coupling (DEM to CFD; for two way, add ShirgaonkarIBTorques) - AYCOCK" << endl; 
    }

    // ***************************** //

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudIB::~cfdemCloudIB()
{
    delete angularVelocities_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::cfdemCloudIB::getDEMdata()
{
    cfdemCloud::getDEMdata();
    if(rotationCoupling_) 
    {
        Info << "PARTICLE ROTATION COMMUNICATED FROM DEM TO CFD - AYCOCK" << endl; 
        dataExchangeM().getData("omega","vector-atom",angularVelocities_); 
    } else 
    {
        Info << "=== cfdemCloudIB::getDEMdata() === particle rotation not considered in CFD" << endl;
    }
}

bool Foam::cfdemCloudIB::reAllocArrays() const
{
    if(cfdemCloud::reAllocArrays())
    {
        dataExchangeM().allocateArray(angularVelocities_,0,3);
    }
    return true;
}

bool Foam::cfdemCloudIB::evolve()
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    // ***************************************************** //
    // The function "dataExchangeM().couple()" calls LIGGGHTS
    // We want to call "dataExchangeM().couple()" with 'run 0' to initialize the solver first. 
    // Then, we want to call "dataExchangeM().couple()" at the *end* of the timestep for each timestep.
    // (Aycock)
    // ***************************************************** //
    if (1) //dataExchangeM().couple()) // Delay the DEM solver from running
    { 
        // Only call DEM solver first on the first timestep (LIGGGHTS needs to initialize the particle variables; 'run 0')
        if(!haveEvolvedOnce_) {
            Info << "Performeing data exchange" << endl;
            dataExchangeM().couple();
            getDEMdata();
            locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
            voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_);
        }
        // Aycock

        Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
        doCouple=true;

//        Info << "skipLagrangeToEulerMapping_: " << skipLagrangeToEulerMapping_ 
//             << " haveEvolvedOnce_: " << haveEvolvedOnce_ << endl;
        // if(!skipLagrangeToEulerMapping_ || !haveEvolvedOnce_)
        // {
        //     if(verbose_) Info << "- getDEMdata()" << endl;
        //     getDEMdata();
        //     Info << "nr particles = " << numberOfParticles() << endl;
        // 
        //     // search cellID of particles
        //     // THERE IS A BUG IN locateM!!! (engineSearchIB; fixed as of 03/02/2016; Aycock)
        //     if(verbose_) Info << "- findCell()" << endl;
        //     locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
        //     if(verbose_) Info << "findCell done." << endl;

        //     // set void fraction field
        //     if(verbose_) Info << "- setvoidFraction()" << endl;
        //     voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_);
        //     if(verbose_) Info << "setvoidFraction done." << endl;
        // }
        
        // set particles forces
        if(verbose_) Info << "- setForce(forces_)" << endl;

        for(int index = 0;index <  numberOfParticles_; ++index){
            for(int i=0;i<3;i++){
                impForces_[index][i] = 0;
                expForces_[index][i] = 0;
                DEMForces_[index][i] = 0;
                DEMTorques_[index][i] = 0; // Aycock
            }
        }
        // Calculate the particle forces
        for (int i=0;i<nrForceModels();i++) forceM(i).setForce();
        if(verbose_) Info << "setForce done." << endl;

        // write DEM data
        if(verbose_) Info << " -giveDEMdata()" << endl;
        // Send the data to the DEM solver
        giveDEMdata();
        
        // Perform DEM AFTER the forces and moments are calculated:
        //Pout << "Data exchange result: " << endl;
        dataExchangeM().couple();

        // Reallocate arrays (in case the number of particles changed):
        //reAllocArrays();

        getDEMdata();
        locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
        voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_);

        haveEvolvedOnce_=true;
    }
    Info << "evolve done." << endl;

    //if(verbose_)    #include "debugInfo.H";

    // do particle IO
    // What does this do?
    IOM().dumpDEMdata();

    return doCouple;
}

void Foam::cfdemCloudIB::calcVelocityCorrection
(
    volScalarField& p,
    volVectorField& U,
    volScalarField& phiIB,
    volScalarField& voidfraction
)
{
    label cellI=0;
    vector uParticle(0,0,0);
    vector rVec(0,0,0);
    vector velRot(0,0,0);
    vector angVel(0,0,0);

    // ********************************************************************* //
    // Impose the 6DOF results on the velocity field:
    // ********************************************************************* //
    for(int index=0; index< numberOfParticles(); index++)
    {
        //if(regionM().inRegion()[index][0])
        //{
            for(int subCell=0;subCell<voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                //Info << "subCell=" << subCell << endl;
                cellI = cellIDs()[index][subCell];

                if (cellI >= 0)
                {
                    // calc particle velocity
                    vector positionI=position(index);
  
                    if(checkPeriodicCells_)
                    {
                        // ***************************** // 
                        // ADD CODE FOR PERIODIC BCs:
                        // ***************************** // 
                        // Generalized for x-dimension:
                        scalar r = radius(index);
                        vector cellPosition = U.mesh().C()[cellI]; 
                        // Expand r to make sure that we capture the boundary cells
                        r *= 1.5; 

                        if(mag(cellPosition - positionI) > r) {

                            const boundBox& globalBb = U.mesh().bounds();
                            vector xMax = globalBb.max();
                            vector xMin = globalBb.min();
                            vector xRange = xMax - xMin;

                            for(int i=0; i<3; i++) {
                                if(positionI[i] > (xMax[i] - r)        && cellPosition[i] < (xMin[i] + r) ) positionI[i] -= xRange[i];
                                else if (positionI[i] < (xMin[i] + r)  && cellPosition[i] > (xMax[i] - r) ) positionI[i] += xRange[i];
                            }
                        }
                        // ***************************** // 
                        // END CHANGES
                        // ***************************** // 
                    }
                       
                    //for(int i=0;i<3;i++) rVec[i]=U.mesh().C()[cellI][i]-position(index)[i];
                    for(int i=0;i<3;i++) rVec[i]=U.mesh().C()[cellI][i]-positionI[i];
                    for(int i=0;i<3;i++) angVel[i]=angularVelocities()[index][i];
                    velRot=angVel^rVec; // Cross product
                    for(int i=0;i<3;i++) uParticle[i] = velocities()[index][i]+velRot[i];

                    // impose field velocity
                    U[cellI]=(1-voidfractions_[index][subCell])*uParticle+voidfractions_[index][subCell]*U[cellI];
                }
            }
        //}
    }

    // ********************************************************************* //
    // Make the velocity field divergence free 
    // ********************************************************************* //
    volScalarField divUCheck = fvc::div(U);
    Info<< endl;
    Info<< "BEFORE CONTINUITY CORRECTION:" << endl;
    Info<< "Maximum per-cell div(U):\t" << max(divUCheck).value() << endl;
    Info<< "Minimum per-cell div(U):\t" << min(divUCheck).value() << endl;

    // Note, 03/18/2016: Iterating seems to cause problems with periodic Segre-Silberberg case! - Aycock
    // phiIB does not make the velocity field divergence-free in one pass, so iterate:
    for (int i = 0; i < 1; i++)
    {
        // Equation for the correction field 'phiIB' that will be used to make the velocity field divergence free
        fvScalarMatrix phiIBEqn
        (
            fvm::laplacian(phiIB) == fvc::div(U) //+ fvc::ddt(voidfraction)
        );

        // set reference value in case it is needed (the pressure?)
        if(phiIB.needReference()) 
        {
            phiIBEqn.setReference(pRefCell_, pRefValue_);
        }
        
        phiIBEqn.solve();

        // Use 'phiIB' to make the velocity field divergence free:
        U = U - fvc::grad(phiIB);
        U.correctBoundaryConditions();

        // Correct the pressure as well
        p = p + phiIB/U.mesh().time().deltaT();  // do we have to  account for rho here?
        //p = p + fvc::ddt(phiIB);  // <- why doesn't this work here? (oscillations during first few timesteps, then ~identical to the above)
        p.correctBoundaryConditions();
    }

    divUCheck = fvc::div(U);
    Info<< "AFTER CONTINUITY CORRECTION:" << endl;
    Info<< "Maximum per-cell div(U) (should be zero!):\t" << max(divUCheck).value() << endl;
    Info<< "Minimum per-cell div(U) (should be zero!):\t" << min(divUCheck).value() << endl;
    Info<< endl;
    // ********************************************************************* //

    if (couplingProperties_.found("checkinterface"))
    {
       Info << "checking no-slip on interface..." << endl;
       //#include "checkInterfaceVelocity.H" //TODO: check carefully!
    }

}

vector Foam::cfdemCloudIB::angularVelocity(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = angularVelocities_[index][i]; 
    return vel;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
