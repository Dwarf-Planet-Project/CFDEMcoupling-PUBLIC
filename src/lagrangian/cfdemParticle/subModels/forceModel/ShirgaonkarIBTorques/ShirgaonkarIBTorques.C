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

#include "error.H"
#include "ShirgaonkarIBTorques.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShirgaonkarIBTorques, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ShirgaonkarIBTorques,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ShirgaonkarIBTorques::ShirgaonkarIBTorques
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    twoDimensional_(false),
    depth_(1),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_)),
    checkPeriodicCells_(false) // Aycock
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "shirgaonkarIBTorques.logDat");
    particleCloud_.probeM().vectorFields_.append("hdtorqe"); 
    particleCloud_.probeM().writeHeader();

    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        depth_ = propsDict_.lookup("depth");
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
        Info << "depth of domain is assumed to be :" << depth_ << endl;
    }
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    particleCloud_.checkCG(false);
    if(propsDict_.found("checkPeriodicCells")) checkPeriodicCells_=true; // Aycock
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ShirgaonkarIBTorques::~ShirgaonkarIBTorques()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ShirgaonkarIBTorques::setForce() const
{
    label cellI;
    vector torque;

    #ifdef comp
        // get viscosity field
        const volScalarField& mufField = particleCloud_.turbulence().mu();
        //volVectorField h = (mufField*fvc::laplacian(U_)-fvc::grad(p_));
        volVectorField h = (mufField*fvc::laplacian(U_)); // Boundary force, viscous contribution only
    #else
        // get viscosity field
        const volScalarField& nufField = particleCloud_.turbulence().nu();
        //volVectorField h = rho_*(nufField*fvc::laplacian(U_)-fvc::grad(p_)); // Boundary force, including viscous and pressure contributions
        volVectorField h = rho_*(nufField*fvc::laplacian(U_)); // Boundary force, viscous contribution only
    #endif

    #include "setupProbeModel.H"

    // For periodic:
    const boundBox& globalBb = particleCloud_.mesh().bounds();
    // END //Aycock

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        Info<< "Initializing torque vector to zeros: " << endl;
        torque=vector::zero;

        for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
        {
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI > -1) // particle Found
            {
                // Calculate the vector from the center of the particle to the current cell:
                vector positionI = particleCloud_.position(index);
                //scalar centreDist = mag(cellCentrePosition-positionCenter);

                // Add check to see if particle position needs to be modified to account for periodic BCs:
                if(checkPeriodicCells_)
                {
                    // ***************************** // 
                    // ADD CODE FOR PERIODIC BCs:
                    // ***************************** // 
                    // Generalized for x-dimension:
                    scalar r = particleCloud_.radius(index);
                    vector cellPosition = h.mesh().C()[cellI];
                    // Expand r to make sure that we capture the boundary cells
                    r *= 1.5;

                    if(mag(cellPosition - positionI) > r) {

                        const boundBox& globalBb = h.mesh().bounds();
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
                // END

                vector rVec(0,0,0);
                for(int i=0;i<3;i++) rVec[i]=h.mesh().C()[cellI][i]-positionI[i]; // calculate the vector from the particle center to the cell center (from cfdemCloudIB)
                vector boundaryForce(0,0,0);
                boundaryForce = h[cellI]*h.mesh().V()[cellI]; // This includes pressure and viscous forces; used for total drag in previous force model by Hager.
                vector localTorque = rVec^boundaryForce; // r x F
                torque += localTorque; 
            }

        }

        // *********************************** //
        // FILE OUTPUT:
        // *********************************** //
        vector forceTemp = torque;
        reduce(forceTemp, sumOp<vector>());
        if(Pstream::master())
        {
            std::ofstream forceFileOutput;
            forceFileOutput.open ("f_ShirgaonkarIBTorques.csv", std::ofstream::app);
            forceFileOutput << rho_.mesh().time().value() << "\t" << forceTemp[0] << "\t" << forceTemp[1] << "\t" << forceTemp[2] << "\n";
            forceFileOutput.close();
        }
        // *********************************** //
        // *********************************** //

        // Pass torque value to DEM data:
        for(int j=0;j<3;j++) DEMTorques()[index][j] += torque[j];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
