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

#include "gravityIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// For csv file output
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gravityIB, 0);

addToRunTimeSelectionTable
(
    forceModel,
    gravityIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gravityIB::gravityIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    twoDimensional_(false),
    particleRho_(readScalar(propsDict_.lookup("particleRho"))),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")), //mod by alice
    voidfractions_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),//mod by alice
    gravityFieldName_(propsDict_.lookup("gravityFieldName")),
    #if defined(version21) || defined(version16ext)
        g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_))
    #elif  defined(version15)
        g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
    #endif
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "gravityIBF.logDat");
    particleCloud_.probeM().vectorFields_.append("gravityIBForce");  //first entry must the be the force
    particleCloud_.probeM().writeHeader();  

    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
    }

    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    treatDEM_=true;
    Info << "accounting for Archimedes only on DEM side!" << endl;
    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gravityIB::~gravityIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gravityIB::setForce() const
{
    vector force;
 
    #include "setupProbeModel.H"
 
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            force=vector::zero;
          
            // NUMERICAL:
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                label cellI = particleCloud_.cellIDs()[index][subCell];
                if (cellI > -1) // particle Found
                {
                    // Calculate the gravitational body force numerically: 
                    force +=  g_.value()*particleRho_*rho_.mesh().V()[cellI]*(1-voidfractions_[cellI]);
                }
            } 

            // ANALYTICAL:
            //if(Pstream::master())
            //{
            //    scalar particleRadius = particleCloud_.radius(index);
            //    scalar particleVolume = 4./3. * M_PI * particleRadius * particleRadius * particleRadius;
            //    force  = particleRho_ * g_.value() * particleVolume; // Assume uniform rho
            //}

            // *********************************** //
            // FILE OUTPUT:
            // *********************************** //
            vector forceTemp = force;
            reduce(forceTemp, sumOp<vector>());
            if(Pstream::master())
            {
                std::ofstream gravityIBfile;
                gravityIBfile.open ("f_gravityIB.csv", std::ofstream::app);
                gravityIBfile << rho_.mesh().time().value() << "\t" << forceTemp[0] << "\t" << forceTemp[1] << "\t" << forceTemp[2] << "\n";
                gravityIBfile.close(); 
            }
            // *********************************** //
            // *********************************** //
            
            //Set value fields and write the probe
            if(probeIt_)
            {
                #include "setupProbeModelfields.H"
                vValues.append(force);           //first entry must the be the force
                particleCloud_.probeM().writeProbe(index, sValues, vValues);
            }

            // set force on particle
            if(twoDimensional_) Warning<<"gravityIB model doesn't work for 2D right now!!\n"<< endl;
            if(!treatDEM_)
            {
                if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += force[j];
                else for(int j=0;j<3;j++) impForces()[index][j] += force[j];
            }
            for(int j=0;j<3;j++) DEMForces()[index][j] += force[j];
        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
