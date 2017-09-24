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

#include "ShirgaonkarIBperiodicGradP.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

//#include "pressureGradientExplicitSource.H" // ???

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShirgaonkarIBperiodicGradP, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ShirgaonkarIBperiodicGradP,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ShirgaonkarIBperiodicGradP::ShirgaonkarIBperiodicGradP
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
    voidFraction_(sm.mesh().lookupObject<volScalarField> ("voidfractionNext"))
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "shirgaonkarIB.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ShirgaonkarIBperiodicGradP::~ShirgaonkarIBperiodicGradP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ShirgaonkarIBperiodicGradP::setForce() const
{

    label cellI;
    vector drag;

    #include "setupProbeModel.H"

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        drag=vector::zero;

        scalar gradP0_ = 0;

        const dictionary& propsDict = rho_.mesh().lookupObject<IOdictionary>("myGradP");
        gradP0_ = readScalar(propsDict.lookup("gradient"));
        Info<< "Gradient read:\t" << gradP0_ << endl;


        for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
        {
            //Info << "subCell=" << subCell << endl;
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI > -1) // particle Found
            {
                // We need to know the flow direction! here, we are assuming it is (1 0 0)
                // Using Gauss law (divergence theorem), we can convert the surface integral into a volume integral, which is more convenient:  
                drag += vector(gradP0_ * rho_.mesh().V()[cellI], 0, 0);
            }

        } 

        // *********************************** //
        // FILE OUTPUT:
        // *********************************** //
        vector forceTemp = drag;
        reduce(forceTemp, sumOp<vector>());
        if(Pstream::master())
        {
            std::ofstream forceFileOutput;
            forceFileOutput.open ("f_ShirgaonkarIBperiodicGradP.csv", std::ofstream::app);
            forceFileOutput << rho_.mesh().time().value() << "\t" << forceTemp[0] << "\t" << forceTemp[1] << "\t" << forceTemp[2] << "\n";
            forceFileOutput.close();
        }
        // *********************************** //
        // *********************************** //

        // ************************************************************* // 
        // Calculate the force from the pressure gradient when using a 
        // momentum source - AYCOCK
        // ************************************************************* //
        //if(calcPeriodicPressureForce_) {
        //    const dictionary& propsDict = h.mesh().lookupObject<IOdictionary>("myGradP");
        //    scalar gradP0_ = readScalar(propsDict.lookup("gradient"));

        //    Info<< "Gradient read:\t" << gradP0_ << endl;
        //    
        //    if(Pstream::master())
        //    {
        //        scalar particleRadius = particleCloud_.radius(index);
        //        scalar particleVolume = 4./3. * M_PI * particleRadius * particleRadius * particleRadius;
 
        //        // We need to know the flow direction! here, we are assuming it is (1 0 0)
        //        // Using Gauss law (divergence theorem), we can convert the surface integral into a volume integral, which is more convenient:
        //        scalar pressureGradientForce = particleVolume * gradP0_ * rho_[0];
        //        Info<< "gradP0_:\t" << gradP0_ << "\t particleVolume:\t" << particleVolume << "\t density:\t" << rho_[0] << endl;
        //        Info<< "Force from pressure gradient ('pressureGradientExplicitSource' momentum source):\t" << pressureGradientForce << endl;
        //        drag += vector(pressureGradientForce, 0,  0);  // Assume uniform rho (incompressible)
        //    }
        //}
        // ************************************************************* //
        // END
        // ************************************************************* //
        
        // set force on particle
        if(twoDimensional_) drag /= depth_;

        //Set value fields and write the probe
        if(probeIt_)
        {
            #include "setupProbeModelfields.H"
            vValues.append(drag);           //first entry must the be the force
            particleCloud_.probeM().writeProbe(index, sValues, vValues);
        }

        if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
        else  for(int j=0;j<3;j++) impForces()[index][j] += drag[j];
        for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
        //Info<< "Drag (master thread): " << drag << endl;

        if(verbose_) Info << "impForces = " << impForces()[index][0]<<","<<impForces()[index][1]<<","<<impForces()[index][2] << endl;
        
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
