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

#include "ShirgaonkarIBviscous.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

//#include "pressureGradientExplicitSource.H" // ???

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShirgaonkarIBviscous, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ShirgaonkarIBviscous,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ShirgaonkarIBviscous::ShirgaonkarIBviscous
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

ShirgaonkarIBviscous::~ShirgaonkarIBviscous()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ShirgaonkarIBviscous::setForce() const
{

    label cellI;
    vector drag;

    #ifdef comp // For compressible flow?
        // get viscosity field
        const volScalarField& mufField = particleCloud_.turbulence().mu();
        //volVectorField h = (mufField*fvc::laplacian(U_)-fvc::grad(p_));
        volVectorField h = mufField*fvc::laplacian(U_);
    #else
        // get viscosity field
        const volScalarField& nufField = particleCloud_.turbulence().nu();
        //volVectorField h = rho_*(nufField*fvc::laplacian(U_)-fvc::grad(p_)); // Pressure and viscous forces
        //volVectorField h = rho_*(-fvc::grad(p_)); // Pressure force only 
        volVectorField h = rho_*(nufField*fvc::laplacian(U_)); // Viscous force only
    #endif

    #include "setupProbeModel.H"

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        drag=vector::zero;

        for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
        {
            cellI = particleCloud_.cellIDs()[index][subCell];

            if (cellI > -1) // particle Found
            {
                drag += h[cellI]*h.mesh().V()[cellI];
                //drag += h[cellI]*h.mesh().V()[cellI] * (1-voidFraction_[cellI]);
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
            forceFileOutput.open ("f_ShirgaonkarIBviscous.csv", std::ofstream::app);
            forceFileOutput << h.mesh().time().value() << "\t" << forceTemp[0] << "\t" << forceTemp[1] << "\t" << forceTemp[2] << "\n";
            forceFileOutput.close();
        }
        // *********************************** //
        // *********************************** //

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
