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


#include "engineSearchIB.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(engineSearchIB, 0);

addToRunTimeSelectionTable
(
    engineSearch,
    engineSearchIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
engineSearchIB::engineSearchIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    engineSearch(dict.subDict(typeName + "Props"),sm),
    propsDict_(dict.subDict(typeName + "Props")),
    zSplit_(readLabel(propsDict_.lookup("zSplit"))),
    xySplit_(readLabel(propsDict_.lookup("xySplit"))),
    checkPeriodicCells_(false)
{

    if(propsDict_.found("checkPeriodicCells")) checkPeriodicCells_=true;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

engineSearchIB::~engineSearchIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label engineSearchIB::findCell
(
    double** const& mask,
    double**& positions,
    double**& cellIDs,
    int size
) const
{

    int numprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const boundBox& globalBb = particleCloud_.mesh().bounds();


    vector position;
    for(int index = 0;index < size; ++index)
    {
        cellIDs[index][0]=-1;
        double radius=particleCloud_.radius(index);
        //if(mask[index][0] && radius > SMALL)
        if(radius > SMALL)
        {
            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];

            // find cell
            label oldID = cellIDs[index][0];

            // Test to see if the particle center is located in the current processor:
            cellIDs[index][0] = findSingleCell(position,oldID);
            //cellIDs[index][0] = particleCloud_.mesh().findCell(position);

            //mod by alice upon from here
            if(cellIDs[index][0]<0)
            {
                // Outline this procedure
                // Either there is a bug, or this method has some limitations
                // that cause some processors to "miss" the particle for a few timesteps
                // when a particle crosses processor boundaries (Aycock)
            	bool foundPos    = 0;
            	int  countPoints = 0;
                // Position vector (center of particle?)
            	vector pos=position;
            	label  altStartPos=-1;
            	label  numberOfPoints = zSplit_*xySplit_-2*(zSplit_-1);
            	label  thetaLevel=0;
            	//scalar theta, phi, thetaSize=180/zSplit_,phiSize=360/xySplit_, factor=M_PI/180.;
                // Fix by Aycock: We need the decimal points!!! Otherwise these values are evaluated
                // as integers.
            	scalar theta, phi, thetaSize=180./zSplit_,phiSize=360./xySplit_, factor=M_PI/180.;  // Aycock
             
                // "Satellite" points on the sphere surface are checked to see if the 
                // particle can be found in the current processor.
            	while(countPoints<numberOfPoints && !foundPos)
                {
                    pos=position;
                    // What's going on here?
                    if(countPoints==0)
                        pos[2]+=radius;
                    else if(countPoints==1)
                        pos[2]-=radius; 
                    // Create the satellite point:
                    else 
                    {
                        thetaLevel=(countPoints-2)/xySplit_;
                        //theta=factor*thetaSize*thetaLevel;
                        theta=factor*thetaSize*(thetaLevel+1); // Aycock
                        phi=factor*phiSize*(countPoints-2-thetaLevel*xySplit_);
                        pos[0]+=radius*sin(theta)*cos(phi);
                        pos[1]+=radius*sin(theta)*sin(phi);
                        pos[2]+=radius*cos(theta);
                    }
                   
                    // Check to see if the satellite point is contained in the current processor:
                    altStartPos=findSingleCell(pos,oldID); 
                    //particleCloud_.mesh().findCell(pos);//

                    //check for periodic domains
                    if(checkPeriodicCells_)
                    {
                        for(int iDir=0;iDir<3;iDir++)
                        {
                            if( pos[iDir] > globalBb.max()[iDir] )
                            {
                                pos[iDir]-=globalBb.max()[iDir]-globalBb.min()[iDir];
                            }
                            else if( pos[iDir] < globalBb.min()[iDir] )
                            {
                                pos[iDir]+=globalBb.max()[iDir]-globalBb.min()[iDir];
                            }
                        }
                  	altStartPos=findSingleCell(pos,oldID);
                        //particleCloud_.mesh().findCell(pos);//
                    }
                    // Was the satellite point found on the current processor?
                    if(altStartPos>=0) foundPos=1;
         
                    // Advance the counter, and check the next satellite point
            	    countPoints++;
            	}
                if(foundPos) cellIDs[index][0]=altStartPos;
            }
        }
    }
    return 1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
