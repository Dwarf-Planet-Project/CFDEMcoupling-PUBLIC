/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])

define(calc, [esyscmd(perl -e 'printf ($1)')])

define(VCOUNT, 0)

define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

convertToMeters 1;

define(xMin,  -5)
define(xMax,  5)

define(yMin,  -2.5)
define(yMax,  2.5)

define(zMin,  -2.5)
define(zMax,  2.5)

define(res, 5)
define(xRes, calc( int( (xMax - xMin) * res) ) )
define(yRes, calc( int( (yMax - yMin) * res) ) )
define(zRes, calc( int( (zMax - zMin) * res) ) )

vertices
(
    (xMin yMin zMin)  
    (xMax yMin zMin)
    (xMax yMax zMin)
    (xMin yMax zMin)
    (xMin yMin zMax)
    (xMax yMin zMax)
    (xMax yMax zMax)
    (xMin yMax zMax)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (xRes yRes zRes) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
);

// ************************************************************************* //
