echo "Starting meshing..."

rm -r 1 2 3 processor* log.*

m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict

echo "Running blockMesh..."

blockMesh > log.blockMesh

echo "Decomposing domain..."

#decomposePar > log.decomposePar

echo "Running snappyHexMesh..."

surfaceFeatureExtract

snappyHexMesh > log.sHM

echo "Reconstructing mesh..."

#reconstructParMesh -latestTime > log.reconstructParMesh

echo "Checking mesh quality..."

checkMesh -latestTime > log.checkMesh

echo "Cleaning up..."

rm -r processor*
 
echo "Done!"
