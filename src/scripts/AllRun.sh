SOLVER=$1
NPROC=$2

## Edit dictionaries
foamDictionary system/decomposeParDict -entry numberOfSubdomains -set $NPROC
foamDictionary system/controlDict -entry application -set $SOLVER

## Run openFoam
decomposePar
mpirun -np $NPROC $SOLVER -parallel > log 
reconstructPar
rm -r processor*
foamToVTK -noZero
