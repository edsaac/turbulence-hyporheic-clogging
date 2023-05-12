FOLDER=$1
NPROCS=$2
SOLVER=$3

runCase(){
	cd $1
	decomposePar
	mpirun -np $2 $3 -parallel > /dev/null 2>&1 
	#mpirun -np $2 $3 -parallel > log
	reconstructPar
	foamToVTK
	rm -r processor*
	cd -
}

runCase $FOLDER $NPROCS $SOLVER
