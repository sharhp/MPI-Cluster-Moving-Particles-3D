#!/bin/bash
#PBS -N assn3
#PBS -q courses
#PBS -l nodes=4:ppn=4
#PBS -l walltime=02:00:00
#merge output and error into a single job_name.number_of_job_in_queue.
#PBS -j oe
#export fabric infiniband related variables
export I_MPI_FABRICS=shm:tmi
export I_MPI_DEVICE=rdma:OpenIB-cma
#change directory to where the job has been submitted from
cd $PBS_O_WORKDIR
#source paths
source /opt/software/intel17_update4/initpaths intel64
#run the job on required number of cores
sort $PBS_NODEFILE > hostfile
mpicc -o src.x src.c
iters=5


# Remove previous results
rm -rf ./hpc/data1/*
rm -rf ./hpc/data2/*

numP=(1 2 4 8 16)
run=(17 16)
numH=(1 1 2 4 8)
path=("./data/data1" "./data/data2")

for r in {0..1}
do
	echo -e "\nWORKING ON " "${path[$r]}" "\n"
	for np in {0..4}
	do	
		for i in `seq 1 $iters`
		do
			if [ ! -z $DEBUG ]
			then
				printf "Running for process size: ${numP[$np]} ($i of $iters iteration(s))...\n"
			fi

			mpiexec -f hostfile -np "${numP[$np]}" ./src.x "${run[$r]}" "${path[$r]}" "hpc"

			if [ ! -z $DEBUG ]
			then
				printf "Done!\n\n"
			fi
		done
	done
done

if [ $# -eq 0 ]
then
    unset DEBUG
fi