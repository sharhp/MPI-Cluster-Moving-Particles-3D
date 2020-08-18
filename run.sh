#!/bin/bash

make clean all

# If no argument specified set no of iterations to 5 by default
if [ $# -eq 0 ]
then
    iters=5
else
    iters=$1
    DEBUG=1
fi

# Remove previous results
rm -rf ./cse/data1/*
rm -rf ./cse/data2/*

# Compile source
mpicc -o src.x src.c

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
				echo "Generating hostfile..."
			fi
			./gen_hostfile.sh -n "${numH[$np]}"
			if [ ! -z $DEBUG ]
			then
				printf "hostfile generated!\n\n"
			fi

			if [ ! -z $DEBUG ]
			then
				printf "Running for process size: ${numP[$np]} ($i of $iters iteration(s))...\n"
			fi

			mpiexec -f hostfile -ppn 2 -np "${numP[$np]}" ./src.x "${run[$r]}" "${path[$r]}" "cse"

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
