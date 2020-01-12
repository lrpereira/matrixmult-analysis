#!/bin/bash

#PBS -N AA_CUDA
#PBS -l walltime=03:00:00
#PBS -l nodes=1:r662:k20:ppn=48

 
source /share/apps/intel/parallel_studio_xe_2019/compilers_and_libraries_2019/linux/bin/compilervars.sh intel64

module purge

module load papi/5.4.1
module load gcc/4.8.2
module load cuda/7.0.28

export CUDA=yes

make clean 
make

StringVal = "dotProductCUDA dotProductBlockCUDA useCUDA"

for i in 32 128 1024 2048; do
	for func in $StringVal; do
		echo "Size=" $i >> "resultadosCUDA.csv"
		echo "\nFunction= " $func  >> "resultadosCUDA.csv"
		./bin/simple $func $i time
	done
done
