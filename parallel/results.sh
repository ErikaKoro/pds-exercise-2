#!/bin/bash

make clean
make mpi_sort_build

echo '*************************************'
echo "Running for Randomly generated data"
echo '*************************************'

for j in 2 4 8 16; do
    echo "127.0.0.1 slots=$j" > hosts
    echo "Running for $j processes"
    for i in 2 4 8 16 32; do
        mpirun -hostfile hosts ./build/mpi_sort.out data/bigbigData.bin $i
    done

    echo
    echo
    echo
done

echo '*************************************'
echo "Running for mnist generated data"
echo '*************************************'

for j in 2 4 8 16; do
    echo "127.0.0.1 slots=$j" > hosts
    echo "Running for $j processes"
    for i in 2 4 8; do
        mpirun -hostfile hosts ./build/mpi_sort.out data/mnist.bin $i
    done

    echo
    echo
    echo
done