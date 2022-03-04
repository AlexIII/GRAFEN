#/bin/bash

mpiexec.mpich -l -machinefile hosts.txt ../src/grafen_rocm
