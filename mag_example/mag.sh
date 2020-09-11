#/bin/bash

mpiexec -l -machinefile hosts.txt ../src/grafen_rocm -dat3D mag_f.dat -Hf 100 -Hfrom -1 -Hto 0 -Hn 1 -l0 63 -dens magModelCube
