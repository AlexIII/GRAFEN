#/bin/bash

mpiexec -l -machinefile hosts.txt ../src/grafen_rocm -dat2D timan.dat -Hf 0.00001 -Hfrom -81 -Hto 0 -Hn 81 -l0 63 -dens resmodel_timan -toRel -DPR 70
