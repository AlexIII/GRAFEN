# High-performance self-demagnetization modeling with GPUs

## Description

This software is supplimentary material for the paper Byzov D., Martyshko P., Chernoskutov A. Computationally Effective Modeling of Self-Demagnetization and Magnetic Field for Bodies of Arbitrary Shape Using Polyhedron Discretization

## Dependencies (Linux + ROCm)

- Ubuntu 20.04 (or anything that can run ROCm)
- mpich (or any other MPI implementation)
- ROCm 3.7+, see the [installation guide](https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html)
- [Boost 1.71.0](https://www.boost.org/users/history/version_1_71_0.html)

Specify paths to MPI in `src/rocm.makefile` and run `cd src && make` to build.

## Hardware requirements

- GPU
  - Nvidia GPU with CUDA compute capability 3.5 or higher (any modern 2015+ Nvidia GPU will do)
  
  OR
  - AMD GPU [with ROCm support](https://github.com/RadeonOpenCompute/ROCm#Hardware-and-Software-Support)
- At least 2GB GPU memory to run the example

## How to run the example

`mpiexec.mpich -l -machinefile hosts.txt src/grafen_rocm [arguments]`

arguments are:

### Ellipsoid model (git tag ellipsoid_example)

ellipEq, ellipPol - Ellipsoid equatorial and polar radii in meters

nl, nB, nR - number of partitioning elements along l, B and R axis of spherical coordinate system (associated with the ellipsoid)

K - magnetic susceptibility of the ellipsoid

HprimeX, HprimeY, HprimeZ - external field in nanotesla

Example:
`-ellipEq 10 -ellipPol 20    -nl 40 -nB 40 -nR 20    -K 2    -HprimeX 14 -HprimeY 14 -HprimeZ 35`

### Well model (git tag well_example)

xlower, xupper, xn, ylower, yupper, yn, zlower, zupper, zn - base cuboid bounds (with dense discretiztion grid) parameters

cxwell, cywell - well center

rwell - well radius

hwell - well height

Kouter - magnetic susceptibility outside of well

Kinner - magnetic susceptibility inside of well

HprimeX, HprimeY, HprimeZ - external field in nanotesla

fieldH - output field plane height (along z axis, relative to z=0)

fieldXfrom, fieldXto, fieldXn, fieldYfrom, fieldYto, fieldYn - output field parameters

### MPI nodes

Cerate `hosts.txt` file. You need to put here hosts that will execute the program. First host is 'root' host - it does not do actual computations. All other hosts perform computations using one GPU per host. You can utilize several GPUs on a single host by putting the same host entry several times in the file.
For example, if you have 2 GPUs on host 192.168.5.1 and 4 GPUs on host 192.168.5.2. Your `host.txt` should be as follows:
```
192.168.5.1
192.168.5.1
192.168.5.1
192.168.5.2
192.168.5.2
192.168.5.2
192.168.5.2
```
Note, that the first host has one extra entry compared to the amount of its GPUs.

## License

This software is distributed under MIT License. © Alexander Chernoskutov, Denis Byzov
