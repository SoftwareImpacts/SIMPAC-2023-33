#!/bin/bash

# Script to run executables built with CMake in the "src/build" directory
export EXECUTABLE="./build/hewes"

# Exporting MPI and OpenMP environment variables
export MPI_NUM_PROCESSES=4
export OMP_NUM_THREADS=1

# Specify the parameters
parameters=(\
    `# ---- General simulation settings -----`\
    /Users/andi/Documents `# output directory`\
    1 `# spatial dimensions of simulation space`\
    1.0e-16 `# CVode relative tolerance`\
    1.0e-16 `# CVode absolute tolerance`\
    13 `# stencil order`\
    300e-6 `# physical length of the simulation box in x-direction`\
    0. `# physical length of the simulation box in y-direction`\
    0. `# physical length of the simulation box in z-direction`\
    6000 `# number of lattice points in x-direction`\
    0 `# number of lattice points in y-direction`\
    0 `# number of lattice points in z-direction`\
    0 `# patches in x-direction`\
    0 `# patches in y-direction`\
    0 `# patches in z-direction`\
    3 `# linear vacuum (0), four-photon (1), six-photon (2), or 4- and 6-photon (3) processes`\
    100.0e-6 `# simulation time in meters`\
    100 `# number of simulation steps performed by CVode`\
    50 `# output step multiples`\
    csv `# output style; binary or csv`\
    `# ---- Parameters for plane waves -----`\
    0 `# use of plane wave(s) yes (1) or no (0)`\
    2 `# number of plane waves; 0, 1, or 2`\
    1e5 `# x-component of 1/lambda normalized wavevector of plane wave 1`\
    0. `# y-component of 1/lambda normalized wavevector of plane wave 1`\
    0. `# z-component of 1/lambda normalized wavevector of plane wave 1`\
    0. `# amplitude of plane wave 1 in x-direction`\
    0. `# amplitude of plane wave 1 in y-direction`\
    0.1 `# amplitude of plane wave 1 in z-direction`\
    0. `# phase shift of plane wave 1 in x-direction`\
    0. `# phase shift of plane wave 1 in y-direction`\
    0. `# phase shift of plane wave 1 in z-direction`\
    -1e6 `# normalized wavevector of plane wave 2 in x-direction`\
    0. `# normalized wavevector of plane wave 2 in y-direction`\
    0. `# normalized wavevector of plane wave 2 in z-direction`\
    0. `# amplitude of plane wave 2 in x-direction`\
    0. `# amplitude of plane wave 2 in y-direction`\
    0.5 `# amplitude of plane wave 2 in z-direction`\
    0. `# phase shift of plane wave 2 in x-direction`\
    0. `# phase shift of plane wave 2 in y-direction`\
    0. `# phase shift of plane wave 2 in z-direction`\
    `# ---- Parameters for 1D Gaussian pulses -----`\
    1 `# use of 1D Gaussian pulse(s) yes (1) or no (0)`\
    2 `# number of 1D Gaussian pulses; 0, 1, or 2`\
    1.0e6 `# x-component of 1/lambda normalized wavevector of 1D Gaussian 1`\
    0. `# y-component of 1/lambda normalized wavevector of 1D Gaussian 1`\
    0. `# z-component of 1/lambda normalized wavevector of 1D Gaussian 1`\
    0. `# amplitude of 1D Gaussian 1 in x-direction`\
    0. `# amplitude of 1D Gaussian 1 in y-direction`\
    0.1 `# amplitude of 1D Gaussian 1 in z-direction`\
    100e-6 `# shift from origin in x-direction of 1D Gaussian 1`\
    0. `# shift from origin in x-direction of 1D Gaussian 1`\
    0. `# shift from origin in x-direction of 1D Gaussian 1`\
    5e-6 `# width of 1D Gaussian 1`\
    0. `# phase shift of 1D Gaussian 1 in x-direction`\
    0. `# phase shift of 1D Gaussian 1 in y-direction`\
    0. `# phase shift of 1D Gaussian 1 in z-direction`\
    -0.2e6 `# x-component of 1/lambda normalized wavevector of 1D Gaussian 2`\
    0. `# y-component of 1/lambda normalized wavevector of 1D Gaussian 2`\
    0. `# z-component of 1/lambda normalized wavevector of 1D Gaussian 2`\
    0. `# amplitude of 1D Gaussian 2 in x-direction`\
    0. `# amplitude of 1D Gaussian 2 in y-direction`\
    0.5 `# amplitude of 1D Gaussian 2 in z-direction`\
    200e-6 `# shift from origin in x-direction of 1D Gaussian 2`\
    0. `# shift from origin in x-direction of 1D Gaussian 2`\
    0. `# shift from origin in x-direction of 1D Gaussian 2`\
    15e-6 `# width of 1D Gaussian 2`\
    0. `# phase shift of 1D Gaussian 2 in x-direction`\
    0. `# phase shift of 1D Gaussian 2 in y-direction`\
    0. `# phase shift of 1D Gaussian 2 in z-direction`\
    `# ---- Parameters for 2D/3D Gaussian pulses (not usable in 1D) -----`\
    0 `# use of 2D/3D Gaussian pulse(s) yes (1) or no (0)`\
    0 `# number of 2D/3D Gaussian pulses; 0, 1, or 2`\
    0 `# center of 2D/3D Gaussian 1 in x-direction`\
    0 `# center of 2D/3D Gaussian 1 in y-direction`\
    0 `# center of 2D/3D Gaussian 1 in z-direction`\
    0 `# x-component of normalized direction from which 2D/3D Gaussian 1 approaches the center`\
    0 `# y-component of normalized direction from which 2D/3D Gaussian 1 approaches the center`\
    0 `# amplitude of 2D/3D Gaussian 1`\
    0 `# polarization rotation from TE-mode (z-axis) in mulitples of pi/4`\
    0 `# taille of 2D/3D Gaussian 1`\
    0 `# Rayleigh length of 2D/3D Gaussian 1`\
    0 `# shift from center in negative propagation direction of 2D/3D Gaussian 1`\
    0 `# beam length of 2D/3D Gaussian 1`\
    0 `# center of 2D/3D Gaussian 2 in x-direction`\
    0 `# center of 2D/3D Gaussian 2 in y-direction`\
    0 `# center of 2D/3D Gaussian 2 in z-direction`\
    0 `# x-component of normalized direction from which 2D/3D Gaussian 2 approaches the center`\
    0 `# y-component of normalized direction from which 2D/3D Gaussian 2 approaches the center`\
    0 `# amplitude of 2D/3D Gaussian 2`\
    0 `# polarization rotation from TE-mode (z-axis) in mulitples of pi/4`\
    0 `# taille of 2D/3D Gaussian 2`\
    0 `# Rayleigh length of 2D/3D Gaussian 2`\
    0 `# shift from center in negative propagation direction of 2D/3D Gaussian 2`\
    0 `# beam length of 2D/3D Gaussian 2`\
)

# Run it
mpirun -np $MPI_NUM_PROCESSES $EXECUTABLE "${parameters[@]}"

#mpirun --use-hwthread-cpus -np $MPI_NUM_PROCESSES $EXECUTABLE "${parameters[@]}"
