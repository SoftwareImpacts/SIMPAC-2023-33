# HEWES: Heisenberg-Euler Weak-Field Expansion Simulator

The Heisenberg-Euler Weak-Field Expansion Simulator is a solver for the all-optical QED vacuum.
It solves the equations of motion for electromagnetic waves in the Heisenberg-Euler effective QED theory in the weak-field expansion with up to six-photon processes.

There is a [paper](https://arxiv.org/abs/2109.08121) that introduces the algorithm and shows remarkable scientific results.  
Check that out before the code if you are interested in this project!

![Harmonic Generation 3D](examples/figures/3d_harmonics.png)


## Contents
- [Preparing the Makefile](#preparing-the-makefile)
- [Short User Manual](#short-user-manual)
    - [Hints for Settings](#hints-for-settings)
    - [Note on Resource Occupation](#note-on-resource-occupation)
    - [Note on Output Analysis](#note-on-output-analysis)
- [Authors](#authors)


## Preparing the Makefile
The following descriptions assume you are using a Unix-like system.  
The `make` utility is used for building and a recent compiler version supporting OpenMP is required.
Features up to the C++20 standard are used.

Additionally required software:

- An MPI implementation.  
While `Intel(R) MPI` has mostly been used for scientific work on high-performance computing systems, the provided Makefile here is by default created for use with the _open_ implementations [`OpenMPI`](https://www.open-mpi.org/) or [`MPICH`](https://www.mpich.org).
While some useful Intel(R) processor specific optimizations and compiler options are not available with the latter, they are easier to get and set up on a personal device simply via the corresponding package manager. 

- The [`SUNDIALS`](https://computing.llnl.gov/projects/sundials) package with the [`CVode`](https://computing.llnl.gov/projects/sundials/cvode) solver.  
Version 6 is required.
The code is presumably compliant with the upcoming version 7.  
For the installation of `SUNDIALS`, [`CMake`](https://cmake.org/) is required.
Follow the installation guide and do not forget to enable MPI and specify the directory of the `mpicxx` wrapper for use of the MPI-based `NVECTOR_PARALLEL` module.
Make sure to edit the `SUNDIALS` include and library directories in the provided Makefile.


A minimal [Makefile](src/Makefile) template is provided.
Further compiler options might be beneficial, depending on the used system and software; e.g., higher vectorization and register usage instructions.


## Short User Manual
You have full control over all high-level simulation settings via the [`main.cpp`](src/main.cpp) file.

- First, specify the path you want the output data to go via the variable `outputDirectory`.

- Second, decide if you want to simulate in 1D, 2D, or 3D and uncomment only that full section.  
You can then specify
    - the relative and absolute integration tolerances of the CVode solver.  
    Recommended values are between 1e-12 and 1e-18.
    - the order of accuracy of the numerical scheme via the stencil order.  
    You can choose an integer in the range 1-13.
    - the physical side lengths of the grid in meters.
    - the number of lattice points per dimension.
    - the slicing of the lattice into patches (only for 2D and 3D simulations, automatic in 1D) – this determines the number of patches and therefore the required distinct processing units for MPI.  
    The total number of processes is given by the product of patches in any dimension.  
    Note: In the 3D case you better insure that every patch is cubic in terms of lattice points. This is decisive for computational efficiency.
    - whether to have periodic or vanishing boundary values (currently has to be chosen periodic).
    - whether you want to simulate on top of the linear vacuum only 4-photon processes (1), 6-photon processes (2), both (3), or none (0) – the linear Maxwell case.
    - the total time of the simulation in units c=1, i.e., the distance propagated by the light waves in meters.
    - the number of time steps that will be solved stepwise by CVode.   
    In order to keep interpolation errors small do not choose this number too small.
    - the multiple of steps at which you want the data to be written to disk.  
    - the output format. It can be 'c' for comma separated values (csv), or 'b' for binary.
    For csv format the name of the files written to the output directory is of the form `{step_number}_{process_number}.csv`.
    For binary output all data per step is written into one file and the step number is the name of
    the file.
    - which electromagnetic waveform(s) you want to propagate.  
    You can choose between a plane wave (not much physical content, but useful for checks) and implementations of Gaussians in 1D, 2D, and 3D.
    Their parameters can be tuned.  
    A description of the wave implementations is given in [ref.pdf](docs/ref.pdf).
    Note that the 3D Gaussians, as they are implemented up to now, should be propagated in the xy-plane.
    More waveform implementations will follow in subsequent versions of the code.

A doxygen-generated complete code reference is provided with [ref.pdf](docs/ref.pdf).

- Third, in the `src` directory, build the executable `Simulation` via the `make` command.

- Forth, run the simulation.  
Make sure to use `src` as working directory as the code uses a relative path to log the configuration in `main.cpp`.  
You determine the number of processes via the MPI execution command.
Note that in 2D and 3D simulations this number has to coincide with the actual number of patches, as described above.  
Here, the simulation would be executed distributed over four processes:
```bash
mpirun -np 4 ./Simulation
```

- Monitor stdout and stderr.
The unique simulation identifier number (starting timestep = name of data directory), the process steps, and the used wall times per step are printed on stdout.
Errors are printed on stderr.  
**Note**: Convergence of the employed CVode solver can not be guaranteed and issues of this kind can hardly be predicted.
On top, they are even system dependent.
Piece of advice: Only pass decimal numbers for the grid settings and initial conditions.  
CVode warnings and errors are reported on stdout and stderr.  
A `config.txt`file containing the relevant part of `main.cpp` is written to the output directory in order to save the simulation settings of each particular run.

You can remove the object files and the executable via `make clean`.


### Note on Simulation Settings
You may want to start with two Gaussian pulses in 1D colliding head-on in a pump-probe setup.
For this event, specify a high-frequency probe pulse with a low amplitude and a low-frequency pump pulse with a high frequency.
Both frequencies should be chosen to be below a forth of the Nyquist frequency to avoid nonphysical dispersion effects.
The wavelengths should neither be chosen too large (bulky wave) on a fine patchwork of narrow patches.
Their communication might be problematic with too small halo layer depths.
You would observe a blurring over time.
The amplitudes need be below 1 – the critical field strength – for the weak-field expansion to be valid.  
You can then investigate the arising of higher harmonics in frequency space via a Fourier analysis.
The signals from the higher harmonics can be highlighted by subtracting the results of the same simulation in the linear Maxwell vacuum.
You will be left with the nonlinear effects.  
Choosing the probe pulse to be polarized with an angle to the polarization of the pump you may observe a fractional polarization flip of the probe due to their nonlinear interaction.  
Decide beforehand which steps you need to be written to disk for your analysis.

Example scenarios of colliding Gaussians are preconfigured for any dimension.


### Note on Resource Occupation
The computational load depends mostly on the grid size and resolution.
The order of accuracy of the numerical scheme and CVode are rather secondary except for simulations running on many processing units, as the communication load is dependent on the stencil order.  
Simulations in 1D are relatively cheap and can easily be run on a modern laptop within minutes.
The output size per step is less than a megabyte.  
Simulations in 2D with about one million grid points are still feasible for a personal machine but might take about an hour of time to finish.
The output size per step is in the range of some dozen megabytes.  
Sensible simulations in 3D require large memory resources and therefore need to be run on distributed systems.
Even hundreds of cores can be kept busy for many hours or days.
The output size quickly amounts to dozens of gigabytes for just a single state.


### Note on Output Analysis
The field data are either written in csv format to one file per MPI process, the ending of which (after an underscore) corresponds to the process number, as described above.
This is not an elegant solution, but a portable way that also works fast and is
straightforward to analyze.  
Or, the option recommended for many larger write operations, in binary format with a single file per
output step.
Raw bytes are written to the files as they are in memory.
This option is more performant and achieved with MPI IO.
However, there is no guarantee of portability; postprocessing/conversion is required.
The step number is the file name.  
A `SimResults` folder is created in the chosen output directory if it does not exist and a folder named after the starting timestep of the simulation (in the form `yy-mm-dd_hh-MM-ss`) is created where the output files are written into.
There are six columns in the csv files, corresponding to the six components of the electromagnetic field: $`E_x`$, $`E_y`$, $`E_z`$, $`B_x`$, $`B_y`$, $`B_z`$.
Each row corresponds to one lattice point.  
Postprocessing is required to read-in the files in order.
A Python [module](examples/get_field_data.py) taking care of this is provided.  
Likewise, a Python [module](examples/get_binary_field_data.py) is provided to read the binary
data of a selected field component into a numpy array – its portability, however, cannot be guaranteed.  
The process numbers first align along dimension 1 until the number of patches is that direction is reached, then continue on dimension two and finally fill dimension 3.
For example, for a 3D simulation on 4x4x4=64 cores, the field data is divided over the patches as follows:
<pre>
z=1                          z=2                         z=3            z=4
                                                         ...            ...
x                            x 
  ^                            ^
1 | 0  4  8 12               1 |16 20 24 28
2 | 1  5  9 13               2 |17 21 25 29
3 | 2  6 10 14               3 |18 22 26 30
4 | 3  7 11 15               4 |19 23 27 31
   –––––––––––––>               –––––––––––––>
    1  2  3  4   y               1  2  3  4   y
</pre>
The axes denote the physical dimensions that are each divided into 4 sectors in this example.
The numbers inside the 4x4 squares indicate the process number, which is the number of the patch and also the number at the end of the corresponding output csv file.
The ordering of the array within a patch follows the standard C convention and can be reshaped in 2D and 3D to the actual size of the path.


More information describing settings and analysis procedures used for actual scientific results are given in the open-access [paper](https://arxiv.org/abs/2109.08121).  
Some example Python analysis scripts can be found in the [examples](examples).
The [first steps](examples/first_steps) demonstrate how the simulated data is accurately read-in from disk to numpy arrays using the provided [get field data module](examples/get_field_data.py).
[Harmonic generation](examples/harmonic_generation) in various forms is sketched as one application showing nonlinear quantum vacuum effects.
Analyses of 3D simulations are more involved due to large volumes of data.
Visualization requires tools like Paraview; examples are shown
[here](examples/3d_paraview_visualizations).
There is however _no simulation data provided_ as it would make the repository size unnecessarily large.



## Authors
- Arnau Pons Domenech
- Hartmut Ruhl (hartmut.ruhl@physik.uni-muenchen.de)
- Andreas Lindner (and.lindner@physik.uni-muenchen.de)
- Baris Ölmez (b.oelmez@physik.uni-muenchen.de)

