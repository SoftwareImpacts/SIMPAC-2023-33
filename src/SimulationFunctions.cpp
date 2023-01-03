/////////////////////////////////////////////////////////////////////
/// @file SimulationFunctions.cpp
/// @brief Implementation of the complete simulation functions for
/// 1D, 2D, and 3D, as called in the main function
/////////////////////////////////////////////////////////////////////

#include "SimulationFunctions.h"

/** Calculate and print the elapsed wall time */
inline void timer(double &t1, double &t2) {
  printf("Elapsed time:  %fs\n", (t2 - t1));
}

// Instantiate and preliminarily initialize the time evolver
// non-const statics to be defined in actual simulation process
int *TimeEvolution::c = nullptr;
void (*TimeEvolution::TimeEvolver)(LatticePatch *, N_Vector, N_Vector,
                                   int *) = nonlinear1DProp;

/** Conduct the complete 1D simulation process */
void Sim1D(const std::array<sunrealtype,2> CVodeTol, const int StencilOrder,
        const sunrealtype phys_dim, const sunindextype disc_dim,
        const bool periodic, int *interactions,
        const sunrealtype endTime, const int numberOfSteps,
        const std::string outputDirectory, const int outputStep,
        const char outputStyle,
        const std::vector<planewave> &planes,
        const std::vector<gaussian1D> &gaussians) {

  // MPI data
  int myPrc = 0, nPrc = 1;
#if defined(_MPI)
  double ts = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &nPrc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPrc);

  // Check feasibility of the patchwork decomposition
  if (myPrc == 0) {
      if (disc_dim % nPrc != 0) {
        errorKill("The number of lattice points must be "
                  "divisible by the number of processes.");
      }
  }
#elif !defined(_MPI) && defined(_OPENMP)
    double ts = omp_get_wtime();
#endif

  // Initialize the simulation, set up the cartesian communicator
  std::array<int, 3> patches = {nPrc, 1, 1};
  Simulation sim(patches[0], patches[1], patches[2], StencilOrder, periodic);

  // Configure the patchwork
  sim.setPhysicalDimensionsOfLattice(phys_dim,1,1);
  sim.setDiscreteDimensionsOfLattice(disc_dim,1,1);
  sim.initializePatchwork(patches[0], patches[1], patches[2]);

  // Add em-waves
  for (const auto &gauss : gaussians)
    sim.icsettings.addGauss1D(gauss.k, gauss.p, gauss.x0, gauss.phig,
                              gauss.phi);
  for (const auto &plane : planes)
    sim.icsettings.addPlaneWave1D(plane.k, plane.p, plane.phi);

  // Check that the patchwork is ready and set the initial conditions
  sim.start();
  sim.addPeriodicICLayerInX();

  // Initialize CVode with abs and rel tolerances
  sim.initializeCVODEobject(CVodeTol[0], CVodeTol[1]);

  // Configure the time evolution function
  TimeEvolution::c = interactions;
  TimeEvolution::TimeEvolver = nonlinear1DProp;

  // Configure the output
  sim.outputManager.generateOutputFolder(outputDirectory);
  if (!myPrc) {
    std::cout << "Simulation code: " << sim.outputManager.getSimCode()
        << std::endl;
  }
  sim.outputManager.set_outputStyle(outputStyle);

  //MPI_Barrier(MPI_COMM_WORLD);  // insure correct output
  //sim.outAllFieldData(0);  // output of initial state
  // Conduct the propagation in space and time
  for (int step = 1; step <= numberOfSteps; step++) {
    sim.advanceToTime(endTime / numberOfSteps * step);
    if (step % outputStep == 0) {
#if defined(_MPI)
      MPI_Barrier(MPI_COMM_WORLD);  // insure correct output
#endif
      sim.outAllFieldData(step);
    }
#if defined(_MPI)
    double tn = MPI_Wtime();
#elif !defined(_MPI) && defined(_OPENMP)
    double tn = omp_get_wtime();
#endif
    if (!myPrc) {
      std::cout << "\rStep " << step << "\t\t" << std::flush;
#if defined(_MPI) || defined(_OPENMP)
      timer(ts, tn);
#endif
    }
  }

  return;
}

/** Conduct the complete 2D simulation process */
void Sim2D(const std::array<sunrealtype,2> CVodeTol, int const StencilOrder,
        const std::array<sunrealtype,2> phys_dims,
        const std::array<sunindextype,2> disc_dims,
        const std::array<int,2> patches, const bool periodic, int *interactions,
        const sunrealtype endTime, const int numberOfSteps,
        const std::string outputDirectory, const int outputStep,
        const char outputStyle,
        const std::vector<planewave> &planes,
        const std::vector<gaussian2D> &gaussians) {

  // MPI data
  int myPrc = 0, nPrc = 1; // Get process rank and number of processes
#if defined(_MPI)
  double ts = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD,
                &myPrc); // Return process rank, number \in [1,nPrc]
  MPI_Comm_size(MPI_COMM_WORLD,
                &nPrc); // Return number of processes (communicator size)

  // Check feasibility of the patchwork decomposition
  if (myPrc == 0) {
    if (nPrc != patches[0] * patches[1]) {
      errorKill(
          "The number of MPI processes must match the number of patches.");
    }
  }
#elif !defined(_MPI) && defined(_OPENMP)
    double ts = omp_get_wtime();
#endif

  // Initialize the simulation, set up the cartesian communicator
  Simulation sim(patches[0], patches[1], 1, StencilOrder, periodic);

  /* // Check that lattice communicator is unique; same as used in patchwork
  generation char cart_comm_name[MPI_MAX_OBJECT_NAME]; int cart_namelen;
  MPI_Comm_get_name(*sim.get_cart_comm(), cart_comm_name, &cart_namelen);
  printf("sim.get_cart_comm gives %s \n", cart_comm_name);
  */

  // Configure the patchwork
  sim.setPhysicalDimensionsOfLattice(phys_dims[0],
                                     phys_dims[1],
                                     1); // spacing of the lattice
  sim.setDiscreteDimensionsOfLattice(
      disc_dims[0], disc_dims[1], 1); // Spacing equivalence to points
  sim.initializePatchwork(patches[0], patches[1], 1);

  // Add em-waves
  for (const auto &gauss : gaussians)
    sim.icsettings.addGauss2D(gauss.x0, gauss.axis, gauss.amp, gauss.phip,
                              gauss.w0, gauss.zr, gauss.ph0, gauss.phA);
  for (const auto &plane : planes)
    sim.icsettings.addPlaneWave2D(plane.k, plane.p, plane.phi);

  // Check that the patchwork is ready and set the initial conditions
  sim.start(); // Check if the lattice is set up, set initial field
               // configuration
  sim.addPeriodicICLayerInXY(); // insure periodicity in propagation directions

  // Initialize CVode with rel and abs tolerances
  sim.initializeCVODEobject(CVodeTol[0], CVodeTol[1]);

  // Configure the time evolution function
  TimeEvolution::c = interactions;
  TimeEvolution::TimeEvolver = nonlinear2DProp;

  // Configure the output
  sim.outputManager.generateOutputFolder(outputDirectory);
  if (!myPrc) {
    std::cout << "Simulation code: " << sim.outputManager.getSimCode()
        << std::endl;
  }
  sim.outputManager.set_outputStyle(outputStyle);

  //MPI_Barrier(MPI_COMM_WORLD);  // insure correct output
  //sim.outAllFieldData(0);  // output of initial state
  // Conduct the propagation in space and time
  for (int step = 1; step <= numberOfSteps; step++) {
    sim.advanceToTime(endTime / numberOfSteps * step);
    if (step % outputStep == 0) {
#if defined(_MPI)
      MPI_Barrier(MPI_COMM_WORLD);  // insure correct output
#endif
      sim.outAllFieldData(step);
    }
#if defined(_MPI)
    double tn = MPI_Wtime();
#elif !defined(_MPI) && defined(_OPENMP)
    double tn = omp_get_wtime();
#endif
    if (!myPrc) {
      std::cout << "\rStep " << step << "\t\t" << std::flush;
#if defined(_MPI) || defined(_OPENMP)
      timer(ts, tn);
#endif
    }
  }

  return;
}

/** Conduct the complete 3D simulation process */
void Sim3D(const std::array<sunrealtype,2> CVodeTol, const int StencilOrder,
        const std::array<sunrealtype,3> phys_dims,
        const std::array<sunindextype,3> disc_dims,
        const std::array<int,3> patches,
        const bool periodic, int *interactions, const sunrealtype endTime,
        const int numberOfSteps, const std::string outputDirectory,
        const int outputStep, const char outputStyle,
        const std::vector<planewave> &planes,
        const std::vector<gaussian3D> &gaussians) {

  // MPI data
  int myPrc = 0, nPrc = 1; // Get process rank and numer of process
#if defined(_MPI)
  double ts = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD,
                &myPrc); // rank of the process inside the world communicator
  MPI_Comm_size(MPI_COMM_WORLD,
                &nPrc); // Size of the communicator is the number of processes

  // Check feasibility of the patchwork decomposition
  if (myPrc == 0) {
    if (nPrc != patches[0] * patches[1] * patches[2]) {
      errorKill(
          "The number of MPI processes must match the number of patches.");
    }
    if ( ( disc_dims[0] / patches[0] != disc_dims[1] / patches[1] ) |
         ( disc_dims[0] / patches[0] != disc_dims[2] / patches[2] ) ) {
         std::clog
          << "\nWarning: Patches should be cubic in terms of the lattice "
             "points for the computational efficiency of larger simulations.\n";
    }
  }
#elif !defined(_MPI) && defined(_OPENMP)
    double ts = omp_get_wtime();
#endif

  // Initialize the simulation, set up the cartesian communicator
  Simulation sim(patches[0], patches[1], patches[2],
                 StencilOrder, periodic); // Simulation object with slicing

  // Create the SUNContext object associated with the thread of execution
  sim.setPhysicalDimensionsOfLattice(phys_dims[0], phys_dims[1],
                                     phys_dims[2]); // spacing of the box
  sim.setDiscreteDimensionsOfLattice(
      disc_dims[0], disc_dims[1],
      disc_dims[2]); // Spacing equivalence to points
  sim.initializePatchwork(patches[0], patches[1], patches[2]);
  //sim.initializeGhostCells();

  // Add em-waves
  for (const auto &plane : planes)
    sim.icsettings.addPlaneWave3D(plane.k, plane.p, plane.phi);
  for (const auto &gauss : gaussians)
    sim.icsettings.addGauss3D(gauss.x0, gauss.axis, gauss.amp, gauss.phip,
                              gauss.w0, gauss.zr, gauss.ph0, gauss.phA);

  // Check that the patchwork is ready and set the initial conditions
  sim.start();

  // Initialize CVode with abs and rel tolerances
  sim.initializeCVODEobject(CVodeTol[0], CVodeTol[1]);

  // Configure the time evolution function
  TimeEvolution::c = interactions;
  TimeEvolution::TimeEvolver = nonlinear3DProp;

  // Configure the output
  sim.outputManager.generateOutputFolder(outputDirectory);
  if (!myPrc) {
    std::cout << "Simulation code: " << sim.outputManager.getSimCode()
        << std::endl;
  }
  sim.outputManager.set_outputStyle(outputStyle);

  //MPI_Barrier(MPI_COMM_WORLD);  // insure correct output
  //sim.outAllFieldData(0);  // output of initial state
  // Conduct the propagation in space and time
  for (int step = 1; step <= numberOfSteps; step++) {
    sim.advanceToTime(endTime / numberOfSteps * step);
    if (step % outputStep == 0) {
#if defined(_MPI)
      MPI_Barrier(MPI_COMM_WORLD);  // insure correct output
#endif
      sim.outAllFieldData(step);
    }
#if defined(_MPI)
    double tn = MPI_Wtime();
#elif !defined(_MPI) && defined(_OPENMP)
    double tn = omp_get_wtime();
#endif
    if (!myPrc) {
      std::cout << "\rStep " << step << "\t\t" << std::flush;
#if defined(_MPI) || defined(_OPENMP)
      timer(ts, tn);
#endif
    }
  }
  return;
}
