/////////////////////////////////////////////////////////////////////////////
/// @file SimulationClass.cpp
/// @brief Interface to the whole Simulation procedure:
/// from wave settings over lattice construction, time evolution and outputs
/// (also all relevant CVODE steps are performed here)
/////////////////////////////////////////////////////////////////////////////

#include "SimulationClass.h"

#include <math.h>

/// Along with the simulation object, create the cartesian communicator and
/// SUNContext object
Simulation::Simulation(const int nx, const int ny, const int nz,
        const int StencilOrder, const bool periodicity) :
    lattice(StencilOrder){
  statusFlags = 0;
  t = 0;
  // Initialize the cartesian communicator
  lattice.initializeCommunicator(nx, ny, nz, periodicity);

  // Create the SUNContext object associated with the thread of execution
  int retval = 0;
  retval = SUNContext_Create(&lattice.comm, &lattice.sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);
}

/// Free the CVode solver memory and Sundials context object with the finish of
/// the simulation
Simulation::~Simulation() {
  // Free solver memory
  if (statusFlags & CvodeObjectSetUp) {
    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    SUNContext_Free(&lattice.sunctx);
  }
}

/// Set the discrete dimensions, the number of points per dimension
void Simulation::setDiscreteDimensionsOfLattice(const sunindextype nx,
        const sunindextype ny, const sunindextype nz) {
  checkNoFlag(LatticePatchworkSetUp);
  lattice.setDiscreteDimensions(nx, ny, nz);
  statusFlags |= LatticeDiscreteSetUp;
}

/// Set the physical dimensions with lenghts in micro meters
void Simulation::setPhysicalDimensionsOfLattice(const sunrealtype lx,
        const sunrealtype ly, const sunrealtype lz) {
  checkNoFlag(LatticePatchworkSetUp);
  lattice.setPhysicalDimensions(lx, ly, lz);
  statusFlags |= LatticePhysicalSetUp;
}

/// Check that the lattice dimensions are set up and generate the patchwork
void Simulation::initializePatchwork(const int nx, const int ny,
        const int nz) {
  checkFlag(LatticeDiscreteSetUp);
  checkFlag(LatticePhysicalSetUp);

  // Generate the patchwork
  generatePatchwork(lattice, latticePatch, nx, ny, nz);
  latticePatch.initializeBuffers();

  statusFlags |= LatticePatchworkSetUp;
}

/// Configure CVODE
void Simulation::initializeCVODEobject(const sunrealtype reltol,
        const sunrealtype abstol) {
  checkFlag(SimulationStarted);

  // CVode settings return value
  int retval = 0;

  // Create CVODE object -- returns a pointer to the cvode memory structure
  // with Adams method (Adams-Moulton formula) solver chosen for non-stiff ODE
  cvode_mem = CVodeCreate(CV_ADAMS, lattice.sunctx);

  // Specify user data and attach it to the main cvode memory block
  retval = CVodeSetUserData(
      cvode_mem,
      &latticePatch); // patch contains the user data as used in CVRhsFn
  if (check_retval(&retval, "CVodeSetUserData", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);

  // Initialize CVODE solver
  retval = CVodeInit(cvode_mem, TimeEvolution::f, 0,
                     latticePatch.u); // allocate memory, CVRhsFn f, t_i=0, u
                                      // contains the initial values
  if (check_retval(&retval, "CVodeInit", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);

  // Create fixed point nonlinear solver object (suitable for non-stiff ODE) and
  // attach it to CVode
  NLS = SUNNonlinSol_FixedPoint(latticePatch.u, 0, lattice.sunctx);
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if (check_retval(&retval, "CVodeSetNonlinearSolver", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);

  // Anderson damping factor
  retval = SUNNonlinSolSetDamping_FixedPoint(NLS,1);
  if (check_retval(&retval, "SUNNonlinSolSetDamping_FixedPoint", 1,
              lattice.my_prc)) MPI_Abort(lattice.comm, 1);

  // Specify integration tolerances -- a scalar relative tolerance and scalar
  // absolute tolerance
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);

  // Specify the maximum number of steps to be taken by the solver in its
  // attempt to reach the next tout
  retval = CVodeSetMaxNumSteps(cvode_mem, 10000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);

  // maximum number of warnings for too small h
  retval = CVodeSetMaxHnilWarns(cvode_mem,3);
  if (check_retval(&retval, "CVodeSetMaxHnilWarns", 1, lattice.my_prc))
    MPI_Abort(lattice.comm, 1);

  statusFlags |= CvodeObjectSetUp;
}

/// Check if the lattice patchwork is set up and set the initial conditions
void Simulation::start() {
  checkFlag(LatticeDiscreteSetUp);
  checkFlag(LatticePhysicalSetUp);
  checkFlag(LatticePatchworkSetUp);
  checkNoFlag(SimulationStarted);
  checkNoFlag(CvodeObjectSetUp);
  setInitialConditions();
  statusFlags |= SimulationStarted;
}

/// Set initial conditions: Fill the lattice points with the initial field
/// values
void Simulation::setInitialConditions() {
  const sunrealtype dx = latticePatch.getDelta(1);
  const sunrealtype dy = latticePatch.getDelta(2);
  const sunrealtype dz = latticePatch.getDelta(3);
  const int nx = latticePatch.discreteSize(1);
  const int ny = latticePatch.discreteSize(2);
  const sunrealtype x0 = latticePatch.origin(1);
  const sunrealtype y0 = latticePatch.origin(2);
  const sunrealtype z0 = latticePatch.origin(3);
  int px = 0, py = 0, pz = 0;
  // space coordinates
  for (int i = 0; i < latticePatch.discreteSize() * 6; i += 6) {
    px = (i / 6) % nx;
    py = ((i / 6) / nx) % ny;
    pz = ((i / 6) / nx) / ny;
    // Call the `eval` function to fill the lattice points with the field data
    icsettings.eval(static_cast<sunrealtype>(px) * dx + x0,
            static_cast<sunrealtype>(py) * dy + y0,
            static_cast<sunrealtype>(pz) * dz + z0, &latticePatch.uData[i]);
  }
  return;
}

/// Use parameters to add periodic IC layers
void Simulation::addInitialConditions(const int xm, const int ym,
        const int zm /* zm=0 always */ ) {
  const sunrealtype dx = latticePatch.getDelta(1);
  const sunrealtype dy = latticePatch.getDelta(2);
  const sunrealtype dz = latticePatch.getDelta(3);
  const int nx = latticePatch.discreteSize(1);
  const int ny = latticePatch.discreteSize(2);
  // Correct for demanded displacement, rest as for setInitialConditions
  const sunrealtype x0 = latticePatch.origin(1) + xm*lattice.get_tot_lx();
  const sunrealtype y0 = latticePatch.origin(2) + ym*lattice.get_tot_ly();
  const sunrealtype z0 = latticePatch.origin(3) + zm*lattice.get_tot_lz();
  int px = 0, py = 0, pz = 0;
  for (int i = 0; i < latticePatch.discreteSize() * 6; i += 6) {
    px = (i / 6) % nx;
    py = ((i / 6) / nx) % ny;
    pz = ((i / 6) / nx) / ny;
    icsettings.add(static_cast<sunrealtype>(px) * dx + x0,
            static_cast<sunrealtype>(py) * dy + y0,
            static_cast<sunrealtype>(pz) * dz + z0, &latticePatch.uData[i]);
  }
  return;
}

/// Add initial conditions in one dimension
void Simulation::addPeriodicICLayerInX() {
  addInitialConditions(-1, 0, 0);
  addInitialConditions(1, 0, 0);
  return;
}

/// Add initial conditions in two dimensions
void Simulation::addPeriodicICLayerInXY() {
  addInitialConditions(-1, -1, 0);
  addInitialConditions(-1, 0, 0);
  addInitialConditions(-1, 1, 0);
  addInitialConditions(0, 1, 0);
  addInitialConditions(0, -1, 0);
  addInitialConditions(1, -1, 0);
  addInitialConditions(1, 0, 0);
  addInitialConditions(1, 1, 0);
  return;
}

/// Advance the solution in time -> integrate the ODE over an interval t
void Simulation::advanceToTime(const sunrealtype &tEnd) {
  checkFlag(SimulationStarted);
  int flag = 0;
  flag = CVode(cvode_mem, tEnd, latticePatch.u, &t,
               CV_NORMAL); // CV_NORMAL: internal steps to reach tEnd, then
                           // interpolate to return latticePatch.u, return time
                           // reached by the solver as t
  if (flag != CV_SUCCESS)
    printf("CVode failed, flag=%d.\n", flag);
}

/// Write specified simulations steps to disk
void Simulation::outAllFieldData(const int & state) {
  checkFlag(SimulationStarted);
  outputManager.outUState(state, lattice, latticePatch);
}

/// Check the presence configuration flags
void Simulation::checkFlag(unsigned int flag) const {
  if (!(statusFlags & flag)) {
    string errorMessage;
    switch (flag) {
    case LatticeDiscreteSetUp:
      errorMessage = "The discrete size of the Simulation has not been set up";
      break;
    case LatticePhysicalSetUp:
      errorMessage = "The physical size of the Simulation has not been set up";
      break;
    case LatticePatchworkSetUp:
      errorMessage = "The patchwork for the Simulation has not been set up";
      break;
    case CvodeObjectSetUp:
      errorMessage = "The CVODE object has not been initialized";
      break;
    case SimulationStarted:
      errorMessage = "The Simulation has not been started";
      break;
    default:
      errorMessage = "Uppss, you've made a non-standard error, sadly I can't "
                     "help you there";
      break;
    }
    errorKill(errorMessage);
  }
  return;
}

/// Check the absence of configuration flags
void Simulation::checkNoFlag(unsigned int flag) const {
  if ((statusFlags & flag)) {
    string errorMessage;
    switch (flag) {
    case LatticeDiscreteSetUp:
      errorMessage =
          "The discrete size of the Simulation has already been set up";
      break;
    case LatticePhysicalSetUp:
      errorMessage =
          "The physical size of the Simulation has already been set up";
      break;
    case LatticePatchworkSetUp:
      errorMessage = "The patchwork for the Simulation has already been set up";
      break;
    case CvodeObjectSetUp:
      errorMessage = "The CVODE object has already been initialized";
      break;
    case SimulationStarted:
      errorMessage = "The simulation has already started, some changes are no "
                     "longer possible";
      break;
    default:
      errorMessage = "Uppss, you've made a non-standard error, sadly I can't "
                     "help you there";
      break;
    }
    errorKill(errorMessage);
  }
  return;
}
