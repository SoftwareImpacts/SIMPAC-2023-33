/////////////////////////////////////////////////////////////////////////////
/// @file SimulationClass.h
/// @brief Class for the Simulation object calling all functionality:
/// from wave settings over lattice construction, time evolution and outputs
/// initialization of the CVode object
/////////////////////////////////////////////////////////////////////////////

// Include Guard
#ifndef SIMULATIONCLASS
#define SIMULATIONCLASS

/* access to the fixed point SUNNonlinearSolver */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

// project subfile headers
#include "ICSetters.h"
#include "LatticePatch.h"
#include "Outputters.h"
#include "TimeEvolutionFunctions.h"

using namespace std;

/// simulation checking flags
enum SimulationOptions {
  LatticeDiscreteSetUp = 0x01,
  LatticePhysicalSetUp = 0x02,
  LatticePatchworkSetUp = 0x04, // not used anymore
  CvodeObjectSetUp = 0x08,
  SimulationStarted = 0x10
  /*OPT_B = 0x02,
  OPT_C = 0x04,
  OPT_D = 0x08,
  OPT_E = 0x10,
  OPT_F = 0x20,*/
};

/** @brief Simulation class to instantiate the whole walkthrough of a Simulation
 */
class Simulation {
private:
  /// Lattice object
  Lattice lattice;
  /// LatticePatch object
  LatticePatch latticePatch;
  /// current time of the simulation
  sunrealtype t;
  /// char for checking simulation flags
  unsigned char statusFlags;

public:
  /// IC Setter object
  ICSetter icsettings;
  /// Output Manager object
  OutputManager outputManager;
  /// Pointer to CVode memory object – public to avoid cross library errors
  void *cvode_mem;
  /// constructor function for the creation of the cartesian communicator
  Simulation(const int nx, const int ny, const int nz, const int StencilOrder,
          const bool periodicity);
  /// destructor function freeing CVode memory and Sundials context
  ~Simulation();
  /// Reference to the cartesian communicator of the lattice -> for debugging
  MPI_Comm *get_cart_comm() { return &lattice.comm; };
  /// function to set discrete dimensions of the lattice
  void setDiscreteDimensionsOfLattice(const sunindextype _tot_nx,
          const sunindextype _tot_ny, const sunindextype _tot_nz);
  /// function to set physical dimensions of the lattice
  void setPhysicalDimensionsOfLattice(const sunrealtype lx, const sunrealtype ly,
                                      const sunrealtype lz);
  /// function to initialize the Patchwork
  void initializePatchwork(const int nx, const int ny, const int nz);
  /// function to initialize the CVODE object with all requirements
  void initializeCVODEobject(const sunrealtype reltol,
                             const sunrealtype abstol);
  /// function to start the simulation for time iteration
  void start();
  /// functions to set the initial field configuration onto the lattice
  void setInitialConditions();
  /// functions to add initial periodic field configurations
  void addInitialConditions(const int xm, const int ym, const int zm = 0);
  /// function to add a periodic IC Layer in one dimension
  void addPeriodicICLayerInX();
  /// function to add periodic IC Layers in two dimensions
  void addPeriodicICLayerInXY();
  /// function to advance solution in time with CVODE
  void advanceToTime(const sunrealtype &tEnd);
  /// function to generate Output of the whole field at a given time
  void outAllFieldData(const int & state);
  /// function to check that a flag has been set and if not print an error
  // message and cause an abort on all ranks
  void checkFlag(unsigned int flag) const;
  /// function to check that if flag has not been set and if print an error
  // message and cause an abort on all ranks
  void checkNoFlag(unsigned int flag) const;
};

// End of Includeguard
#endif
