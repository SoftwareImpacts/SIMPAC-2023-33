//////////////////////////////////////////////////////////////////////////
/// @file LatticePatch.h
/// @brief Declaration of the lattice and lattice patches
//////////////////////////////////////////////////////////////////////////

// Include Guard
#ifndef LATTICEPATCH
#define LATTICEPATCH

// IO
#include <iomanip>
#include <iostream>
#include <sstream>

// string, container, algorithm
#include <string>
//#include <string_view>
#include <array>
#include <vector>
#include <algorithm>

// MPI & OpenMP
#include <mpi.h>
#include <omp.h>

// Sundials
#include <cvode/cvode.h>              /* prototypes for CVODE fcts. */
#include <nvector/nvector_parallel.h> /* definition of N_Vector and macros */
#include <sundials/sundials_types.h>  /* definition of type sunrealtype */

// stencils
#include "DerivationStencils.h"

using namespace std;

// lattice construction checking flags
enum LatticeOptions {
    FLatticeDimensionSet = 0x01,  // 1
    /*OPT_B = 0x02,  // 2
    OPT_C = 0x04,  // 4
    OPT_D = 0x08,  // 8
    OPT_E = 0x10,  // 16
    OPT_F = 0x20,*/  // 32
};

/** @brief Lattice class for the construction of the enveloping discrete
 * simulation space */
class Lattice {
private:
  /// physical size of the lattice in x-direction
  sunrealtype tot_lx;
  /// physical size of the lattice in y-direction
  sunrealtype tot_ly;
  /// physical size of the lattice in z-direction
  sunrealtype tot_lz;
  /// number of points in x-direction
  sunindextype tot_nx;
  /// number of points in y-direction
  sunindextype tot_ny;
  /// number of points in z-direction
  sunindextype tot_nz;
  /// total number of lattice points
  sunindextype tot_noP;
  /// dimension of each data point -> set once and for all
  static constexpr int dataPointDimension = 6;
  /// number of lattice points times data dimension of each point
  sunindextype tot_noDP;
  /// physical distance between lattice points in x-direction
  sunrealtype dx;
  /// physical distance between lattice points in y-direction
  sunrealtype dy;
  /// physical distance between lattice points in z-direction
  sunrealtype dz;
  /// stencil order
  const int stencilOrder;
  /// required width of ghost layers (depends on the stencil order)
  const int ghostLayerWidth;
  /// char for checking if lattice flags are set
  unsigned char statusFlags;

public:
  /// number of MPI processes
  int n_prc;
  /// number of MPI process
  int my_prc;
  /// personal communicator of the lattice
  MPI_Comm comm;
  /// function to create and deploy the cartesian communicator
  void initializeCommunicator(const int nx, const int ny,
          const int nz, const bool per);
  /// default construction
  Lattice(const int StO);
  /// SUNContext object
  SUNContext sunctx;
  /// component function for resizing the discrete dimensions of the lattice
  void setDiscreteDimensions(const sunindextype _nx,
          const sunindextype _ny, const sunindextype _nz);
  /// component function for resizing the physical size of the lattice
  void setPhysicalDimensions(const sunrealtype _lx,
          const sunrealtype _ly, const sunrealtype _lz);
  ///@{
  /** getter function */
  [[nodiscard]] const sunrealtype &get_tot_lx() const { return tot_lx; }
  [[nodiscard]] const sunrealtype &get_tot_ly() const { return tot_ly; }
  [[nodiscard]] const sunrealtype &get_tot_lz() const { return tot_lz; }
  [[nodiscard]] const sunindextype &get_tot_nx() const { return tot_nx; }
  [[nodiscard]] const sunindextype &get_tot_ny() const { return tot_ny; }
  [[nodiscard]] const sunindextype &get_tot_nz() const { return tot_nz; }
  [[nodiscard]] const sunindextype &get_tot_noP() const { return tot_noP; }
  [[nodiscard]] const sunindextype &get_tot_noDP() const { return tot_noDP; }
  [[nodiscard]] const sunrealtype &get_dx() const { return dx; }
  [[nodiscard]] const sunrealtype &get_dy() const { return dy; }
  [[nodiscard]] const sunrealtype &get_dz() const { return dz; }
  [[nodiscard]] constexpr int get_dataPointDimension() const {
    return dataPointDimension;
  }
  [[nodiscard]] const int &get_stencilOrder() const { return stencilOrder; }
  [[nodiscard]] const int &get_ghostLayerWidth() const {
    return ghostLayerWidth;
  }
  ///@}
};

/// lattice patch construction checking flags
enum LatticePatchOptions {
  FLatticePatchSetUp = 0x01,
  TranslocationLookupSetUp = 0x02,
  GhostLayersInitialized = 0x04,
  BuffersInitialized = 0x08
  /*OPT_D = 0x08,
  OPT_E = 0x10,
  OPT_F = 0x20,*/
};

/** @brief LatticePatch class for the construction of the patches in the
 * enveloping lattice */
class LatticePatch {
private:
  /// origin of the patch in physical space; x-coordinate
  sunrealtype x0;
  /// origin of the patch in physical space; y-coordinate
  sunrealtype y0;
  /// origin of the patch in physical space; z-coordinate
  sunrealtype z0;
  /// inner position of lattice-patch in the lattice patchwork; x-points
  sunindextype LIx;
  /// inner position of lattice-patch in the lattice patchwork; y-points
  sunindextype LIy;
  /// inner position of lattice-patch in the lattice patchwork; z-points
  sunindextype LIz;
  /// physical size of the lattice-patch in the x-dimension
  sunrealtype lx;
  /// physical size of the lattice-patch in the y-dimension
  sunrealtype ly;
  /// physical size of the lattice-patch in the z-dimension
  sunrealtype lz;
  /// number of points in the lattice patch in the x-dimension
  sunindextype nx;
  /// number of points in the lattice patch in the y-dimension
  sunindextype ny;
  /// number of points in the lattice patch in the z-dimension
  sunindextype nz;
  /// physical distance between lattice points in x-direction
  sunrealtype dx;
  /// physical distance between lattice points in y-direction
  sunrealtype dy;
  /// physical distance between lattice points in z-direction
  sunrealtype dz;
  /// pointer to the enveloping lattice
  const Lattice *envelopeLattice;
  ///@{
  /** translocation lookup table */
  vector<int> uTox, uToy, uToz, xTou, yTou, zTou;
  ///@}
  /// aid (auxilliarly) vector including ghost cells to compute the derivatives
  vector<sunrealtype> uAux;
  ///@{
  /** buffer to save spatial derivative values */
  vector<sunrealtype> buffX, buffY, buffZ;
  ///@}
  ///@{
  /** buffer for passing ghost cell data */
  vector<sunrealtype> ghostCellLeft, ghostCellRight, ghostCellLeftToSend,
      ghostCellRightToSend, ghostCellsToSend, ghostCells;
  ///@}
  ///@{
  /** ghost cell translocation lookup table */
  vector<int> lgcTox, rgcTox, lgcToy, rgcToy, lgcToz, rgcToz;
  ///@}
  /** char for checking flags */
  unsigned char statusFlags;
  ///@{
  /** rotate and translocate an input array according to a lookup into an output
   * array */
  inline void rotateToX(sunrealtype *outArray, const sunrealtype *inArray,
                        const vector<int> &lookup);
  inline void rotateToY(sunrealtype *outArray, const sunrealtype *inArray,
                        const vector<int> &lookup);
  inline void rotateToZ(sunrealtype *outArray, const sunrealtype *inArray,
                        const vector<int> &lookup);
  ///@}
public:
  /// ID of the LatticePatch, corresponds to process number
  // (required solely for debugging)
  int ID;
  /// N_Vector for saving field components u=(E,B) in lattice points
  N_Vector u;
  /// N_Vector for saving temporal derivatives of the field data
  N_Vector du;
  /// pointer to field data
  sunrealtype *uData;
  /// pointer to auxiliary data vector
  sunrealtype *uAuxData;
  /// pointer to time-derivative data
  sunrealtype *duData;
  ///@{
  /** pointer to halo data */
  sunrealtype *gCLData, *gCRData;
  ///@}
  /// pointer to spatial derivative data buffers
  array<sunrealtype *, 3> buffData;
  /// constructor setting up a default first lattice patch
  LatticePatch();
  /// destructor freeing parallel vectors
  ~LatticePatch();
  /// friend function for creating the patchwork slicing of the overall lattice
  friend int generatePatchwork(const Lattice &envelopeLattice,
                               LatticePatch &patchToMold, const int DLx,
                               const int DLy, const int DLz);
  /// function to get the discrete size of the LatticePatch
  // (0 direction corresponds to total)
  int discreteSize(int dir=0) const;
  /// function to get the origin of the patch
  sunrealtype origin(const int dir) const;
  /// function to get distance between points
  sunrealtype getDelta(const int dir) const;
  /// function to fill out the lookup tables
  // for translocation and de-translocation of data point
  void generateTranslocationLookup();
  /// function to rotate u into Z-matrix eigenraum
  // and make it the primary lattice direction of dir
  void rotateIntoEigen(const int dir);
  /// function to derotate uAux into dudata lattice direction of x
  void derotate(int dir, sunrealtype *buffOut);
  /// initialize ghost cells for halo exchange
  void initializeGhostLayer();
  /// initialize buffers to save derivatives
  void initializeBuffers();
  /// function to exchange ghost cells in uAux for the derivative
  void exchangeGhostCells(const int dir);
  /// function to derive the centered values in uAux and save them noncentered
  void derive(const int dir);
  /// function to check if a flag has been set and if not abort
  void checkFlag(unsigned int flag) const;
};

// helper function for error messages
void errorKill(const string & errorMessage);

// helper function to check for CVode success
int check_retval(void *returnvalue, const char *funcname, int opt, int id);

// End of Includeguard
#endif
