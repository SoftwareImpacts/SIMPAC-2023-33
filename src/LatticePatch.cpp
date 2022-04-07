////////////////////////////////////////////////////////////////////////////////
/// @file LatticePatch.cpp
/// @brief Costruction of the overall envelope lattice and the lattice patches
////////////////////////////////////////////////////////////////////////////////

#include "LatticePatch.h"

#include <math.h>

///////////////////////////////////////////////////////
//// Implementation of Lattice component functions ////
///////////////////////////////////////////////////////

/// Initialize the cartesian communicator
void Lattice::initializeCommunicator(const int nx, const int ny,
        const int nz, const bool per) {
  const int dims[3] = {nz, ny, nx};
  const int periods[3] = {static_cast<int>(per), static_cast<int>(per),
                    static_cast<int>(per)};
  // Create the cartesian communicator for MPI_COMM_WORLD
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &comm);
  // Set MPI variables of the lattice
  MPI_Comm_size(comm, &(n_prc));
  MPI_Comm_rank(comm, &(my_prc));
  // Associate name to the communicator to identify it -> for debugging and
  // nicer error messages
  constexpr char lattice_comm_name[] = "Lattice";
  MPI_Comm_set_name(comm, lattice_comm_name);

  // Test if process naming is the same for both communicators
  /*
  int MYPRC;
  MPI_Comm_rank(MPI_COMM_WORLD,&MYPRC);
  cout<<"\r"<<my_prc<<"\t"<<MYPRC<<endl;
  */
}

/// Construct the lattice and set the stencil order
Lattice::Lattice(const int StO) : stencilOrder(StO),
    ghostLayerWidth(StO/2+1) {
  statusFlags = 0;
}

/// Set the number of points in each dimension of the lattice
void Lattice::setDiscreteDimensions(const sunindextype _nx,
        const sunindextype _ny, const sunindextype _nz) {
  // copy the given data for number of points
  tot_nx = _nx;
  tot_ny = _ny;
  tot_nz = _nz;
  // compute the resulting number of points and datapoints
  tot_noP = tot_nx * tot_ny * tot_nz;
  tot_noDP = dataPointDimension * tot_noP;
  // compute the new Delta, the physical resolution
  dx = tot_lx / tot_nx;
  dy = tot_ly / tot_ny;
  dz = tot_lz / tot_nz;
}

/// Set the physical size of the lattice
void Lattice::setPhysicalDimensions(const sunrealtype _lx,
        const sunrealtype _ly, const sunrealtype _lz) {
  tot_lx = _lx;
  tot_ly = _ly;
  tot_lz = _lz;
  // calculate physical distance between points
  dx = tot_lx / tot_nx;
  dy = tot_ly / tot_ny;
  dz = tot_lz / tot_nz;
  statusFlags |= FLatticeDimensionSet;
}

////////////////////////////////////////////////////////////
//// Implementation of LatticePatch component functions ////
////////////////////////////////////////////////////////////

/// Construct the lattice patch
LatticePatch::LatticePatch() {
  // set default origin coordinates to (0,0,0)
  x0 = y0 = z0 = 0;
  // set default position in Lattice-Patchwork to (0,0,0)
  LIx = LIy = LIz = 0;
  // set default physical lentgth for lattice patch to (0,0,0)
  lx = ly = lz = 0;
  // set default discrete length for lattice patch to (0,1,1)
  /* This is done in this manner as even in 1D simulations require a 1 point
   * width */
  nx = 0;
  ny = nz = 1;

  // u is not initialized as it wouldn't make any sense before the dimensions
  // are set idem for the enveloping lattice

  // set default statusFlags to non set
  statusFlags = 0;
}

/// Destruct the patch and thereby destroy the NVectors
LatticePatch::~LatticePatch() {
  // Deallocate memory for solution vector
  if (statusFlags & FLatticePatchSetUp) {
    // Destroy data vectors
    N_VDestroy_Parallel(u);
    N_VDestroy_Parallel(du);
  }
}

/// Set up the patchwork
int generatePatchwork(const Lattice &envelopeLattice, LatticePatch &patchToMold,
                      const int DLx, const int DLy, const int DLz) {
  // Retrieve the ghost layer depth
  const int gLW = envelopeLattice.get_ghostLayerWidth();
  // Retrieve the data point dimension
  const int dPD = envelopeLattice.get_dataPointDimension();
  // MPI process/patch
  const int my_prc = envelopeLattice.my_prc;
  // Determine thicknes of the slice
  const sunindextype tot_NOXP = envelopeLattice.get_tot_nx(); // total points of lattice
  const sunindextype tot_NOYP = envelopeLattice.get_tot_ny();
  const sunindextype tot_NOZP = envelopeLattice.get_tot_nz();
  // position of the patch in the lattice of patches – process associated to
  // position
  const sunindextype LIx = my_prc % DLx;
  const sunindextype LIy = (my_prc / DLx) % DLy;
  const sunindextype LIz = (my_prc / DLx) / DLy;
  // Determine the number of points in the patch and first absolute points in
  // each dimension
  const sunindextype local_NOXP = tot_NOXP / DLx;
  const sunindextype local_NOYP = tot_NOYP / DLy;
  const sunindextype local_NOZP = tot_NOZP / DLz;
  // absolute positions of the first point in each dimension
  const sunindextype firstXPoint = local_NOXP * LIx;
  const sunindextype firstYPoint = local_NOYP * LIy;
  const sunindextype firstZPoint = local_NOZP * LIz;
  // total number of points in the patch
  const sunindextype local_NODP = dPD * local_NOXP * local_NOYP * local_NOZP;

  // Set patch up with above derived quantities
  patchToMold.dx = envelopeLattice.get_dx();
  patchToMold.dy = envelopeLattice.get_dy();
  patchToMold.dz = envelopeLattice.get_dz();
  patchToMold.x0 = firstXPoint * patchToMold.dx;
  patchToMold.y0 = firstYPoint * patchToMold.dy;
  patchToMold.z0 = firstZPoint * patchToMold.dz;
  patchToMold.LIx = LIx;
  patchToMold.LIy = LIy;
  patchToMold.LIz = LIz;
  patchToMold.nx = local_NOXP;
  patchToMold.ny = local_NOYP;
  patchToMold.nz = local_NOZP;
  patchToMold.lx = patchToMold.nx * patchToMold.dx;
  patchToMold.ly = patchToMold.ny * patchToMold.dy;
  patchToMold.lz = patchToMold.nz * patchToMold.dz;
  /* Create and allocate memory for parallel vectors with defined local and
   * global lenghts *
   * (-> CVode problem sizes Nlocal and N)
   * for field data and temporal derivatives and set extra pointers to them */
  patchToMold.u =
      N_VNew_Parallel(envelopeLattice.comm, local_NODP,
                      envelopeLattice.get_tot_noDP(), envelopeLattice.sunctx);
  patchToMold.uData = NV_DATA_P(patchToMold.u);
  patchToMold.du =
      N_VNew_Parallel(envelopeLattice.comm, local_NODP,
                      envelopeLattice.get_tot_noDP(), envelopeLattice.sunctx);
  patchToMold.duData = NV_DATA_P(patchToMold.du);
  // Allocate space for auxiliary uAux so that the lattice and all possible
  // directions of ghost Layers fit
  const int s1 = patchToMold.nx, s2 = patchToMold.ny, s3 = patchToMold.nz;
  const int s_min = min(s1, min(s2, s3));
  patchToMold.uAux.resize(s1 * s2 * s3 / s_min * (s_min + 2 * gLW) * dPD);
  patchToMold.uAuxData = &patchToMold.uAux[0];
  patchToMold.envelopeLattice = &envelopeLattice;
  // Set patch "name" to process number -> only for debugging
  // patchToMold.ID=my_prc;
  // set flag
  patchToMold.statusFlags = FLatticePatchSetUp;
  patchToMold.generateTranslocationLookup();
  return 0;
}

/// Return the discrete size of the patch: number of lattice patch points in
/// specified dimension
int LatticePatch::discreteSize(int dir) const {
  switch (dir) {
  case 0:
    return nx * ny * nz;
  case 1:
    return nx;
  case 2:
    return ny;
  case 3:
    return nz;
  // case 4: return uAux.size(); // for debugging
  default:
    return -1;
  }
}

/// Return the physical origin of the patch in a dimension
sunrealtype LatticePatch::origin(const int dir) const {
  switch (dir) {
  case 1:
    return x0;
  case 2:
    return y0;
  case 3:
    return z0;
  default:
    errorKill("LatticePatch::origin function called with wrong dir parameter");
    return -1;
  }
}

/// Return the distance between points in the patch in a dimension
sunrealtype LatticePatch::getDelta(const int dir) const {
  switch (dir) {
  case 1:
    return dx;
  case 2:
    return dy;
  case 3:
    return dz;
  default:
    errorKill(
        "LatticePatch::getDelta function called with wrong dir parameter");
    return -1;
  }
}

/** To avoid cache misses:
 * create vectors to translate u vector into space coordinates and vice versa
 * and same for left and right ghost layers to space */
void LatticePatch::generateTranslocationLookup() {
  // Check that the lattice has been set up
  checkFlag(FLatticeDimensionSet);
  // lenghts for auxilliary layers, including ghost layers
  const int gLW = envelopeLattice->get_ghostLayerWidth();
  const int mx = nx + 2 * gLW;
  const int my = ny + 2 * gLW;
  const int mz = nz + 2 * gLW;
  // sizes for lookup vectors
  // generate u->uAux
  uTox.resize(nx * ny * nz);
  uToy.resize(nx * ny * nz);
  uToz.resize(nx * ny * nz);
  // generate uAux->u with length including halo
  xTou.resize(mx * ny * nz);
  yTou.resize(nx * my * nz);
  zTou.resize(nx * ny * mz);
  // variables for cartesian position in the 3D discrete lattice
  int px = 0, py = 0, pz = 0;
  for (int i = 0; i < uToy.size(); i++) { // loop over all points in the patch
    // calulate cartesian coordinates
    px = i % nx;
    py = (i / nx) % ny;
    pz = (i / nx) / ny;
    // fill lookups extended by halos (useful for y and z direction)
    uTox[i] = (px + gLW) + py * mx +
              pz * mx * ny; // unroll (de-flatten) cartesian dimension
    xTou[px + py * mx + pz * mx * ny] =
        i; // match cartesian point to u location
    uToy[i] = (py + gLW) + pz * my + px * my * nz;
    yTou[py + pz * my + px * my * nz] = i;
    uToz[i] = (pz + gLW) + px * mz + py * mz * nx;
    zTou[pz + px * mz + py * mz * nx] = i;
  }
  // same for ghost layer lookup tables
  lgcTox.resize(gLW * ny * nz);
  rgcTox.resize(gLW * ny * nz);
  for (int i = 0; i < lgcTox.size(); i++) {
    px = i % gLW;
    py = (i / gLW) % ny;
    pz = (i / gLW) / ny;
    lgcTox[i] = px + py * mx + pz * mx * ny;
    rgcTox[i] = px + nx + gLW + py * mx + pz * mx * ny;
  }
  lgcToy.resize(gLW * nx * nz);
  rgcToy.resize(gLW * nx * nz);
  for (int i = 0; i < lgcToy.size(); i++) {
    px = i % nx;
    py = (i / nx) % gLW;
    pz = (i / nx) / gLW;
    lgcToy[i] = py + pz * my + px * my * nz;
    rgcToy[i] = py + ny + gLW + pz * my + px * my * nz;
  }
  lgcToz.resize(gLW * nx * ny);
  rgcToz.resize(gLW * nx * ny);
  for (int i = 0; i < lgcToz.size(); i++) {
    px = i % nx;
    py = (i / nx) % ny;
    pz = (i / nx) / ny;
    lgcToz[i] = pz + px * mz + py * mz * nx;
    rgcToz[i] = pz + nz + gLW + px * mz + py * mz * nx;
  }
  statusFlags |= TranslocationLookupSetUp;
}

/** Rotate into eigenraum along R matrices of paper using below rotation
 * functions
 * -> uAuxData gets the rotated left-halo-, inner-patch-, right-halo-data */
void LatticePatch::rotateIntoEigen(const int dir) {
  // Check that the lattice, ghost layers as well as the translocation lookups
  // have been set up;
  checkFlag(FLatticePatchSetUp);
  checkFlag(TranslocationLookupSetUp);
  checkFlag(GhostLayersInitialized); // this check is only after call to
                                     // exchange ghost cells
  switch (dir) {
  case 1:
    rotateToX(uAuxData, gCLData, lgcTox);
    rotateToX(uAuxData, uData, uTox);
    rotateToX(uAuxData, gCRData, rgcTox);
    break;
  case 2:
    rotateToY(uAuxData, gCLData, lgcToy);
    rotateToY(uAuxData, uData, uToy);
    rotateToY(uAuxData, gCRData, rgcToy);
    break;
  case 3:
    rotateToZ(uAuxData, gCLData, lgcToz);
    rotateToZ(uAuxData, uData, uToz);
    rotateToZ(uAuxData, gCRData, rgcToz);
    break;
  default:
    errorKill("Tried to rotate into the wrong direction");
    break;
  }
}

/// Rotate halo and inner-patch data vectors with rotation matrix Rx into
/// eigenspace of Z matrix and write to auxiliary vector
inline void LatticePatch::rotateToX(sunrealtype *outArray,
                                    const sunrealtype *inArray,
                                    const vector<int> &lookup) {
  int ii = 0, target = 0;
#pragma omp simd // safelen(6)
  for (int i = 0; i < lookup.size(); i++) {
    // get correct u-vector and spatial indices along previously defined lookup
    // tables
    target = envelopeLattice->get_dataPointDimension() * lookup[i];
    ii = envelopeLattice->get_dataPointDimension() * i;
    outArray[target + 0] = -inArray[1 + ii] + inArray[5 + ii];
    outArray[target + 1] = inArray[2 + ii] + inArray[4 + ii];
    outArray[target + 2] = inArray[1 + ii] + inArray[5 + ii];
    outArray[target + 3] = -inArray[2 + ii] + inArray[4 + ii];
    outArray[target + 4] = inArray[3 + ii];
    outArray[target + 5] = inArray[ii];
  }
}

/// Rotate halo and inner-patch data vectors with rotation matrix Ry into
/// eigenspace of Z matrix and write to auxiliary vector
inline void LatticePatch::rotateToY(sunrealtype *outArray,
                                    const sunrealtype *inArray,
                                    const vector<int> &lookup) {
  int ii = 0, target = 0;
#pragma omp simd
  for (int i = 0; i < lookup.size(); i++) {
    target = envelopeLattice->get_dataPointDimension() * lookup[i];
    ii = envelopeLattice->get_dataPointDimension() * i;
    outArray[target + 0] = inArray[ii] + inArray[5 + ii];
    outArray[target + 1] = -inArray[2 + ii] + inArray[3 + ii];
    outArray[target + 2] = -inArray[ii] + inArray[5 + ii];
    outArray[target + 3] = inArray[2 + ii] + inArray[3 + ii];
    outArray[target + 4] = inArray[4 + ii];
    outArray[target + 5] = inArray[1 + ii];
  }
}

/// Rotate halo and inner-patch data vectors with rotation matrix Rz into
/// eigenspace of Z matrix and write to auxiliary vector
inline void LatticePatch::rotateToZ(sunrealtype *outArray,
                                    const sunrealtype *inArray,
                                    const vector<int> &lookup) {
  int ii = 0, target = 0;
#pragma omp simd
  for (int i = 0; i < lookup.size(); i++) {
    target = envelopeLattice->get_dataPointDimension() * lookup[i];
    ii = envelopeLattice->get_dataPointDimension() * i;
    outArray[target + 0] = -inArray[ii] + inArray[4 + ii];
    outArray[target + 1] = inArray[1 + ii] + inArray[3 + ii];
    outArray[target + 2] = inArray[ii] + inArray[4 + ii];
    outArray[target + 3] = -inArray[1 + ii] + inArray[3 + ii];
    outArray[target + 4] = inArray[5 + ii];
    outArray[target + 5] = inArray[2 + ii];
  }
}

/// Derotate uAux with transposed rotation matrices and write to derivative
/// buffer – normalization is done here by the factor 1/2
void LatticePatch::derotate(int dir, sunrealtype *buffOut) {
  // Check that the lattice as well as the translocation lookups have been set
  // up;
  checkFlag(FLatticePatchSetUp);
  checkFlag(TranslocationLookupSetUp);
  const int dPD = envelopeLattice->get_dataPointDimension();
  const int gLW = envelopeLattice->get_ghostLayerWidth();
  const int uSize = discreteSize();
  int ii = 0, target = 0;
  switch (dir) {
  case 1:
#pragma omp simd
    for (int i = 0; i < uSize; i++) {
      // get correct indices in u and rotation space
      target = dPD * i;
      ii = dPD * (uTox[i] - gLW);
      buffOut[target + 0] = uAux[5 + ii];
      buffOut[target + 1] = (-uAux[ii] + uAux[2 + ii]) / 2.;
      buffOut[target + 2] = (uAux[1 + ii] - uAux[3 + ii]) / 2.;
      buffOut[target + 3] = uAux[4 + ii];
      buffOut[target + 4] = (uAux[1 + ii] + uAux[3 + ii]) / 2.;
      buffOut[target + 5] = (uAux[ii] + uAux[2 + ii]) / 2.;
    }
    break;
  case 2:
#pragma omp simd
    for (int i = 0; i < uSize; i++) {
      target = dPD * i;
      ii = dPD * (uToy[i] - gLW);
      buffOut[target + 0] = (uAux[ii] - uAux[2 + ii]) / 2.;
      buffOut[target + 1] = uAux[5 + ii];
      buffOut[target + 2] = (-uAux[1 + ii] + uAux[3 + ii]) / 2.;
      buffOut[target + 3] = (uAux[1 + ii] + uAux[3 + ii]) / 2.;
      buffOut[target + 4] = uAux[4 + ii];
      buffOut[target + 5] = (uAux[ii] + uAux[2 + ii]) / 2.;
    }
    break;
  case 3:
#pragma omp simd
    for (int i = 0; i < uSize; i++) {
      target = dPD * i;
      ii = dPD * (uToz[i] - gLW);
      buffOut[target + 0] = (-uAux[ii] + uAux[2 + ii]) / 2.;
      buffOut[target + 1] = (uAux[1 + ii] - uAux[3 + ii]) / 2.;
      buffOut[target + 2] = uAux[5 + ii];
      buffOut[target + 3] = (uAux[1 + ii] + uAux[3 + ii]) / 2.;
      buffOut[target + 4] = (uAux[ii] + uAux[2 + ii]) / 2.;
      buffOut[target + 5] = uAux[4 + ii];
    }
    break;
  default:
    errorKill("Tried to derotate from the wrong direction");
    break;
  }
}

/// Create buffers to save derivative values, optimizing computational load
void LatticePatch::initializeBuffers() {
  // Check that the lattice has been set up
  checkFlag(FLatticeDimensionSet);
  const int dPD = envelopeLattice->get_dataPointDimension();
  buffX.resize(nx * ny * nz * dPD);
  buffY.resize(nx * ny * nz * dPD);
  buffZ.resize(nx * ny * nz * dPD);
  // Set pointers used for propagation functions
  buffData[0] = &buffX[0];
  buffData[1] = &buffY[0];
  buffData[2] = &buffZ[0];
  statusFlags |= BuffersInitialized;
}

/// Perform the ghost cell exchange in a specified direction
void LatticePatch::exchangeGhostCells(const int dir) {
  // Check that the lattice has been set up
  checkFlag(FLatticeDimensionSet);
  checkFlag(FLatticePatchSetUp);
  // Variables to per dimension calculate the halo indices, and distance to
  // other side halo boundary
  int mx = 1, my = 1, mz = 1, distToRight = 1;
  const int gLW = envelopeLattice->get_ghostLayerWidth();
  // In the chosen direction m is set to ghost layer width while the others
  // remain to form the plane
  switch (dir) {
  case 1:
    mx = gLW;
    my = ny;
    mz = nz;
    distToRight = (nx - gLW);
    break;
  case 2:
    mx = nx;
    my = gLW;
    mz = nz;
    distToRight = nx * (ny - gLW);
    break;
  case 3:
    mx = nx;
    my = ny;
    mz = gLW;
    distToRight = nx * ny * (nz - gLW);
    break;
  }
  // total number of exchanged points
  const int dPD = envelopeLattice->get_dataPointDimension();
  const int exchangeSize = mx * my * mz * dPD;
  // provide size of the halos for ghost cells
  ghostCellLeft.resize(exchangeSize);
  ghostCellRight.resize(ghostCellLeft.size());
  ghostCellLeftToSend.resize(ghostCellLeft.size());
  ghostCellRightToSend.resize(ghostCellLeft.size());
  gCLData = &ghostCellLeft[0];
  gCRData = &ghostCellRight[0];
  statusFlags |= GhostLayersInitialized;

  // Initialize running index li for the halo buffers, and index ui of uData for
  // data transfer
  int li = 0, ui = 0;

  for (int iz = 0; iz < mz; iz++) {
    for (int iy = 0; iy < my; iy++) {
      // uData vector start index of halo data to be transferred
      // with each z-step add the whole xy-plane and with y-step the x-range ->
      // iterate all x-ranges
      ui = (iz * nx * ny + iy * nx) * dPD;
      // copy left halo data from uData to buffer, transfer size is given by
      // x-length (not x-range) perhaps faster but more fragile C lib copy
      // operation (contained in cstring header)
      /*
      memcpy(&ghostCellLeftToSend[li],
             &uData[ui],
             sizeof(sunrealtype)*mx*dPD);
      // increase ui by the distance to vis-a-vis boundary and copy right halo
      data to buffer ui+=distToRight*dPD; memcpy(&ghostCellRightToSend[li],
             &uData[ui],
             sizeof(sunrealtype)*mx*dPD);
      */
      // perhaps more safe but slower copy operation (contained in algorithm
      // header) performance highly system dependent
      copy(&uData[ui], &uData[ui + mx * dPD], &ghostCellLeftToSend[li]);
      ui += distToRight * dPD;
      copy(&uData[ui], &uData[ui + mx * dPD], &ghostCellRightToSend[li]);

      // increase halo index by transferred items per y-iteration step
      // (x-length)
      li += mx * dPD;
    }
  }

  /* Send and receive the data to and from neighboring latticePatches */
  // Adjust direction to cartesian communicator
  int dim = 2; // default for dir==1
  if (dir == 2) {
    dim = 1;
  } else if (dir == 3) {
    dim = 0;
  }
  int rank_source = 0, rank_dest = 0;
  MPI_Cart_shift(envelopeLattice->comm, dim, -1, &rank_source,
                 &rank_dest); // s.t. rank_dest is left & v.v.

  // nonblocking Irecv/Isend
  
  MPI_Request requests[4];
  MPI_Irecv(&ghostCellRight[0], exchangeSize, MPI_SUNREALTYPE, rank_source, 1,
  envelopeLattice->comm, &requests[0]);
  MPI_Isend(&ghostCellLeftToSend[0], exchangeSize, MPI_SUNREALTYPE, rank_dest,
  1, envelopeLattice->comm, &requests[1]);
  MPI_Irecv(&ghostCellLeft[0], exchangeSize, MPI_SUNREALTYPE, rank_dest, 2,
  envelopeLattice->comm, &requests[2]);
  MPI_Isend(&ghostCellRightToSend[0], exchangeSize, MPI_SUNREALTYPE,
  rank_source, 2, envelopeLattice->comm, &requests[3]);
  MPI_Waitall(4, requests, MPI_STATUS_IGNORE);
  

  // blocking Sendrecv:
  /*
  MPI_Sendrecv(&ghostCellLeftToSend[0], exchangeSize, MPI_SUNREALTYPE,
               rank_dest, 1, &ghostCellRight[0], exchangeSize, MPI_SUNREALTYPE,
               rank_source, 1, envelopeLattice->comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&ghostCellRightToSend[0], exchangeSize, MPI_SUNREALTYPE,
               rank_source, 2, &ghostCellLeft[0], exchangeSize, MPI_SUNREALTYPE,
               rank_dest, 2, envelopeLattice->comm, MPI_STATUS_IGNORE);
  */
}

/// Check if all flags are set
void LatticePatch::checkFlag(unsigned int flag) const {
  if (!(statusFlags & flag)) {
    string errorMessage;
    switch (flag) {
    case FLatticePatchSetUp:
      errorMessage = "The Lattice patch was not set up please make sure to "
                     "initilize a Lattice topology";
      break;
    case TranslocationLookupSetUp:
      errorMessage = "The translocation lookup tables have not been generated, "
                     "please be sure to run generateTranslocationLookup()";
      break;
    case GhostLayersInitialized:
      errorMessage = "The space for the ghost layers has not been allocated, "
                     "please be sure to run initializeGhostLayer()";
      break;
    case BuffersInitialized:
      errorMessage = "The space for the buffers has not been allocated, please "
                     "be sure to run initializeBuffers()";
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

/// Calculate derivatives in the patch (uAux) in the specified direction
void LatticePatch::derive(const int dir) {
  // ghost layer width
  const int gLW = envelopeLattice->get_ghostLayerWidth();
  // dimensionality of data points -> 6
  const int dPD = envelopeLattice->get_dataPointDimension();
  // total width of patch in given direction including ghost layers at ends
  const int dirWidth = discreteSize(dir) + 2 * gLW;
  // width of patch only in given direction
  const int dirWidthO = discreteSize(dir);
  // size of plane perpendicular to given dimension
  const int perpPlainSize = discreteSize() / discreteSize(dir);
  // physical distance between points in that direction
  sunrealtype dxi = NAN;
  switch (dir) {
  case 1:
    dxi = dx;
    break;
  case 2:
    dxi = dy;
    break;
  case 3:
    dxi = dz;
    break;
  default:
    dxi = 1;
    errorKill("Tried to derive in the wrong direction");
    break;
  }
  // Derive according to chosen stencil accuracy order (which determines also
  // gLW)
  const int order = envelopeLattice->get_stencilOrder();
  switch (order) {
  case 1:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s1b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s1b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s1f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s1f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s1f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s1f(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 2:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s2b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s2b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s2f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s2f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s2c(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s2c(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 3:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s3b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s3b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s3f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s3f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s3f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s3f(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 4:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s4b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s4b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s4f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s4f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s4c(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s4c(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 5:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s5b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s5b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s5f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s5f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s5f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s5f(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 6:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s6b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s6b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s6f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s6f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s6c(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s6c(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 7:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s7b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s7b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s7f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s7f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s7f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s7f(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 8:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s8b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s8b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s8f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s8f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s8c(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s8c(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 9:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s9b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s9b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s9f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s9f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s9f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s9f(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 10:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s10b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s10b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s10f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s10f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s10c(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s10c(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 11:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s11b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s11b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s11f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s11f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s11f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s11f(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 12:
    for (int i = 0; i < perpPlainSize; i++) {
      for (int j = (i * dirWidth + gLW) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        uAux[j + 0 - gLW * dPD] = s12b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s12b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s12f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s12f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s12c(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s12c(&uAux[j + 5]) / dxi;
      }
    }
    break;
  case 13:
    // Iterate through all points in the plane perpendicular to the given
    // direction
    for (int i = 0; i < perpPlainSize; i++) {
      // Iterate through the direction for each perpendicular plane point
      for (int j = (i * dirWidth + gLW /*to shift left by gLW below */) * dPD;
           j < (i * dirWidth + gLW + dirWidthO) * dPD; j += dPD) {
        /* Compute the stencil derivative for any of the six field components
         * with a ghostlayer width adjusted to the order of the finite
         * difference scheme */
        uAux[j + 0 - gLW * dPD] = s13b(&uAux[j + 0]) / dxi;
        uAux[j + 1 - gLW * dPD] = s13b(&uAux[j + 1]) / dxi;
        uAux[j + 2 - gLW * dPD] = s13f(&uAux[j + 2]) / dxi;
        uAux[j + 3 - gLW * dPD] = s13f(&uAux[j + 3]) / dxi;
        uAux[j + 4 - gLW * dPD] = s13f(&uAux[j + 4]) / dxi;
        uAux[j + 5 - gLW * dPD] = s13f(&uAux[j + 5]) / dxi;
      }
    }
    break;

  default:
    errorKill("Please set an existing stencil order");
    break;
  }
}

///////// Helper functions /////////

/// Print a specific error message to stdout
void errorKill(const string & errorMessage) {
  cerr << endl << "Error: " << errorMessage << " Aborting..." << endl;
  MPI_Abort(MPI_COMM_WORLD, 1);
  return;
}

/** Check function return value. From CVode examples.
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
int check_retval(void *returnvalue, const char *funcname, int opt, int id) {
  int *retval = nullptr;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == nullptr) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n", id,
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *)returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with retval = %d\n\n",
              id, funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == nullptr) {
    fprintf(stderr,
            "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n", id,
            funcname);
    return (1);
  }

  return (0);
}
