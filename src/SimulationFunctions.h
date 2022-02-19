//////////////////////////////////////////////////////////////////////////
/// @file SimulationFunctions.h
/// @brief Full simulation functions for 1D, 2D, and 3D used in main.cpp
//////////////////////////////////////////////////////////////////////////

// math
#include <cmath>
//#include <mathimf.h>

// project subfile headers
#include "LatticePatch.h"
#include "SimulationClass.h"
#include "TimeEvolutionFunctions.h"

/****** EM-wave structures ******/

/// plane wave structure
struct planewave {
  vector<sunrealtype> k;   /**< wavevector (normalized to \f$ 1/\lambda \f$) */
  vector<sunrealtype> p;   /**< amplitde & polarization vector */
  vector<sunrealtype> phi; /**< phase shift */
};

/// 1D Gaussian wave structure
struct gaussian1D {
  vector<sunrealtype> k;   /**< wavevector (normalized to \f$ 1/\lambda \f$) */
  vector<sunrealtype> p;   /**< amplitude & polarization vector */
  vector<sunrealtype> x0;  /**< shift from origin */
  sunrealtype phig;        /**< width */
  vector<sunrealtype> phi; /**< phase shift */
};

/// 2D Gaussian wave structure
struct gaussian2D {
  vector<sunrealtype> x0;   /**< center */
  vector<sunrealtype> axis; /**< direction to center */
  sunrealtype amp;          /**< amplitude */
  sunrealtype phip;         /**< polarization rotation */
  sunrealtype w0;           /**< taille */
  sunrealtype zr;           /**< Rayleigh length */
  sunrealtype ph0;          /**< beam center */
  sunrealtype phA;          /**< beam length */
};

/// 3D Gaussian wave structure
struct gaussian3D {
  vector<sunrealtype> x0;   /**< center */
  vector<sunrealtype> axis; /**< direction to center */
  sunrealtype amp;          /**< amplitude */
  sunrealtype phip;         /**< polarization rotation */
  sunrealtype w0;           /**< taille */
  sunrealtype zr;           /**< Rayleigh length */
  sunrealtype ph0;          /**< beam center */
  sunrealtype phA;          /**< beam length */
};

/****** simulation function declarations ******/

/// complete 1D Simulation function
void Sim1D(const array<sunrealtype,2>, const int, const sunrealtype,
        const sunindextype, const bool, int *, const sunrealtype, const int,
        const string, const int, const vector<planewave> &,
        const vector<gaussian1D> &);
/// complete 2D Simulation function
void Sim2D(const array<sunrealtype,2>, const int, const array<sunrealtype,2>,
        const array<sunindextype,2>, const array<int,2>, const bool, int *,
        const sunrealtype, const int, const string, const int,
        const vector<planewave> &, const vector<gaussian2D> &);
/// complete 3D Simulation function
void Sim3D(const array<sunrealtype,2>, const int, const array<sunrealtype,3>,
        const array<sunindextype,3>, const array<int,3>, const bool, int *,
        const sunrealtype, const int, const string, const int,
        const vector<planewave> &, const vector<gaussian3D> &);

/** MPI timer function */
void timer(double &, double &);
