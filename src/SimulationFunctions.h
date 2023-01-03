//////////////////////////////////////////////////////////////////////////
/// @file SimulationFunctions.h
/// @brief Full simulation functions for 1D, 2D, and 3D used in main.cpp
//////////////////////////////////////////////////////////////////////////

#pragma once

// math
#include <cmath>

// project subfile headers
#include "LatticePatch.h"
#include "SimulationClass.h"
#include "TimeEvolutionFunctions.h"

/****** EM-wave structures ******/

/// plane wave structure
struct planewave {
  std::array<sunrealtype, 3> k;   /**< wavevector (normalized to \f$ 1/\lambda \f$) */
  std::array<sunrealtype, 3> p;   /**< amplitde & polarization vector */
  std::array<sunrealtype, 3> phi; /**< phase shift */
};

/// 1D Gaussian wave structure
struct gaussian1D {
  std::array<sunrealtype, 3> k;   /**< wavevector (normalized to \f$ 1/\lambda \f$) */
  std::array<sunrealtype, 3> p;   /**< amplitude & polarization vector */
  std::array<sunrealtype, 3> x0;  /**< shift from origin */
  sunrealtype phig;        /**< width */
  std::array<sunrealtype, 3> phi; /**< phase shift */
};

/// 2D Gaussian wave structure
struct gaussian2D {
  std::array<sunrealtype, 3> x0;   /**< center */
  std::array<sunrealtype, 3> axis; /**< direction from where it comes */
  sunrealtype amp;          /**< amplitude */
  sunrealtype phip;         /**< polarization rotation */
  sunrealtype w0;           /**< taille */
  sunrealtype zr;           /**< Rayleigh length */
  sunrealtype ph0;          /**< beam center */
  sunrealtype phA;          /**< beam length */
};

/// 3D Gaussian wave structure
struct gaussian3D {
  std::array<sunrealtype, 3> x0;   /**< center */
  std::array<sunrealtype, 3> axis; /**< direction from where it comes */
  sunrealtype amp;          /**< amplitude */
  sunrealtype phip;         /**< polarization rotation */
  sunrealtype w0;           /**< taille */
  sunrealtype zr;           /**< Rayleigh length */
  sunrealtype ph0;          /**< beam center */
  sunrealtype phA;          /**< beam length */
};

/****** simulation function declarations ******/

/// complete 1D Simulation function
void Sim1D(const std::array<sunrealtype,2>, const int, const sunrealtype,
        const sunindextype, const bool, int *, const sunrealtype, const int,
        const std::string, const int, const char,
        const std::vector<planewave> &,
        const std::vector<gaussian1D> &);
/// complete 2D Simulation function
void Sim2D(const std::array<sunrealtype,2>, const int,
        const std::array<sunrealtype,2>,
        const std::array<sunindextype,2>, const std::array<int,2>,
        const bool, int *,
        const sunrealtype, const int, const std::string,
        const int, const char,
        const std::vector<planewave> &, const std::vector<gaussian2D> &);
/// complete 3D Simulation function
void Sim3D(const std::array<sunrealtype,2>, const int,
        const std::array<sunrealtype,3>,
        const std::array<sunindextype,3>, const std::array<int,3>,
        const bool, int *,
        const sunrealtype, const int, const std::string,
        const int, const char,
        const std::vector<planewave> &, const std::vector<gaussian3D> &);

/// timer function
void timer(double &, double &);
