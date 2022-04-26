////////////////////////////////////////////////////////
/// @file TimeEvolutionFunctions.h
/// @brief Functions to propagate data vectors in time
/// according to Maxwell's equations, and various
/// orders in the HE weak-field expansion
////////////////////////////////////////////////////////

#pragma once

#include "LatticePatch.h"
#include "SimulationClass.h"

/** @brief monostate TimeEvolution Class to propagate the field data in time in
 * a given order of the HE weak-field expansion */
class TimeEvolution {
public:
  /// choice which processes of the weak field expansion are included
  static int *c;

  /// Pointer to functions for differentiation and time evolution
  static void (*TimeEvolver)(LatticePatch *, N_Vector, N_Vector, int *);

  /// CVODE right hand side function (CVRhsFn) to provide IVP of the ODE
  static int f(sunrealtype t, N_Vector u, N_Vector udot, void *data_loc);
};

/// Maxwell propagation function for 1D -- only for reference
void linear1DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c);
/// HE propagation function for 1D
void nonlinear1DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c);
/// Maxwell propagation function for 2D -- only for reference
void linear2DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c);
/// HE propagation function for 2D
void nonlinear2DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c);
/// Maxwell propagation function for 3D -- only for reference
void linear3DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c);
/// HE propagation function for 3D
void nonlinear3DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c);

