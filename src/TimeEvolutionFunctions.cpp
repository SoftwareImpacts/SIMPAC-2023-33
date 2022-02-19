//////////////////////////////////////////////////////////
/// @file TimeEvolutionFunctions.cpp
/// @brief Implementation of functions to propagate
/// data vectors in time according to Maxwell's equations,
/// and various orders in the HE weak-field expansion
//////////////////////////////////////////////////////////

#include "TimeEvolutionFunctions.h"

#include <math.h>

/// CVode right-hand-side function (CVRhsFn)
int TimeEvolution::f(sunrealtype t, N_Vector u, N_Vector udot, void *data_loc) {
  // Set recover pointer to provided lattice patch where the data resides
  LatticePatch *data = nullptr;
  data = static_cast<LatticePatch *>(data_loc);

  // pointers for update circle
  sunrealtype *udata = nullptr, *dudata = nullptr;
  sunrealtype *originaluData = nullptr, *originalduData = nullptr;

  // Access NVECTOR_PARALLEL argument data with pointers
  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);

  // Store original data location of the patch
  originaluData = data->uData;
  originalduData = data->duData;
  // Point patch data to arguments of f
  data->uData = udata;
  data->duData = dudata;

  // Time-evolve these arguments (the field data) with specific propagator below
  TimeEvolver(data, u, udot, c);

  // Refer patch data back to original location
  data->uData = originaluData;
  data->duData = originalduData;

  return (0);
}

/// only under-the-hood-callable Maxwell propagation in 1D
// unused parameters 2-4 for compliance with CVRhsFn
// same as the respective nonlinear function without nonlinear terms
void linear1DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c) {

  // pointers to temporal and spatial derivative data
  sunrealtype *duData = data->duData;
  sunrealtype *dxData = data->buffData[1 - 1];

  // sequence along any dimension:
  data->exchangeGhostCells(1); // exchange halos
  data->rotateIntoEigen(
      1);          // -> rotate all data to prepare derivative operation
  data->derive(1); // -> perform derivative on it
  data->derotate(
      1, dxData); // -> derotate derivative data to x-space for further use

  int totalNP = data->discreteSize();
  int pp = 0;
  for (int i = 0; i < totalNP; i++) {
    pp = i * 6;
    /*
     simple vacuum Maxwell equations for spatial deriative only in x-direction
     temporal derivative is approximated by spatial derivative according to the
     numerical scheme with Jacobi=0 -> no polarization or magnetization terms
    */
    duData[pp + 0] = 0;
    duData[pp + 1] = -dxData[pp + 5];
    duData[pp + 2] = dxData[pp + 4];
    duData[pp + 3] = 0;
    duData[pp + 4] = dxData[pp + 2];
    duData[pp + 5] = -dxData[pp + 1];
  }
}

/// nonlinear 1D HE propagation
void nonlinear1DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c) {

  // pointer to spatial derivative data sufficient, temporal derivative data
  // provided with udot
  sunrealtype *dxData = data->buffData[1 - 1];

  // same sequence as in the linear case
  data->exchangeGhostCells(1);
  data->rotateIntoEigen(1);
  data->derive(1);
  data->derotate(1, dxData);

  /*
  F and G are nonzero in the nonlinear case,
  polarization and magnetization contributions in Jacobi matrix style
  with derivatives of polarization and magnetization
  w.r.t. E- and B-field
  */
  sunrealtype f = NAN, g = NAN; // em field invariants F, G
  sunrealtype lf = NAN, lff = NAN, lfg = NAN, lg = NAN,
              lgg = NAN; // derivatives of Lagrangian w.r.t. field invariants
  array<sunrealtype, 36> JMM; // Jacobi matrix
  array<sunrealtype, 6> Quad; // array to hold E^2 and B^2 components
  array<sunrealtype, 6> h; // holding temporal derivatives of E and B components
                           // before operating (1+Z)^-1
  sunrealtype pseudoDenom = NAN; // needed for inversion of 1+Z
  sunrealtype *udata = nullptr,
              *dudata = nullptr; // pointers to data and temp. derivative data
  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);
  int totalNP = data->discreteSize(); // number of points in the patch
  for (int pp = 0; pp < totalNP * 6;
       pp += 6) { // loops through all 6dim points in the patch
                  //    for(int ppB=0;ppB<totalNP*6;ppB+=6*6){
                  //      for(int pp=ppB;pp<min(totalNP*6,ppB+6*6);pp+=6){
    /// Calculation of the Jacobi matrix
    // 1. Calculate F and G
    f = 0.5 * ((Quad[0] = udata[pp] * udata[pp]) +
               (Quad[1] = udata[pp + 1] * udata[pp + 1]) +
               (Quad[2] = udata[pp + 2] * udata[pp + 2]) -
               (Quad[3] = udata[pp + 3] * udata[pp + 3]) -
               (Quad[4] = udata[pp + 4] * udata[pp + 4]) -
               (Quad[5] = udata[pp + 5] * udata[pp + 5]));
    g = udata[pp] * udata[pp + 3] + udata[pp + 1] * udata[pp + 4] +
        udata[pp + 2] * udata[pp + 5];
    // 2. Choose process/expansion order and assign derivative values of L
    // w.r.t. F, G
    switch (*c) {
    case 0:
      lf = 0;
      lff = 0;
      lfg = 0;
      lg = 0;
      lgg = 0;
      break;
    case 2:
      lf = 0.000354046449700427580438254 * f * f +
           0.000191775160254398272737387 * g * g;
      lff = 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = 0.0003835503205087965454747749 * f * g;
      lgg = 0.0003835503205087965454747749 * f;
      break;
    case 1:
      lf = 0.000206527095658582755255648 * f;
      lff = 0.000206527095658582755255648;
      lfg = 0;
      lg = 0.0003614224174025198216973841 * g;
      lgg = 0.0003614224174025198216973841;
      break;
    case 3:
      lf = (0.000206527095658582755255648 + 0.000354046449700427580438254 * f) *
               f +
           0.000191775160254398272737387 * g * g;
      lff = 0.000206527095658582755255648 + 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = (0.0003614224174025198216973841 +
            0.0003835503205087965454747749 * f) *
           g;
      lgg = 0.0003614224174025198216973841 + 0.0003835503205087965454747749 * f;
      break;
    default:
      errorKill(
          "You need to specify a correct order in the weak-field expansion.");
    }
    // 3. Assign Jacobi components
    JMM[0] = lf + lff * Quad[0] +
             udata[3 + pp] * (2 * lfg * udata[pp] + lgg * udata[3 + pp]);
    JMM[6] =
        lff * udata[pp] * udata[1 + pp] + lfg * udata[1 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[4 + pp] + lgg * udata[3 + pp] * udata[4 + pp];
    JMM[7] = lf + lff * Quad[1] +
             udata[4 + pp] * (2 * lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[12] =
        lff * udata[pp] * udata[2 + pp] + lfg * udata[2 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[5 + pp] + lgg * udata[3 + pp] * udata[5 + pp];
    JMM[13] = lff * udata[1 + pp] * udata[2 + pp] +
              lfg * udata[2 + pp] * udata[4 + pp] +
              lfg * udata[1 + pp] * udata[5 + pp] +
              lgg * udata[4 + pp] * udata[5 + pp];
    JMM[14] = lf + lff * Quad[2] +
              udata[5 + pp] * (2 * lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[18] = lg + lfg * (Quad[0] - Quad[3 + 0]) +
              (-lff + lgg) * udata[pp] * udata[3 + pp];
    JMM[19] = -(udata[3 + pp] * (lff * udata[1 + pp] + lfg * udata[4 + pp])) +
              udata[pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[20] = -(udata[3 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[21] = -lf + lgg * Quad[0] +
              udata[3 + pp] * (-2 * lfg * udata[pp] + lff * udata[3 + pp]);
    JMM[24] = udata[1 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[4 + pp];
    JMM[25] = lg + lfg * (Quad[1] - Quad[4 + 0]) +
              (-lff + lgg) * udata[1 + pp] * udata[4 + pp];
    JMM[26] = -(udata[4 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[1 + pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[27] = lgg * udata[pp] * udata[1 + pp] +
              lff * udata[3 + pp] * udata[4 + pp] -
              lfg * (udata[1 + pp] * udata[3 + pp] + udata[pp] * udata[4 + pp]);
    JMM[28] = -lf + lgg * Quad[1] +
              udata[4 + pp] * (-2 * lfg * udata[1 + pp] + lff * udata[4 + pp]);
    JMM[30] = udata[2 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[5 + pp];
    JMM[31] = udata[2 + pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]) -
              (lff * udata[1 + pp] + lfg * udata[4 + pp]) * udata[5 + pp];
    JMM[32] = lg + lfg * (Quad[2] - Quad[5 + 0]) +
              (-lff + lgg) * udata[2 + pp] * udata[5 + pp];
    JMM[33] = lgg * udata[pp] * udata[2 + pp] +
              lff * udata[3 + pp] * udata[5 + pp] -
              lfg * (udata[2 + pp] * udata[3 + pp] + udata[pp] * udata[5 + pp]);
    JMM[34] =
        lgg * udata[1 + pp] * udata[2 + pp] +
        lff * udata[4 + pp] * udata[5 + pp] -
        lfg * (udata[2 + pp] * udata[4 + pp] + udata[1 + pp] * udata[5 + pp]);
    JMM[35] = -lf + lgg * Quad[2] +
              udata[5 + pp] * (-2 * lfg * udata[2 + pp] + lff * udata[5 + pp]);
    for (int i = 0; i < 6; i++) {
      for (int j = i + 1; j < 6; j++) {
        JMM[i * 6 + j] = JMM[j * 6 + i];
      }
    }
    // 4. Final values for temporal derivatives of field values
    h[0] = 0;
    h[1] = dxData[pp] * JMM[30] + dxData[1 + pp] * JMM[31] +
           dxData[2 + pp] * JMM[32] + dxData[3 + pp] * JMM[33] +
           dxData[4 + pp] * JMM[34] + dxData[5 + pp] * (-1 + JMM[35]);
    h[2] = -(dxData[pp] * JMM[24]) - dxData[1 + pp] * JMM[25] -
           dxData[2 + pp] * JMM[26] - dxData[3 + pp] * JMM[27] +
           dxData[4 + pp] * (1 - JMM[28]) - dxData[5 + pp] * JMM[29];
    h[3] = 0;
    h[4] = dxData[2 + pp];
    h[5] = -dxData[1 + pp];
    h[0] -= h[3] * JMM[3] + h[4] * JMM[4] + h[5] * JMM[5];
    h[1] -= h[3] * JMM[9] + h[4] * JMM[10] + h[5] * JMM[11];
    h[2] -= h[3] * JMM[15] + h[4] * JMM[16] + h[5] * JMM[17];
    // (1+Z)^-1 applies only to E components
    dudata[pp + 0] =
        h[2] * (-(JMM[2] * (1 + JMM[7])) + JMM[1] * JMM[8]) +
        h[1] * (JMM[2] * JMM[13] - JMM[1] * (1 + JMM[14])) +
        h[0] * (1 - JMM[8] * JMM[13] + JMM[14] + JMM[7] * (1 + JMM[14]));
    dudata[pp + 1] =
        h[2] * (JMM[2] * JMM[6] - (1 + JMM[0]) * JMM[8]) +
        h[1] * (1 - JMM[2] * JMM[12] + JMM[14] + JMM[0] * (1 + JMM[14])) +
        h[0] * (JMM[8] * JMM[12] - JMM[6] * (1 + JMM[14]));
    dudata[pp + 2] =
        h[2] * (1 - JMM[1] * JMM[6] + JMM[7] + JMM[0] * (1 + JMM[7])) +
        h[1] * (JMM[1] * JMM[12] - (1 + JMM[0]) * JMM[13]) +
        h[0] * (-((1 + JMM[7]) * JMM[12]) + JMM[6] * JMM[13]);
    pseudoDenom =
        -((1 + JMM[7]) * (-1 + JMM[2] * JMM[12])) +
        (JMM[2] * JMM[6] - JMM[8]) * JMM[13] + JMM[14] + JMM[7] * JMM[14] +
        JMM[0] * (1 + JMM[7] - JMM[8] * JMM[13] + (1 + JMM[7]) * JMM[14]) -
        JMM[1] * (-(JMM[8] * JMM[12]) + JMM[6] * (1 + JMM[14]));
    dudata[pp + 0] /= pseudoDenom;
    dudata[pp + 1] /= pseudoDenom;
    dudata[pp + 2] /= pseudoDenom;
    dudata[pp + 3] = h[3];
    dudata[pp + 4] = h[4];
    dudata[pp + 5] = h[5];
  }
  return;
}

/// only under-the-hood-callable Maxwell propagation in 2D
// unused parameters 2-4 for compliance with CVRhsFn
// same as the respective nonlinear function without nonlinear terms
void linear2DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c) {

  sunrealtype *duData = data->duData;
  sunrealtype *dxData = data->buffData[1 - 1];
  sunrealtype *dyData = data->buffData[2 - 1];

  data->exchangeGhostCells(1);
  data->rotateIntoEigen(1);
  data->derive(1);
  data->derotate(1, dxData);
  data->exchangeGhostCells(2);
  data->rotateIntoEigen(2);
  data->derive(2);
  data->derotate(2, dyData);

  int totalNP = data->discreteSize();
  int pp = 0;
  for (int i = 0; i < totalNP; i++) {
    pp = i * 6;
    duData[pp + 0] = dyData[pp + 5];
    duData[pp + 1] = -dxData[pp + 5];
    duData[pp + 2] = -dyData[pp + 3] + dxData[pp + 4];
    duData[pp + 3] = -dyData[pp + 2];
    duData[pp + 4] = dxData[pp + 2];
    duData[pp + 5] = dyData[pp + 0] - dxData[pp + 1];
  }
}

/// nonlinear 2D HE propagation
void nonlinear2DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c) {

  sunrealtype *dxData = data->buffData[1 - 1];
  sunrealtype *dyData = data->buffData[2 - 1];

  data->exchangeGhostCells(1);
  data->rotateIntoEigen(1);
  data->derive(1);
  data->derotate(1, dxData);
  data->exchangeGhostCells(2);
  data->rotateIntoEigen(2);
  data->derive(2);
  data->derotate(2, dyData);

  sunrealtype f = NAN, g = NAN;
  sunrealtype lf = NAN, lff = NAN, lfg = NAN, lg = NAN, lgg = NAN;
  array<sunrealtype, 36> JMM;
  array<sunrealtype, 6> Quad;
  array<sunrealtype, 6> h;
  sunrealtype pseudoDenom = NAN;
  sunrealtype *udata = nullptr, *dudata = nullptr;
  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);
  int totalNP = data->discreteSize();
  for (int pp = 0; pp < totalNP * 6; pp += 6) {
    // 1
    f = 0.5 * ((Quad[0] = udata[pp] * udata[pp]) +
               (Quad[1] = udata[pp + 1] * udata[pp + 1]) +
               (Quad[2] = udata[pp + 2] * udata[pp + 2]) -
               (Quad[3] = udata[pp + 3] * udata[pp + 3]) -
               (Quad[4] = udata[pp + 4] * udata[pp + 4]) -
               (Quad[5] = udata[pp + 5] * udata[pp + 5]));
    g = udata[pp] * udata[pp + 3] + udata[pp + 1] * udata[pp + 4] +
        udata[pp + 2] * udata[pp + 5];
    // 2
    switch (*c) {
    case 0:
      lf = 0;
      lff = 0;
      lfg = 0;
      lg = 0;
      lgg = 0;
      break;
    case 2:
      lf = 0.000354046449700427580438254 * f * f +
           0.000191775160254398272737387 * g * g;
      lff = 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = 0.0003835503205087965454747749 * f * g;
      lgg = 0.0003835503205087965454747749 * f;
      break;
    case 1:
      lf = 0.000206527095658582755255648 * f;
      lff = 0.000206527095658582755255648;
      lfg = 0;
      lg = 0.0003614224174025198216973841 * g;
      lgg = 0.0003614224174025198216973841;
      break;
    case 3:
      lf = (0.000206527095658582755255648 + 0.000354046449700427580438254 * f) *
               f +
           0.000191775160254398272737387 * g * g;
      lff = 0.000206527095658582755255648 + 0.000708092899400855160876508 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = (0.000361422417402519821697384 + 0.000383550320508796545474775 * f) *
           g;
      lgg = 0.000361422417402519821697384 + 0.000383550320508796545474775 * f;
      break;
    default:
      errorKill(
          "You need to specify a correct order in the weak-field expansion.");
    }
    // 3
    JMM[0] = lf + lff * Quad[0] +
             udata[3 + pp] * (2 * lfg * udata[pp] + lgg * udata[3 + pp]);
    JMM[6] =
        lff * udata[pp] * udata[1 + pp] + lfg * udata[1 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[4 + pp] + lgg * udata[3 + pp] * udata[4 + pp];
    JMM[7] = lf + lff * Quad[1] +
             udata[4 + pp] * (2 * lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[12] =
        lff * udata[pp] * udata[2 + pp] + lfg * udata[2 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[5 + pp] + lgg * udata[3 + pp] * udata[5 + pp];
    JMM[13] = lff * udata[1 + pp] * udata[2 + pp] +
              lfg * udata[2 + pp] * udata[4 + pp] +
              lfg * udata[1 + pp] * udata[5 + pp] +
              lgg * udata[4 + pp] * udata[5 + pp];
    JMM[14] = lf + lff * Quad[2] +
              udata[5 + pp] * (2 * lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[18] = lg + lfg * (Quad[0] - Quad[3 + 0]) +
              (-lff + lgg) * udata[pp] * udata[3 + pp];
    JMM[19] = -(udata[3 + pp] * (lff * udata[1 + pp] + lfg * udata[4 + pp])) +
              udata[pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[20] = -(udata[3 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[21] = -lf + lgg * Quad[0] +
              udata[3 + pp] * (-2 * lfg * udata[pp] + lff * udata[3 + pp]);
    JMM[24] = udata[1 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[4 + pp];
    JMM[25] = lg + lfg * (Quad[1] - Quad[4 + 0]) +
              (-lff + lgg) * udata[1 + pp] * udata[4 + pp];
    JMM[26] = -(udata[4 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[1 + pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[27] = lgg * udata[pp] * udata[1 + pp] +
              lff * udata[3 + pp] * udata[4 + pp] -
              lfg * (udata[1 + pp] * udata[3 + pp] + udata[pp] * udata[4 + pp]);
    JMM[28] = -lf + lgg * Quad[1] +
              udata[4 + pp] * (-2 * lfg * udata[1 + pp] + lff * udata[4 + pp]);
    JMM[30] = udata[2 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[5 + pp];
    JMM[31] = udata[2 + pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]) -
              (lff * udata[1 + pp] + lfg * udata[4 + pp]) * udata[5 + pp];
    JMM[32] = lg + lfg * (Quad[2] - Quad[5 + 0]) +
              (-lff + lgg) * udata[2 + pp] * udata[5 + pp];
    JMM[33] = lgg * udata[pp] * udata[2 + pp] +
              lff * udata[3 + pp] * udata[5 + pp] -
              lfg * (udata[2 + pp] * udata[3 + pp] + udata[pp] * udata[5 + pp]);
    JMM[34] =
        lgg * udata[1 + pp] * udata[2 + pp] +
        lff * udata[4 + pp] * udata[5 + pp] -
        lfg * (udata[2 + pp] * udata[4 + pp] + udata[1 + pp] * udata[5 + pp]);
    JMM[35] = -lf + lgg * Quad[2] +
              udata[5 + pp] * (-2 * lfg * udata[2 + pp] + lff * udata[5 + pp]);
    // 4
    for (int i = 0; i < 6; i++) {
      for (int j = i + 1; j < 6; j++) {
        JMM[i * 6 + j] = JMM[j * 6 + i];
      }
    }
    h[0] = 0;
    h[1] = dxData[pp] * JMM[30] + dxData[1 + pp] * JMM[31] +
           dxData[2 + pp] * JMM[32] + dxData[3 + pp] * JMM[33] +
           dxData[4 + pp] * JMM[34] + dxData[5 + pp] * (-1 + JMM[35]);
    h[2] = -(dxData[pp] * JMM[24]) - dxData[1 + pp] * JMM[25] -
           dxData[2 + pp] * JMM[26] - dxData[3 + pp] * JMM[27] +
           dxData[4 + pp] * (1 - JMM[28]) - dxData[5 + pp] * JMM[29];
    h[3] = 0;
    h[4] = dxData[2 + pp];
    h[5] = -dxData[1 + pp];
    h[0] += -(dyData[pp] * JMM[30]) - dyData[1 + pp] * JMM[31] -
            dyData[2 + pp] * JMM[32] - dyData[3 + pp] * JMM[33] -
            dyData[4 + pp] * JMM[34] + dyData[5 + pp] * (1 - JMM[35]);
    h[1] += 0;
    h[2] += dyData[pp] * JMM[18] + dyData[1 + pp] * JMM[19] +
            dyData[2 + pp] * JMM[20] + dyData[3 + pp] * (-1 + JMM[21]) +
            dyData[4 + pp] * JMM[22] + dyData[5 + pp] * JMM[23];
    h[3] += -dyData[2 + pp];
    h[4] += 0;
    h[5] += dyData[pp];
    h[0] -= h[3] * JMM[3] + h[4] * JMM[4] + h[5] * JMM[5];
    h[1] -= h[3] * JMM[9] + h[4] * JMM[10] + h[5] * JMM[11];
    h[2] -= h[3] * JMM[15] + h[4] * JMM[16] + h[5] * JMM[17];
    dudata[pp + 0] =
        h[2] * (-(JMM[2] * (1 + JMM[7])) + JMM[1] * JMM[8]) +
        h[1] * (JMM[2] * JMM[13] - JMM[1] * (1 + JMM[14])) +
        h[0] * (1 - JMM[8] * JMM[13] + JMM[14] + JMM[7] * (1 + JMM[14]));
    dudata[pp + 1] =
        h[2] * (JMM[2] * JMM[6] - (1 + JMM[0]) * JMM[8]) +
        h[1] * (1 - JMM[2] * JMM[12] + JMM[14] + JMM[0] * (1 + JMM[14])) +
        h[0] * (JMM[8] * JMM[12] - JMM[6] * (1 + JMM[14]));
    dudata[pp + 2] =
        h[2] * (1 - JMM[1] * JMM[6] + JMM[7] + JMM[0] * (1 + JMM[7])) +
        h[1] * (JMM[1] * JMM[12] - (1 + JMM[0]) * JMM[13]) +
        h[0] * (-((1 + JMM[7]) * JMM[12]) + JMM[6] * JMM[13]);
    pseudoDenom =
        -((1 + JMM[7]) * (-1 + JMM[2] * JMM[12])) +
        (JMM[2] * JMM[6] - JMM[8]) * JMM[13] + JMM[14] + JMM[7] * JMM[14] +
        JMM[0] * (1 + JMM[7] - JMM[8] * JMM[13] + (1 + JMM[7]) * JMM[14]) -
        JMM[1] * (-(JMM[8] * JMM[12]) + JMM[6] * (1 + JMM[14]));
    dudata[pp + 0] /= pseudoDenom;
    dudata[pp + 1] /= pseudoDenom;
    dudata[pp + 2] /= pseudoDenom;
    dudata[pp + 3] = h[3];
    dudata[pp + 4] = h[4];
    dudata[pp + 5] = h[5];
  }
  return;
}

/// only under-the-hood-callable Maxwell propagation in 3D
// unused parameters 2-4 for compliance with CVRhsFn
// same as the respective nonlinear function without nonlinear terms
void linear3DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c) {

  sunrealtype *duData = data->duData;
  sunrealtype *dxData = data->buffData[1 - 1];
  sunrealtype *dyData = data->buffData[2 - 1];
  sunrealtype *dzData = data->buffData[3 - 1];

  data->exchangeGhostCells(1);
  data->rotateIntoEigen(1);
  data->derive(1);
  data->derotate(1, dxData);
  data->exchangeGhostCells(2);
  data->rotateIntoEigen(2);
  data->derive(2);
  data->derotate(2, dyData);
  data->exchangeGhostCells(3);
  data->rotateIntoEigen(3);
  data->derive(3);
  data->derotate(3, dzData);

  int totalNP = data->discreteSize();
  int pp = 0;
  for (int i = 0; i < totalNP; i++) {
    pp = i * 6;
    duData[pp + 0] = dyData[pp + 5] - dzData[pp + 4];
    duData[pp + 1] = dzData[pp + 3] - dxData[pp + 5];
    duData[pp + 2] = dxData[pp + 4] - dyData[pp + 3];
    duData[pp + 3] = -dyData[pp + 2] + dzData[pp + 1];
    duData[pp + 4] = -dzData[pp + 0] + dxData[pp + 2];
    duData[pp + 5] = -dxData[pp + 1] + dyData[pp + 0];
  }
}

/// nonlinear 3D HE propagation
void nonlinear3DProp(LatticePatch *data, N_Vector u, N_Vector udot, int *c) {

  sunrealtype *dxData = data->buffData[1 - 1];
  sunrealtype *dyData = data->buffData[2 - 1];
  sunrealtype *dzData = data->buffData[3 - 1];

  data->exchangeGhostCells(1);
  data->rotateIntoEigen(1);
  data->derive(1);
  data->derotate(1,dxData);
  data->exchangeGhostCells(2);
  data->rotateIntoEigen(2);
  data->derive(2);
  data->derotate(2,dyData);
  data->exchangeGhostCells(3);
  data->rotateIntoEigen(3);
  data->derive(3);
  data->derotate(3,dzData);
  
  sunrealtype f = NAN, g = NAN;
  sunrealtype lf = NAN, lff = NAN, lfg = NAN, lg = NAN, lgg = NAN;
  array<sunrealtype, 36> JMM;
  array<sunrealtype, 6> Quad;
  array<sunrealtype, 6> h;
  sunrealtype pseudoDenom = NAN;
  sunrealtype *udata = nullptr, *dudata = nullptr;
  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);
  int totalNP = data->discreteSize();
  for (int pp = 0; pp < totalNP * 6; pp += 6) {
    // 1
    f = 0.5 * ((Quad[0] = udata[pp] * udata[pp]) +
               (Quad[1] = udata[pp + 1] * udata[pp + 1]) +
               (Quad[2] = udata[pp + 2] * udata[pp + 2]) -
               (Quad[3] = udata[pp + 3] * udata[pp + 3]) -
               (Quad[4] = udata[pp + 4] * udata[pp + 4]) -
               (Quad[5] = udata[pp + 5] * udata[pp + 5]));
    g = udata[pp] * udata[pp + 3] + udata[pp + 1] * udata[pp + 4] +
        udata[pp + 2] * udata[pp + 5];
    // 2
    switch (*c) {
    case 0:
      lf = 0;
      lff = 0;
      lfg = 0;
      lg = 0;
      lgg = 0;
      break;
    case 2:
      lf = 0.000354046449700427580438254 * f * f +
           0.000191775160254398272737387 * g * g;
      lff = 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = 0.0003835503205087965454747749 * f * g;
      lgg = 0.0003835503205087965454747749 * f;
      break;
    case 1:
      lf = 0.000206527095658582755255648 * f;
      lff = 0.000206527095658582755255648;
      lfg = 0;
      lg = 0.0003614224174025198216973841 * g;
      lgg = 0.0003614224174025198216973841;
      break;
    case 3:
      lf = (0.000206527095658582755255648 + 0.000354046449700427580438254 * f) *
               f +
           0.000191775160254398272737387 * g * g;
      lff = 0.000206527095658582755255648 + 0.000708092899400855160876508 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = (0.000361422417402519821697384 + 0.000383550320508796545474775 * f) *
           g;
      lgg = 0.000361422417402519821697384 + 0.000383550320508796545474775 * f;
      break;
    default:
      errorKill(
          "You need to specify a correct order in the weak-field expansion.");
    }
    // 3
    JMM[0] = lf + lff * Quad[0] +
             udata[3 + pp] * (2 * lfg * udata[pp] + lgg * udata[3 + pp]);
    JMM[6] =
        lff * udata[pp] * udata[1 + pp] + lfg * udata[1 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[4 + pp] + lgg * udata[3 + pp] * udata[4 + pp];
    JMM[7] = lf + lff * Quad[1] +
             udata[4 + pp] * (2 * lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[12] =
        lff * udata[pp] * udata[2 + pp] + lfg * udata[2 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[5 + pp] + lgg * udata[3 + pp] * udata[5 + pp];
    JMM[13] = lff * udata[1 + pp] * udata[2 + pp] +
              lfg * udata[2 + pp] * udata[4 + pp] +
              lfg * udata[1 + pp] * udata[5 + pp] +
              lgg * udata[4 + pp] * udata[5 + pp];
    JMM[14] = lf + lff * Quad[2] +
              udata[5 + pp] * (2 * lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[18] = lg + lfg * (Quad[0] - Quad[3 + 0]) +
              (-lff + lgg) * udata[pp] * udata[3 + pp];
    JMM[19] = -(udata[3 + pp] * (lff * udata[1 + pp] + lfg * udata[4 + pp])) +
              udata[pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[20] = -(udata[3 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[21] = -lf + lgg * Quad[0] +
              udata[3 + pp] * (-2 * lfg * udata[pp] + lff * udata[3 + pp]);
    JMM[24] = udata[1 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[4 + pp];
    JMM[25] = lg + lfg * (Quad[1] - Quad[4 + 0]) +
              (-lff + lgg) * udata[1 + pp] * udata[4 + pp];
    JMM[26] = -(udata[4 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[1 + pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[27] = lgg * udata[pp] * udata[1 + pp] +
              lff * udata[3 + pp] * udata[4 + pp] -
              lfg * (udata[1 + pp] * udata[3 + pp] + udata[pp] * udata[4 + pp]);
    JMM[28] = -lf + lgg * Quad[1] +
              udata[4 + pp] * (-2 * lfg * udata[1 + pp] + lff * udata[4 + pp]);
    JMM[30] = udata[2 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[5 + pp];
    JMM[31] = udata[2 + pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]) -
              (lff * udata[1 + pp] + lfg * udata[4 + pp]) * udata[5 + pp];
    JMM[32] = lg + lfg * (Quad[2] - Quad[5 + 0]) +
              (-lff + lgg) * udata[2 + pp] * udata[5 + pp];
    JMM[33] = lgg * udata[pp] * udata[2 + pp] +
              lff * udata[3 + pp] * udata[5 + pp] -
              lfg * (udata[2 + pp] * udata[3 + pp] + udata[pp] * udata[5 + pp]);
    JMM[34] =
        lgg * udata[1 + pp] * udata[2 + pp] +
        lff * udata[4 + pp] * udata[5 + pp] -
        lfg * (udata[2 + pp] * udata[4 + pp] + udata[1 + pp] * udata[5 + pp]);
    JMM[35] = -lf + lgg * Quad[2] +
              udata[5 + pp] * (-2 * lfg * udata[2 + pp] + lff * udata[5 + pp]);
    // 4
    for (int i = 0; i < 6; i++) {
      for (int j = i + 1; j < 6; j++) {
        JMM[i * 6 + j] = JMM[j * 6 + i];
      }
    }
    h[0] = 0;
    h[1] = dxData[pp] * JMM[30] + dxData[1 + pp] * JMM[31] +
           dxData[2 + pp] * JMM[32] + dxData[3 + pp] * JMM[33] +
           dxData[4 + pp] * JMM[34] + dxData[5 + pp] * (-1 + JMM[35]);
    h[2] = -(dxData[pp] * JMM[24]) - dxData[1 + pp] * JMM[25] -
           dxData[2 + pp] * JMM[26] - dxData[3 + pp] * JMM[27] +
           dxData[4 + pp] * (1 - JMM[28]) - dxData[5 + pp] * JMM[29];
    h[3] = 0;
    h[4] = dxData[2 + pp];
    h[5] = -dxData[1 + pp];
    h[0] += -(dyData[pp] * JMM[30]) - dyData[1 + pp] * JMM[31] -
            dyData[2 + pp] * JMM[32] - dyData[3 + pp] * JMM[33] -
            dyData[4 + pp] * JMM[34] + dyData[5 + pp] * (1 - JMM[35]);
    h[1] += 0;
    h[2] += dyData[pp] * JMM[18] + dyData[1 + pp] * JMM[19] +
            dyData[2 + pp] * JMM[20] + dyData[3 + pp] * (-1 + JMM[21]) +
            dyData[4 + pp] * JMM[22] + dyData[5 + pp] * JMM[23];
    h[3] += -dyData[2 + pp];
    h[4] += 0;
    h[5] += dyData[pp];
    h[0] += dzData[pp] * JMM[24] + dzData[1 + pp] * JMM[25] +
            dzData[2 + pp] * JMM[26] + dzData[3 + pp] * JMM[27] +
            dzData[4 + pp] * (-1 + JMM[28]) + dzData[5 + pp] * JMM[29];
    h[1] += -(dzData[pp] * JMM[18]) - dzData[1 + pp] * JMM[19] -
            dzData[2 + pp] * JMM[20] + dzData[3 + pp] * (1 - JMM[21]) -
            dzData[4 + pp] * JMM[22] - dzData[5 + pp] * JMM[23];
    h[2] += 0;
    h[3] += dzData[1 + pp];
    h[4] += -dzData[pp];
    h[5] += 0;
    h[0] -= h[3] * JMM[3] + h[4] * JMM[4] + h[5] * JMM[5];
    h[1] -= h[3] * JMM[9] + h[4] * JMM[10] + h[5] * JMM[11];
    h[2] -= h[3] * JMM[15] + h[4] * JMM[16] + h[5] * JMM[17];
    dudata[pp + 0] =
        h[2] * (-(JMM[2] * (1 + JMM[7])) + JMM[1] * JMM[8]) +
        h[1] * (JMM[2] * JMM[13] - JMM[1] * (1 + JMM[14])) +
        h[0] * (1 - JMM[8] * JMM[13] + JMM[14] + JMM[7] * (1 + JMM[14]));
    dudata[pp + 1] =
        h[2] * (JMM[2] * JMM[6] - (1 + JMM[0]) * JMM[8]) +
        h[1] * (1 - JMM[2] * JMM[12] + JMM[14] + JMM[0] * (1 + JMM[14])) +
        h[0] * (JMM[8] * JMM[12] - JMM[6] * (1 + JMM[14]));
    dudata[pp + 2] =
        h[2] * (1 - JMM[1] * JMM[6] + JMM[7] + JMM[0] * (1 + JMM[7])) +
        h[1] * (JMM[1] * JMM[12] - (1 + JMM[0]) * JMM[13]) +
        h[0] * (-((1 + JMM[7]) * JMM[12]) + JMM[6] * JMM[13]);
    pseudoDenom =
        -((1 + JMM[7]) * (-1 + JMM[2] * JMM[12])) +
        (JMM[2] * JMM[6] - JMM[8]) * JMM[13] + JMM[14] + JMM[7] * JMM[14] +
        JMM[0] * (1 + JMM[7] - JMM[8] * JMM[13] + (1 + JMM[7]) * JMM[14]) -
        JMM[1] * (-(JMM[8] * JMM[12]) + JMM[6] * (1 + JMM[14]));
    dudata[pp + 0] /= pseudoDenom;
    dudata[pp + 1] /= pseudoDenom;
    dudata[pp + 2] /= pseudoDenom;
    dudata[pp + 3] = h[3];
    dudata[pp + 4] = h[4];
    dudata[pp + 5] = h[5];
  }
  return;
}
