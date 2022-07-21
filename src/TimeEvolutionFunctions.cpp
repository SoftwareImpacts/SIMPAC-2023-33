//////////////////////////////////////////////////////////
/// @file TimeEvolutionFunctions.cpp
/// @brief Implementation of functions to propagate
/// data vectors in time according to Maxwell's equations,
/// and various orders in the HE weak-field expansion
//////////////////////////////////////////////////////////

#include "TimeEvolutionFunctions.h"

/// CVode right-hand-side function (CVRhsFn)
int TimeEvolution::f(sunrealtype t, N_Vector u, N_Vector udot, void *data_loc) {

  // Set recover pointer to provided lattice patch where the field data resides
  LatticePatch *data = static_cast<LatticePatch *>(data_loc);

  // update circle
  // Access data with pointers
  sunrealtype *udata = NV_DATA_P(u),
              *dudata = NV_DATA_P(udot);

  // Store original data location of the patch
  sunrealtype *originaluData = data->uData,
              *originalduData = data->duData;

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

/// only under-the-hood-callable Maxwell propagation in 1D;
// unused parameters 2-4 for compliance with CVRhsFn;
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

  const sunindextype totalNP = data->discreteSize();
  sunindextype pp = 0;
  for (sunindextype i = 0; i < totalNP; i++) {
    pp = i * 6;
    /*
     simple vacuum Maxwell equations for spatial deriative only in x-direction;
     temporal derivative is approximated by spatial derivative according to the
     numerical scheme without polarization or magnetization terms
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

  // pointer to spatial derivative data sufficient; temporal derivative data
  // provided with udot
  sunrealtype *dxData = data->buffData[1 - 1];

  // same sequence as in the linear case
  data->exchangeGhostCells(1);
  data->rotateIntoEigen(1);
  data->derive(1);
  data->derotate(1, dxData);

  /*
  F and G are nonzero in the nonlinear case,
  polarization and magnetization derivatives
  w.r.t. E- and B-field go into the e.o.m.
  */
  static sunrealtype f, g; // em field invariants F, G
  // derivatives of HE Lagrangian w.r.t. field invariants
  static sunrealtype lf, lff, lfg, lg, lgg;
  // matrix to hold derivatives of polarization and magnetization
  static std::array<sunrealtype, 21> JMM;
  // array to hold E^2 and B^2 components
  static std::array<sunrealtype, 6> Quad;
  // array to hold intermediate temp. derivatives of E and B
  static std::array<sunrealtype, 6> h;
  // determinant needed for explicit matrix inversion
  static sunrealtype detC = nan("0x12345");
  // pointers to field values and their temp. derivatives
  sunrealtype *udata = NV_DATA_P(u),
              *dudata = NV_DATA_P(udot);

  // number of points in the patch
  const sunindextype totalNP = data->discreteSize();
  for (sunindextype pp = 0; pp < totalNP * 6;
       pp += 6) { // loop over all 6dim points in the patch
    // em field Lorentz invariants F and G
    f = 0.5 * ((Quad[0] = udata[pp] * udata[pp]) +
               (Quad[1] = udata[pp + 1] * udata[pp + 1]) +
               (Quad[2] = udata[pp + 2] * udata[pp + 2]) -
               (Quad[3] = udata[pp + 3] * udata[pp + 3]) -
               (Quad[4] = udata[pp + 4] * udata[pp + 4]) -
               (Quad[5] = udata[pp + 5] * udata[pp + 5]));
    g = udata[pp] * udata[pp + 3] + udata[pp + 1] * udata[pp + 4] +
        udata[pp + 2] * udata[pp + 5];
    // process/expansion order and corresponding derivative values of L
    // w.r.t. F, G
    switch (*c) {
    case 0:  // linear Maxwell vacuum
      lf = 0;
      lff = 0;
      lfg = 0;
      lg = 0;
      lgg = 0;
      break;
    case 1:  // only 4-photon processes
      lf = 0.000206527095658582755255648 * f;
      lff = 0.000206527095658582755255648;
      lfg = 0;
      lg = 0.0003614224174025198216973841 * g;
      lgg = 0.0003614224174025198216973841;
      break;
    case 2:  // only 6-photon processes
      lf = 0.000354046449700427580438254 * f * f +
           0.000191775160254398272737387 * g * g;
      lff = 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = 0.0003835503205087965454747749 * f * g;
      lgg = 0.0003835503205087965454747749 * f;
      break;
    case 3:  // 4- and 6-photon processes
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

    // derivatives of polarization and magnetization w.r.t. E and B
    // Jpx(Ex)
    JMM[0] = lf + lff * Quad[0] +
             udata[3 + pp] * (2 * lfg * udata[pp] + lgg * udata[3 + pp]);
    // Jpx(Ey)
    JMM[1] =
        lff * udata[pp] * udata[1 + pp] + lfg * udata[1 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[4 + pp] + lgg * udata[3 + pp] * udata[4 + pp];
    // Jpy(Ey)
    JMM[2] = lf + lff * Quad[1] +
             udata[4 + pp] * (2 * lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    // Jpx(Ez) = Jpz(Ex)
    JMM[3] =
        lff * udata[pp] * udata[2 + pp] + lfg * udata[2 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[5 + pp] + lgg * udata[3 + pp] * udata[5 + pp];
    // Jpy(Ez) = Jpz(Ey)
    JMM[4] = lff * udata[1 + pp] * udata[2 + pp] +
              lfg * udata[2 + pp] * udata[4 + pp] +
              lfg * udata[1 + pp] * udata[5 + pp] +
              lgg * udata[4 + pp] * udata[5 + pp];
    // Jpz(Ez)
    JMM[5] = lf + lff * Quad[2] +
              udata[5 + pp] * (2 * lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    // Jpx(Bx) = Jmx(Ex)
    JMM[6] = lg + lfg * (Quad[0] - Quad[3 + 0]) +
              (-lff + lgg) * udata[pp] * udata[3 + pp];
    // Jpy(Bx) = Jmx(Ey)
    JMM[7] = -(udata[3 + pp] * (lff * udata[1 + pp] + lfg * udata[4 + pp])) +
              udata[pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    // Jpz(Bx) = Jmx(Ez)
    JMM[8] = -(udata[3 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    // Jmx(Bx)
    JMM[9] = -lf + lgg * Quad[0] +
              udata[3 + pp] * (-2 * lfg * udata[pp] + lff * udata[3 + pp]);
    // Jpx(By) = Jmy(Ex)
    JMM[10] = udata[1 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[4 + pp];
    // Jpy(By) = Jmy(Ey)
    JMM[11] = lg + lfg * (Quad[1] - Quad[4 + 0]) +
              (-lff + lgg) * udata[1 + pp] * udata[4 + pp];
    // Jpz(By) = Jmy(Ez)
    JMM[12] = -(udata[4 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[1 + pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    // Jmx(By) = Jmy(Bx)
    JMM[13] = lgg * udata[pp] * udata[1 + pp] +
              lff * udata[3 + pp] * udata[4 + pp] -
              lfg * (udata[1 + pp] * udata[3 + pp] + udata[pp] * udata[4 + pp]);
    // Jmy(By)
    JMM[14] = -lf + lgg * Quad[1] +
              udata[4 + pp] * (-2 * lfg * udata[1 + pp] + lff * udata[4 + pp]);
    // Jmz(Ex) = Jpx(Bz)
    JMM[15] = udata[2 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[5 + pp];
    // Jmz(Ey) = Jpy(Bz)
    JMM[16] = udata[2 + pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]) -
              (lff * udata[1 + pp] + lfg * udata[4 + pp]) * udata[5 + pp];
    // Jpz(Bz) = Jmz(Ez)
    JMM[17] = lg + lfg * (Quad[2] - Quad[5 + 0]) +
              (-lff + lgg) * udata[2 + pp] * udata[5 + pp];
    // Jmz(Bx) = Jmx(Bz)
    JMM[18] = lgg * udata[pp] * udata[2 + pp] +
              lff * udata[3 + pp] * udata[5 + pp] -
              lfg * (udata[2 + pp] * udata[3 + pp] + udata[pp] * udata[5 + pp]);
    // Jmy(Bz) = Jmz(By)
    JMM[19] =
        lgg * udata[1 + pp] * udata[2 + pp] +
        lff * udata[4 + pp] * udata[5 + pp] -
        lfg * (udata[2 + pp] * udata[4 + pp] + udata[1 + pp] * udata[5 + pp]);
    // Jmz(Bz)
    JMM[20] = -lf + lgg * Quad[2] +
              udata[5 + pp] * (-2 * lfg * udata[2 + pp] + lff * udata[5 + pp]);

    // apply Z
    // top bock: -QJm(E)*E, Q-QJm(B)*B
    h[0] = 0;
    h[1] = dxData[pp] * JMM[15] + dxData[1 + pp] * JMM[16] +
           dxData[2 + pp] * JMM[17] + dxData[3 + pp] * JMM[18] +
           dxData[4 + pp] * JMM[19] + dxData[5 + pp] * (-1 + JMM[20]);
    h[2] = -(dxData[pp] * JMM[10]) - dxData[1 + pp] * JMM[11] -
           dxData[2 + pp] * JMM[12] - dxData[3 + pp] * JMM[13] +
           dxData[4 + pp] * (1 - JMM[14]) - dxData[5 + pp] * JMM[19];
    // bottom blocks: -Q*E
    h[3] = 0;
    h[4] = dxData[2 + pp];
    h[5] = -dxData[1 + pp];
    // (1+A)^-1 applies only to E components
    // -Jp(B)*B
    h[0] -= h[3] * JMM[6] + h[4] * JMM[10] + h[5] * JMM[15];
    h[1] -= h[3] * JMM[7] + h[4] * JMM[11] + h[5] * JMM[16];
    h[2] -= h[3] * JMM[8] + h[4] * JMM[12] + h[5] * JMM[17];
    // apply C^-1 explicitly, with C=1+Jp(E)
    dudata[pp + 0] =
        h[2] * (-(JMM[3] * (1 + JMM[2])) + JMM[1] * JMM[4]) +
        h[1] * (JMM[3] * JMM[4] - JMM[1] * (1 + JMM[5])) +
        h[0] * (1 - JMM[4] * JMM[4] + JMM[5] + JMM[2] * (1 + JMM[5]));
    dudata[pp + 1] =
        h[2] * (JMM[3] * JMM[1] - (1 + JMM[0]) * JMM[4]) +
        h[1] * (1 - JMM[3] * JMM[3] + JMM[5] + JMM[0] * (1 + JMM[5])) +
        h[0] * (JMM[4] * JMM[3] - JMM[1] * (1 + JMM[5]));
    dudata[pp + 2] =
        h[2] * (1 - JMM[1] * JMM[1] + JMM[2] + JMM[0] * (1 + JMM[2])) +
        h[1] * (JMM[1] * JMM[3] - (1 + JMM[0]) * JMM[4]) +
        h[0] * (-((1 + JMM[2]) * JMM[3]) + JMM[1] * JMM[4]);
    detC =  // determinant of C
        -((1 + JMM[2]) * (-1 + JMM[3] * JMM[3])) +
        (JMM[3] * JMM[1] - JMM[4]) * JMM[4] + JMM[5] + JMM[2] * JMM[5] +
        JMM[0] * (1 + JMM[2] - JMM[4] * JMM[4] + (1 + JMM[2]) * JMM[5]) -
        JMM[1] * (-(JMM[4] * JMM[3]) + JMM[1] * (1 + JMM[5]));
    dudata[pp + 0] /= detC;
    dudata[pp + 1] /= detC;
    dudata[pp + 2] /= detC;
    dudata[pp + 3] = h[3];
    dudata[pp + 4] = h[4];
    dudata[pp + 5] = h[5];
  }
  return;
}

/// only under-the-hood-callable Maxwell propagation in 2D;
// unused parameters 2-4 for compliance with CVRhsFn;
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

  const sunindextype totalNP = data->discreteSize();
  sunindextype pp = 0;
  for (sunindextype i = 0; i < totalNP; i++) {
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

  static sunrealtype f, g;
  static sunrealtype lf, lff, lfg, lg, lgg;
  static std::array<sunrealtype, 21> JMM;
  static std::array<sunrealtype, 6> Quad;
  static std::array<sunrealtype, 6> h;
  static sunrealtype detC;
  sunrealtype *udata = NV_DATA_P(u),
              *dudata = NV_DATA_P(udot);

  const sunindextype totalNP = data->discreteSize();
  for (sunindextype pp = 0; pp < totalNP * 6; pp += 6) {
    f = 0.5 * ((Quad[0] = udata[pp] * udata[pp]) +
               (Quad[1] = udata[pp + 1] * udata[pp + 1]) +
               (Quad[2] = udata[pp + 2] * udata[pp + 2]) -
               (Quad[3] = udata[pp + 3] * udata[pp + 3]) -
               (Quad[4] = udata[pp + 4] * udata[pp + 4]) -
               (Quad[5] = udata[pp + 5] * udata[pp + 5]));
    g = udata[pp] * udata[pp + 3] + udata[pp + 1] * udata[pp + 4] +
        udata[pp + 2] * udata[pp + 5];
    switch (*c) {
    case 0:
      lf = 0;
      lff = 0;
      lfg = 0;
      lg = 0;
      lgg = 0;
      break;
    case 1:
      lf = 0.000206527095658582755255648 * f;
      lff = 0.000206527095658582755255648;
      lfg = 0;
      lg = 0.0003614224174025198216973841 * g;
      lgg = 0.0003614224174025198216973841;
      break;
    case 2:
      lf = 0.000354046449700427580438254 * f * f +
           0.000191775160254398272737387 * g * g;
      lff = 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = 0.0003835503205087965454747749 * f * g;
      lgg = 0.0003835503205087965454747749 * f;
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

    JMM[0] = lf + lff * Quad[0] +
             udata[3 + pp] * (2 * lfg * udata[pp] + lgg * udata[3 + pp]);
    JMM[1] =
        lff * udata[pp] * udata[1 + pp] + lfg * udata[1 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[4 + pp] + lgg * udata[3 + pp] * udata[4 + pp];
    JMM[2] = lf + lff * Quad[1] +
             udata[4 + pp] * (2 * lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[3] =
        lff * udata[pp] * udata[2 + pp] + lfg * udata[2 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[5 + pp] + lgg * udata[3 + pp] * udata[5 + pp];
    JMM[4] = lff * udata[1 + pp] * udata[2 + pp] +
              lfg * udata[2 + pp] * udata[4 + pp] +
              lfg * udata[1 + pp] * udata[5 + pp] +
              lgg * udata[4 + pp] * udata[5 + pp];
    JMM[5] = lf + lff * Quad[2] +
              udata[5 + pp] * (2 * lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[6] = lg + lfg * (Quad[0] - Quad[3 + 0]) +
              (-lff + lgg) * udata[pp] * udata[3 + pp];
    JMM[7] = -(udata[3 + pp] * (lff * udata[1 + pp] + lfg * udata[4 + pp])) +
              udata[pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[8] = -(udata[3 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[9] = -lf + lgg * Quad[0] +
              udata[3 + pp] * (-2 * lfg * udata[pp] + lff * udata[3 + pp]);
    JMM[10] = udata[1 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[4 + pp];
    JMM[11] = lg + lfg * (Quad[1] - Quad[4 + 0]) +
              (-lff + lgg) * udata[1 + pp] * udata[4 + pp];
    JMM[12] = -(udata[4 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[1 + pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[13] = lgg * udata[pp] * udata[1 + pp] +
              lff * udata[3 + pp] * udata[4 + pp] -
              lfg * (udata[1 + pp] * udata[3 + pp] + udata[pp] * udata[4 + pp]);
    JMM[14] = -lf + lgg * Quad[1] +
              udata[4 + pp] * (-2 * lfg * udata[1 + pp] + lff * udata[4 + pp]);
    JMM[15] = udata[2 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[5 + pp];
    JMM[16] = udata[2 + pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]) -
              (lff * udata[1 + pp] + lfg * udata[4 + pp]) * udata[5 + pp];
    JMM[17] = lg + lfg * (Quad[2] - Quad[5 + 0]) +
              (-lff + lgg) * udata[2 + pp] * udata[5 + pp];
    JMM[18] = lgg * udata[pp] * udata[2 + pp] +
              lff * udata[3 + pp] * udata[5 + pp] -
              lfg * (udata[2 + pp] * udata[3 + pp] + udata[pp] * udata[5 + pp]);
    JMM[19] =
        lgg * udata[1 + pp] * udata[2 + pp] +
        lff * udata[4 + pp] * udata[5 + pp] -
        lfg * (udata[2 + pp] * udata[4 + pp] + udata[1 + pp] * udata[5 + pp]);
    JMM[20] = -lf + lgg * Quad[2] +
              udata[5 + pp] * (-2 * lfg * udata[2 + pp] + lff * udata[5 + pp]);

    h[0] = 0;
    h[1] = dxData[pp] * JMM[15] + dxData[1 + pp] * JMM[16] +
           dxData[2 + pp] * JMM[17] + dxData[3 + pp] * JMM[18] +
           dxData[4 + pp] * JMM[19] + dxData[5 + pp] * (-1 + JMM[20]);
    h[2] = -(dxData[pp] * JMM[10]) - dxData[1 + pp] * JMM[11] -
           dxData[2 + pp] * JMM[12] - dxData[3 + pp] * JMM[13] +
           dxData[4 + pp] * (1 - JMM[14]) - dxData[5 + pp] * JMM[19];
    h[3] = 0;
    h[4] = dxData[2 + pp];
    h[5] = -dxData[1 + pp];
    h[0] += -(dyData[pp] * JMM[15]) - dyData[1 + pp] * JMM[16] -
            dyData[2 + pp] * JMM[17] - dyData[3 + pp] * JMM[18] -
            dyData[4 + pp] * JMM[19] + dyData[5 + pp] * (1 - JMM[20]);
    h[1] += 0;
    h[2] += dyData[pp] * JMM[6] + dyData[1 + pp] * JMM[7] +
            dyData[2 + pp] * JMM[8] + dyData[3 + pp] * (-1 + JMM[9]) +
            dyData[4 + pp] * JMM[13] + dyData[5 + pp] * JMM[18];
    h[3] += -dyData[2 + pp];
    h[4] += 0;
    h[5] += dyData[pp];
    h[0] -= h[3] * JMM[6] + h[4] * JMM[10] + h[5] * JMM[15];
    h[1] -= h[3] * JMM[7] + h[4] * JMM[11] + h[5] * JMM[16];
    h[2] -= h[3] * JMM[8] + h[4] * JMM[12] + h[5] * JMM[17];
    dudata[pp + 0] =
        h[2] * (-(JMM[3] * (1 + JMM[2])) + JMM[1] * JMM[4]) +
        h[1] * (JMM[3] * JMM[4] - JMM[1] * (1 + JMM[5])) +
        h[0] * (1 - JMM[4] * JMM[4] + JMM[5] + JMM[2] * (1 + JMM[5]));
    dudata[pp + 1] =
        h[2] * (JMM[3] * JMM[1] - (1 + JMM[0]) * JMM[4]) +
        h[1] * (1 - JMM[3] * JMM[3] + JMM[5] + JMM[0] * (1 + JMM[5])) +
        h[0] * (JMM[4] * JMM[3] - JMM[1] * (1 + JMM[5]));
    dudata[pp + 2] =
        h[2] * (1 - JMM[1] * JMM[1] + JMM[2] + JMM[0] * (1 + JMM[2])) +
        h[1] * (JMM[1] * JMM[3] - (1 + JMM[0]) * JMM[4]) +
        h[0] * (-((1 + JMM[2]) * JMM[3]) + JMM[1] * JMM[4]);
    detC =
        -((1 + JMM[2]) * (-1 + JMM[3] * JMM[3])) +
        (JMM[3] * JMM[1] - JMM[4]) * JMM[4] + JMM[5] + JMM[2] * JMM[5] +
        JMM[0] * (1 + JMM[2] - JMM[4] * JMM[4] + (1 + JMM[2]) * JMM[5]) -
        JMM[1] * (-(JMM[4] * JMM[3]) + JMM[1] * (1 + JMM[5]));
    dudata[pp + 0] /= detC;
    dudata[pp + 1] /= detC;
    dudata[pp + 2] /= detC;
    dudata[pp + 3] = h[3];
    dudata[pp + 4] = h[4];
    dudata[pp + 5] = h[5];
  }
  return;
}

/// only under-the-hood-callable Maxwell propagation in 3D;
// unused parameters 2-4 for compliance with CVRhsFn;
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

  const sunindextype totalNP = data->discreteSize();
  sunindextype pp = 0;
  for (sunindextype i = 0; i < totalNP; i++) {
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

  static sunrealtype f, g;
  static sunrealtype lf, lff, lfg, lg, lgg;
  static std::array<sunrealtype, 21> JMM;
  static std::array<sunrealtype, 6> Quad;
  static std::array<sunrealtype, 6> h;
  static sunrealtype detC = nan("0x12345");
  sunrealtype *udata = NV_DATA_P(u),
              *dudata = NV_DATA_P(udot);

  const sunindextype totalNP = data->discreteSize();
  for (sunindextype pp = 0; pp < totalNP * 6; pp += 6) {
    f = 0.5 * ((Quad[0] = udata[pp] * udata[pp]) +
               (Quad[1] = udata[pp + 1] * udata[pp + 1]) +
               (Quad[2] = udata[pp + 2] * udata[pp + 2]) -
               (Quad[3] = udata[pp + 3] * udata[pp + 3]) -
               (Quad[4] = udata[pp + 4] * udata[pp + 4]) -
               (Quad[5] = udata[pp + 5] * udata[pp + 5]));
    g = udata[pp] * udata[pp + 3] + udata[pp + 1] * udata[pp + 4] +
        udata[pp + 2] * udata[pp + 5];
    switch (*c) {
    case 0:
      lf = 0;
      lff = 0;
      lfg = 0;
      lg = 0;
      lgg = 0;
      break;
    case 1:
      lf = 0.000206527095658582755255648 * f;
      lff = 0.000206527095658582755255648;
      lfg = 0;
      lg = 0.0003614224174025198216973841 * g;
      lgg = 0.0003614224174025198216973841;
      break;
    case 2:
      lf = 0.000354046449700427580438254 * f * f +
           0.000191775160254398272737387 * g * g;
      lff = 0.0007080928994008551608765075 * f;
      lfg = 0.0003835503205087965454747749 * g;
      lg = 0.0003835503205087965454747749 * f * g;
      lgg = 0.0003835503205087965454747749 * f;
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

    JMM[0] = lf + lff * Quad[0] +
             udata[3 + pp] * (2 * lfg * udata[pp] + lgg * udata[3 + pp]);
    JMM[1] =
        lff * udata[pp] * udata[1 + pp] + lfg * udata[1 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[4 + pp] + lgg * udata[3 + pp] * udata[4 + pp];
    JMM[2] = lf + lff * Quad[1] +
             udata[4 + pp] * (2 * lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[3] =
        lff * udata[pp] * udata[2 + pp] + lfg * udata[2 + pp] * udata[3 + pp] +
        lfg * udata[pp] * udata[5 + pp] + lgg * udata[3 + pp] * udata[5 + pp];
    JMM[4] = lff * udata[1 + pp] * udata[2 + pp] +
              lfg * udata[2 + pp] * udata[4 + pp] +
              lfg * udata[1 + pp] * udata[5 + pp] +
              lgg * udata[4 + pp] * udata[5 + pp];
    JMM[5] = lf + lff * Quad[2] +
              udata[5 + pp] * (2 * lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[6] = lg + lfg * (Quad[0] - Quad[3 + 0]) +
              (-lff + lgg) * udata[pp] * udata[3 + pp];
    JMM[7] = -(udata[3 + pp] * (lff * udata[1 + pp] + lfg * udata[4 + pp])) +
              udata[pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]);
    JMM[8] = -(udata[3 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[9] = -lf + lgg * Quad[0] +
              udata[3 + pp] * (-2 * lfg * udata[pp] + lff * udata[3 + pp]);
    JMM[10] = udata[1 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[4 + pp];
    JMM[11] = lg + lfg * (Quad[1] - Quad[4 + 0]) +
              (-lff + lgg) * udata[1 + pp] * udata[4 + pp];
    JMM[12] = -(udata[4 + pp] * (lff * udata[2 + pp] + lfg * udata[5 + pp])) +
              udata[1 + pp] * (lfg * udata[2 + pp] + lgg * udata[5 + pp]);
    JMM[13] = lgg * udata[pp] * udata[1 + pp] +
              lff * udata[3 + pp] * udata[4 + pp] -
              lfg * (udata[1 + pp] * udata[3 + pp] + udata[pp] * udata[4 + pp]);
    JMM[14] = -lf + lgg * Quad[1] +
              udata[4 + pp] * (-2 * lfg * udata[1 + pp] + lff * udata[4 + pp]);
    JMM[15] = udata[2 + pp] * (lfg * udata[pp] + lgg * udata[3 + pp]) -
              (lff * udata[pp] + lfg * udata[3 + pp]) * udata[5 + pp];
    JMM[16] = udata[2 + pp] * (lfg * udata[1 + pp] + lgg * udata[4 + pp]) -
              (lff * udata[1 + pp] + lfg * udata[4 + pp]) * udata[5 + pp];
    JMM[17] = lg + lfg * (Quad[2] - Quad[5 + 0]) +
              (-lff + lgg) * udata[2 + pp] * udata[5 + pp];
    JMM[18] = lgg * udata[pp] * udata[2 + pp] +
              lff * udata[3 + pp] * udata[5 + pp] -
              lfg * (udata[2 + pp] * udata[3 + pp] + udata[pp] * udata[5 + pp]);
    JMM[19] =
        lgg * udata[1 + pp] * udata[2 + pp] +
        lff * udata[4 + pp] * udata[5 + pp] -
        lfg * (udata[2 + pp] * udata[4 + pp] + udata[1 + pp] * udata[5 + pp]);
    JMM[20] = -lf + lgg * Quad[2] +
              udata[5 + pp] * (-2 * lfg * udata[2 + pp] + lff * udata[5 + pp]);

    h[0] = 0;
    h[1] = dxData[pp] * JMM[15] + dxData[1 + pp] * JMM[16] +
           dxData[2 + pp] * JMM[17] + dxData[3 + pp] * JMM[18] +
           dxData[4 + pp] * JMM[19] + dxData[5 + pp] * (-1 + JMM[20]);
    h[2] = -(dxData[pp] * JMM[10]) - dxData[1 + pp] * JMM[11] -
           dxData[2 + pp] * JMM[12] - dxData[3 + pp] * JMM[13] +
           dxData[4 + pp] * (1 - JMM[14]) - dxData[5 + pp] * JMM[19];
    h[3] = 0;
    h[4] = dxData[2 + pp];
    h[5] = -dxData[1 + pp];
    h[0] += -(dyData[pp] * JMM[15]) - dyData[1 + pp] * JMM[16] -
            dyData[2 + pp] * JMM[17] - dyData[3 + pp] * JMM[18] -
            dyData[4 + pp] * JMM[19] + dyData[5 + pp] * (1 - JMM[20]);
    h[1] += 0;
    h[2] += dyData[pp] * JMM[6] + dyData[1 + pp] * JMM[7] +
            dyData[2 + pp] * JMM[8] + dyData[3 + pp] * (-1 + JMM[9]) +
            dyData[4 + pp] * JMM[13] + dyData[5 + pp] * JMM[18];
    h[3] += -dyData[2 + pp];
    h[4] += 0;
    h[5] += dyData[pp];
    h[0] += dzData[pp] * JMM[10] + dzData[1 + pp] * JMM[11] +
            dzData[2 + pp] * JMM[12] + dzData[3 + pp] * JMM[13] +
            dzData[4 + pp] * (-1 + JMM[14]) + dzData[5 + pp] * JMM[19];
    h[1] += -(dzData[pp] * JMM[6]) - dzData[1 + pp] * JMM[7] -
            dzData[2 + pp] * JMM[8] + dzData[3 + pp] * (1 - JMM[9]) -
            dzData[4 + pp] * JMM[13] - dzData[5 + pp] * JMM[18];
    h[2] += 0;
    h[3] += dzData[1 + pp];
    h[4] += -dzData[pp];
    h[5] += 0;
    h[0] -= h[3] * JMM[6] + h[4] * JMM[10] + h[5] * JMM[15];
    h[1] -= h[3] * JMM[7] + h[4] * JMM[11] + h[5] * JMM[16];
    h[2] -= h[3] * JMM[8] + h[4] * JMM[12] + h[5] * JMM[17];
    dudata[pp + 0] =
        h[2] * (-(JMM[3] * (1 + JMM[2])) + JMM[1] * JMM[4]) +
        h[1] * (JMM[3] * JMM[4] - JMM[1] * (1 + JMM[5])) +
        h[0] * (1 - JMM[4] * JMM[4] + JMM[5] + JMM[2] * (1 + JMM[5]));
    dudata[pp + 1] =
        h[2] * (JMM[3] * JMM[1] - (1 + JMM[0]) * JMM[4]) +
        h[1] * (1 - JMM[3] * JMM[3] + JMM[5] + JMM[0] * (1 + JMM[5])) +
        h[0] * (JMM[4] * JMM[3] - JMM[1] * (1 + JMM[5]));
    dudata[pp + 2] =
        h[2] * (1 - JMM[1] * JMM[1] + JMM[2] + JMM[0] * (1 + JMM[2])) +
        h[1] * (JMM[1] * JMM[3] - (1 + JMM[0]) * JMM[4]) +
        h[0] * (-((1 + JMM[2]) * JMM[3]) + JMM[1] * JMM[4]);
    detC =
        -((1 + JMM[2]) * (-1 + JMM[3] * JMM[3])) +
        (JMM[3] * JMM[1] - JMM[4]) * JMM[4] + JMM[5] + JMM[2] * JMM[5] +
        JMM[0] * (1 + JMM[2] - JMM[4] * JMM[4] + (1 + JMM[2]) * JMM[5]) -
        JMM[1] * (-(JMM[4] * JMM[3]) + JMM[1] * (1 + JMM[5]));
    dudata[pp + 0] /= detC;
    dudata[pp + 1] /= detC;
    dudata[pp + 2] /= detC;
    dudata[pp + 3] = h[3];
    dudata[pp + 4] = h[4];
    dudata[pp + 5] = h[5];
  }
  return;
}
