/////////////////////////////////////////////////////////////////
/// @file DerivationStencils.h
/// @brief Definition of derivation stencils from order 1 to 13
/////////////////////////////////////////////////////////////////

#pragma once

#include <sundials/sundials_types.h> /* definition of type sunrealtype */

////////////////////////////////////////////////////////
// Stencils with variable dPD -- data point dimension //
////////////////////////////////////////////////////////

// Downwind (forward) dfferentiating
const inline sunrealtype s1f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 1.0 * udata[-1 * dPD] + udata[0];
}
// Upwind (backward) differentiating
const inline sunrealtype s1b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 1.0 * udata[0] + udata[1 * dPD];
}

const inline sunrealtype s2f(sunrealtype const *udata, const int dPD) {
  return 1.0 / 2.0 * udata[-2 * dPD] - 2.0 / 1.0 * udata[-1 * dPD] +
         3.0 / 2.0 * udata[0];
}
const inline sunrealtype s2c(sunrealtype const *udata, const int dPD) {
  return -1.0 / 2.0 * udata[-1 * dPD] + 0 + 1.0 / 2.0 * udata[1 * dPD];
}
const inline sunrealtype s2b(sunrealtype const *udata, const int dPD) {
  return -3.0 / 2.0 * udata[0] + 2.0 / 1.0 * udata[1 * dPD] -
         1.0 / 2.0 * udata[2 * dPD];
}
const inline sunrealtype s3f(sunrealtype const *udata, const int dPD) {
  return 1.0 / 6.0 * udata[-2 * dPD] - 1.0 / 1.0 * udata[-1 * dPD] +
         1.0 / 2.0 * udata[0] + 1.0 / 3.0 * udata[1 * dPD];
}
const inline sunrealtype s3b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 3.0 * udata[-1 * dPD] - 1.0 / 2.0 * udata[0] + udata[1 * dPD] -
         1.0 / 6.0 * udata[2 * dPD];
}
const inline sunrealtype s4f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 12.0 * udata[-3 * dPD] + 1.0 / 2.0 * udata[-2 * dPD] -
         3.0 / 2.0 * udata[-1 * dPD] + 5.0 / 6.0 * udata[0] +
         1.0 / 4.0 * udata[1 * dPD];
}
const inline sunrealtype s4c(sunrealtype const *udata, const int dPD) {
  return 1.0 / 12.0 * udata[-2 * dPD] - 2.0 / 3.0 * udata[-1 * dPD] + 0 +
         2.0 / 3.0 * udata[1 * dPD] - 1.0 / 12.0 * udata[2 * dPD];
}
const inline sunrealtype s4b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 4.0 * udata[-1 * dPD] - 5.0 / 6.0 * udata[0] +
         3.0 / 2.0 * udata[1 * dPD] - 1.0 / 2.0 * udata[2 * dPD] +
         1.0 / 12.0 * udata[3 * dPD];
}
const inline sunrealtype s5f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 30.0 * udata[-3 * dPD] + 1.0 / 4.0 * udata[-2 * dPD] -
         1.0 / 1.0 * udata[-1 * dPD] + 1.0 / 3.0 * udata[0] +
         1.0 / 2.0 * udata[1 * dPD] - 1.0 / 20.0 * udata[2 * dPD];
}
const inline sunrealtype s5b(sunrealtype const *udata, const int dPD) {
  return 1.0 / 20.0 * udata[-2 * dPD] - 1.0 / 2.0 * udata[-1 * dPD] -
         1.0 / 3.0 * udata[0] + udata[1 * dPD] - 1.0 / 4.0 * udata[2 * dPD] +
         1.0 / 30.0 * udata[3 * dPD];
}
const inline sunrealtype s6f(sunrealtype const *udata, const int dPD) {
  return 1.0 / 60.0 * udata[-4 * dPD] - 2.0 / 15.0 * udata[-3 * dPD] +
         1.0 / 2.0 * udata[-2 * dPD] - 4.0 / 3.0 * udata[-1 * dPD] +
         7.0 / 12.0 * udata[0] + 2.0 / 5.0 * udata[1 * dPD] -
         1.0 / 30.0 * udata[2 * dPD];
}
const inline sunrealtype s6c(sunrealtype const *udata, const int dPD) {
  return -1.0 / 60.0 * udata[-3 * dPD] + 3.0 / 20.0 * udata[-2 * dPD] -
         3.0 / 4.0 * udata[-1 * dPD] + 0 + 3.0 / 4.0 * udata[1 * dPD] -
         3.0 / 20.0 * udata[2 * dPD] + 1.0 / 60.0 * udata[3 * dPD];
}
const inline sunrealtype s6b(sunrealtype const *udata, const int dPD) {
  return 1.0 / 30.0 * udata[-2 * dPD] - 2.0 / 5.0 * udata[-1 * dPD] -
         7.0 / 12.0 * udata[0] + 4.0 / 3.0 * udata[1 * dPD] -
         1.0 / 2.0 * udata[2 * dPD] + 2.0 / 15.0 * udata[3 * dPD] -
         1.0 / 60.0 * udata[4 * dPD];
}
const inline sunrealtype s7f(sunrealtype const *udata, const int dPD) {
  return 1.0 / 140.0 * udata[-4 * dPD] - 1.0 / 15.0 * udata[-3 * dPD] +
         3.0 / 10.0 * udata[-2 * dPD] - 1.0 / 1.0 * udata[-1 * dPD] +
         1.0 / 4.0 * udata[0] + 3.0 / 5.0 * udata[1 * dPD] -
         1.0 / 10.0 * udata[2 * dPD] + 1.0 / 105.0 * udata[3 * dPD];
}
const inline sunrealtype s7b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 105.0 * udata[-3 * dPD] + 1.0 / 10.0 * udata[-2 * dPD] -
         3.0 / 5.0 * udata[-1 * dPD] - 1.0 / 4.0 * udata[0] + udata[1 * dPD] -
         3.0 / 10.0 * udata[2 * dPD] + 1.0 / 15.0 * udata[3 * dPD] -
         1.0 / 140.0 * udata[4 * dPD];
}
const inline sunrealtype s8f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 280.0 * udata[-5 * dPD] + 1.0 / 28.0 * udata[-4 * dPD] -
         1.0 / 6.0 * udata[-3 * dPD] + 1.0 / 2.0 * udata[-2 * dPD] -
         5.0 / 4.0 * udata[-1 * dPD] + 9.0 / 20.0 * udata[0] +
         1.0 / 2.0 * udata[1 * dPD] - 1.0 / 14.0 * udata[2 * dPD] +
         1.0 / 168.0 * udata[3 * dPD];
}
const inline sunrealtype s8c(sunrealtype const *udata, const int dPD) {
  return 1.0 / 280.0 * udata[-4 * dPD] - 4.0 / 105.0 * udata[-3 * dPD] +
         1.0 / 5.0 * udata[-2 * dPD] - 4.0 / 5.0 * udata[-1 * dPD] + 0 +
         4.0 / 5.0 * udata[1 * dPD] - 1.0 / 5.0 * udata[2 * dPD] +
         4.0 / 105.0 * udata[3 * dPD] - 1.0 / 280.0 * udata[4 * dPD];
}
const inline sunrealtype s8b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 168.0 * udata[-3 * dPD] + 1.0 / 14.0 * udata[-2 * dPD] -
         1.0 / 2.0 * udata[-1 * dPD] - 9.0 / 20.0 * udata[0] +
         5.0 / 4.0 * udata[1 * dPD] - 1.0 / 2.0 * udata[2 * dPD] +
         1.0 / 6.0 * udata[3 * dPD] - 1.0 / 28.0 * udata[4 * dPD] +
         1.0 / 280.0 * udata[5 * dPD];
}
const inline sunrealtype s9f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 630.0 * udata[-5 * dPD] + 1.0 / 56.0 * udata[-4 * dPD] -
         2.0 / 21.0 * udata[-3 * dPD] + 1.0 / 3.0 * udata[-2 * dPD] -
         1.0 / 1.0 * udata[-1 * dPD] + 1.0 / 5.0 * udata[0] +
         2.0 / 3.0 * udata[1 * dPD] - 1.0 / 7.0 * udata[2 * dPD] +
         1.0 / 42.0 * udata[3 * dPD] - 1.0 / 504.0 * udata[4 * dPD];
}
const inline sunrealtype s9b(sunrealtype const *udata, const int dPD) {
  return 1.0 / 504.0 * udata[-4 * dPD] - 1.0 / 42.0 * udata[-3 * dPD] +
         1.0 / 7.0 * udata[-2 * dPD] - 2.0 / 3.0 * udata[-1 * dPD] -
         1.0 / 5.0 * udata[0] + udata[1 * dPD] - 1.0 / 3.0 * udata[2 * dPD] +
         2.0 / 21.0 * udata[3 * dPD] - 1.0 / 56.0 * udata[4 * dPD] +
         1.0 / 630.0 * udata[5 * dPD];
}
const inline sunrealtype s10f(sunrealtype const *udata, const int dPD) {
  return 1.0 / 1260.0 * udata[-6 * dPD] - 1.0 / 105.0 * udata[-5 * dPD] +
         3.0 / 56.0 * udata[-4 * dPD] - 4.0 / 21.0 * udata[-3 * dPD] +
         1.0 / 2.0 * udata[-2 * dPD] - 6.0 / 5.0 * udata[-1 * dPD] +
         11.0 / 30.0 * udata[0] + 4.0 / 7.0 * udata[1 * dPD] -
         3.0 / 28.0 * udata[2 * dPD] + 1.0 / 63.0 * udata[3 * dPD] -
         1.0 / 840.0 * udata[4 * dPD];
}
const inline sunrealtype s10c(sunrealtype const *udata, const int dPD) {
  return -1.0 / 1260.0 * udata[-5 * dPD] + 5.0 / 504.0 * udata[-4 * dPD] -
         5.0 / 84.0 * udata[-3 * dPD] + 5.0 / 21.0 * udata[-2 * dPD] -
         5.0 / 6.0 * udata[-1 * dPD] + 0 + 5.0 / 6.0 * udata[1 * dPD] -
         5.0 / 21.0 * udata[2 * dPD] + 5.0 / 84.0 * udata[3 * dPD] -
         5.0 / 504.0 * udata[4 * dPD] + 1.0 / 1260.0 * udata[5 * dPD];
}
const inline sunrealtype s10b(sunrealtype const *udata, const int dPD) {
  return 1.0 / 840.0 * udata[-4 * dPD] - 1.0 / 63.0 * udata[-3 * dPD] +
         3.0 / 28.0 * udata[-2 * dPD] - 4.0 / 7.0 * udata[-1 * dPD] -
         11.0 / 30.0 * udata[0] + 6.0 / 5.0 * udata[1 * dPD] -
         1.0 / 2.0 * udata[2 * dPD] + 4.0 / 21.0 * udata[3 * dPD] -
         3.0 / 56.0 * udata[4 * dPD] + 1.0 / 105.0 * udata[5 * dPD] -
         1.0 / 1260.0 * udata[6 * dPD];
}
const inline sunrealtype s11f(sunrealtype const *udata, const int dPD) {
  return 1.0 / 2772.0 * udata[-6 * dPD] - 1.0 / 210.0 * udata[-5 * dPD] +
         5.0 / 168.0 * udata[-4 * dPD] - 5.0 / 42.0 * udata[-3 * dPD] +
         5.0 / 14.0 * udata[-2 * dPD] - 1.0 / 1.0 * udata[-1 * dPD] +
         1.0 / 6.0 * udata[0] + 5.0 / 7.0 * udata[1 * dPD] -
         5.0 / 28.0 * udata[2 * dPD] + 5.0 / 126.0 * udata[3 * dPD] -
         1.0 / 168.0 * udata[4 * dPD] + 1.0 / 2310.0 * udata[5 * dPD];
}
const inline sunrealtype s11b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 2310.0 * udata[-5 * dPD] + 1.0 / 168.0 * udata[-4 * dPD] -
         5.0 / 126.0 * udata[-3 * dPD] + 5.0 / 28.0 * udata[-2 * dPD] -
         5.0 / 7.0 * udata[-1 * dPD] - 1.0 / 6.0 * udata[0] + udata[1 * dPD] -
         5.0 / 14.0 * udata[2 * dPD] + 5.0 / 42.0 * udata[3 * dPD] -
         5.0 / 168.0 * udata[4 * dPD] + 1.0 / 210.0 * udata[5 * dPD] -
         1.0 / 2772.0 * udata[6 * dPD];
}
const inline sunrealtype s12f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 5544.0 * udata[-7 * dPD] + 1.0 / 396.0 * udata[-6 * dPD] -
         1.0 / 60.0 * udata[-5 * dPD] + 5.0 / 72.0 * udata[-4 * dPD] -
         5.0 / 24.0 * udata[-3 * dPD] + 1.0 / 2.0 * udata[-2 * dPD] -
         7.0 / 6.0 * udata[-1 * dPD] + 13.0 / 42.0 * udata[0] +
         5.0 / 8.0 * udata[1 * dPD] - 5.0 / 36.0 * udata[2 * dPD] +
         1.0 / 36.0 * udata[3 * dPD] - 1.0 / 264.0 * udata[4 * dPD] +
         1.0 / 3960.0 * udata[5 * dPD];
}
const inline sunrealtype s12c(sunrealtype const *udata, const int dPD) {
  return 1.0 / 5544.0 * udata[-6 * dPD] - 1.0 / 385.0 * udata[-5 * dPD] +
         1.0 / 56.0 * udata[-4 * dPD] - 5.0 / 63.0 * udata[-3 * dPD] +
         15.0 / 56.0 * udata[-2 * dPD] - 6.0 / 7.0 * udata[-1 * dPD] + 0 +
         6.0 / 7.0 * udata[1 * dPD] - 15.0 / 56.0 * udata[2 * dPD] +
         5.0 / 63.0 * udata[3 * dPD] - 1.0 / 56.0 * udata[4 * dPD] +
         1.0 / 385.0 * udata[5 * dPD] - 1.0 / 5544.0 * udata[6 * dPD];
}
const inline sunrealtype s12b(sunrealtype const *udata, const int dPD) {
  return -1.0 / 3960.0 * udata[-5 * dPD] + 1.0 / 264.0 * udata[-4 * dPD] -
         1.0 / 36.0 * udata[-3 * dPD] + 5.0 / 36.0 * udata[-2 * dPD] -
         5.0 / 8.0 * udata[-1 * dPD] - 13.0 / 42.0 * udata[0] +
         7.0 / 6.0 * udata[1 * dPD] - 1.0 / 2.0 * udata[2 * dPD] +
         5.0 / 24.0 * udata[3 * dPD] - 5.0 / 72.0 * udata[4 * dPD] +
         1.0 / 60.0 * udata[5 * dPD] - 1.0 / 396.0 * udata[6 * dPD] +
         1.0 / 5544.0 * udata[7 * dPD];
}
const inline sunrealtype s13f(sunrealtype const *udata, const int dPD) {
  return -1.0 / 12012.0 * udata[-7 * dPD] + 1.0 / 792.0 * udata[-6 * dPD] -
         1.0 / 110.0 * udata[-5 * dPD] + 1.0 / 24.0 * udata[-4 * dPD] -
         5.0 / 36.0 * udata[-3 * dPD] + 3.0 / 8.0 * udata[-2 * dPD] -
         1.0 / 1.0 * udata[-1 * dPD] + 1.0 / 7.0 * udata[0] +
         3.0 / 4.0 * udata[1 * dPD] - 5.0 / 24.0 * udata[2 * dPD] +
         1.0 / 18.0 * udata[3 * dPD] - 1.0 / 88.0 * udata[4 * dPD] +
         1.0 / 660.0 * udata[5 * dPD] - 1.0 / 10296.0 * udata[6 * dPD];
}
const inline sunrealtype s13b(sunrealtype const *udata, const int dPD) {
  return 1.0 / 10296.0 * udata[-6 * dPD] - 1.0 / 660.0 * udata[-5 * dPD] +
         1.0 / 88.0 * udata[-4 * dPD] - 1.0 / 18.0 * udata[-3 * dPD] +
         5.0 / 24.0 * udata[-2 * dPD] - 3.0 / 4.0 * udata[-1 * dPD] -
         1.0 / 7.0 * udata[0] + udata[1 * dPD] - 3.0 / 8.0 * udata[2 * dPD] +
         5.0 / 36.0 * udata[3 * dPD] - 1.0 / 24.0 * udata[4 * dPD] +
         1.0 / 110.0 * udata[5 * dPD] - 1.0 / 792.0 * udata[6 * dPD] +
         1.0 / 12012.0 * udata[7 * dPD];
}

//////////////////////////////////
// Stencils with dPD fixed to 6 //
//////////////////////////////////

const inline sunrealtype s1f(sunrealtype const *udata) { return s1f(udata, 6); }
const inline sunrealtype s1b(sunrealtype const *udata) { return s1b(udata, 6); }
const inline sunrealtype s2f(sunrealtype const *udata) { return s2f(udata, 6); }
const inline sunrealtype s2c(sunrealtype const *udata) { return s2c(udata, 6); }
const inline sunrealtype s2b(sunrealtype const *udata) { return s2b(udata, 6); }
const inline sunrealtype s3f(sunrealtype const *udata) { return s3f(udata, 6); }
const inline sunrealtype s3b(sunrealtype const *udata) { return s3b(udata, 6); }
const inline sunrealtype s4f(sunrealtype const *udata) { return s4f(udata, 6); }
const inline sunrealtype s4c(sunrealtype const *udata) { return s4c(udata, 6); }
const inline sunrealtype s4b(sunrealtype const *udata) { return s4b(udata, 6); }
const inline sunrealtype s5f(sunrealtype const *udata) { return s5f(udata, 6); }
const inline sunrealtype s5b(sunrealtype const *udata) { return s5b(udata, 6); }
const inline sunrealtype s6f(sunrealtype const *udata) { return s6f(udata, 6); }
const inline sunrealtype s6c(sunrealtype const *udata) { return s6c(udata, 6); }
const inline sunrealtype s6b(sunrealtype const *udata) { return s6b(udata, 6); }
const inline sunrealtype s7f(sunrealtype const *udata) { return s7f(udata, 6); }
const inline sunrealtype s7b(sunrealtype const *udata) { return s7b(udata, 6); }
const inline sunrealtype s8f(sunrealtype const *udata) { return s8f(udata, 6); }
const inline sunrealtype s8c(sunrealtype const *udata) { return s8c(udata, 6); }
const inline sunrealtype s8b(sunrealtype const *udata) { return s8b(udata, 6); }
const inline sunrealtype s9f(sunrealtype const *udata) { return s9f(udata, 6); }
const inline sunrealtype s9b(sunrealtype const *udata) { return s9b(udata, 6); }
const inline sunrealtype s10f(sunrealtype const *udata)
{ return s10f(udata, 6); }
const inline sunrealtype s10c(sunrealtype const *udata)
{ return s10c(udata, 6); }
const inline sunrealtype s10b(sunrealtype const *udata)
{ return s10b(udata, 6); }
const inline sunrealtype s11f(sunrealtype const *udata)
{ return s11f(udata, 6); }
const inline sunrealtype s11b(sunrealtype const *udata)
{ return s11b(udata, 6); }
const inline sunrealtype s12f(sunrealtype const *udata)
{ return s12f(udata, 6); }
const inline sunrealtype s12c(sunrealtype const *udata)
{ return s12c(udata, 6); }
const inline sunrealtype s12b(sunrealtype const *udata)
{ return s12b(udata, 6); }
const inline sunrealtype s13f(sunrealtype const *udata)
{ return s13f(udata, 6); }
const inline sunrealtype s13b(sunrealtype const *udata)
{ return s13b(udata, 6); }

