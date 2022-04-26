/////////////////////////////////////////////////////////////////
/// @file DerivationStencils.h
/// @brief Definition of derivation stencils from order 1 to 13
/////////////////////////////////////////////////////////////////

#pragma once

#include <sundials/sundials_types.h> /* definition of type sunrealtype */

///////////////////////////////////////////////////////
// Stencils with variable nx -- data point dimension //
///////////////////////////////////////////////////////

// Downwind (forward) dfferentiating
inline sunrealtype s1f(sunrealtype *udata, int nx) {
  return -1.0 / 1.0 * udata[-1 * nx] + udata[0];
}
// Upwind (backward) differentiating
inline sunrealtype s1b(sunrealtype *udata, int nx) {
  return -1.0 / 1.0 * udata[0] + udata[1 * nx];
}

inline sunrealtype s2f(const sunrealtype *udata, int nx) {
  return 1.0 / 2.0 * udata[-2 * nx] - 2.0 / 1.0 * udata[-1 * nx] +
         3.0 / 2.0 * udata[0];
}
inline sunrealtype s2c(const sunrealtype *udata, int nx) {
  return -1.0 / 2.0 * udata[-1 * nx] + 0 + 1.0 / 2.0 * udata[1 * nx];
}
inline sunrealtype s2b(const sunrealtype *udata, int nx) {
  return -3.0 / 2.0 * udata[0] + 2.0 / 1.0 * udata[1 * nx] -
         1.0 / 2.0 * udata[2 * nx];
}
inline sunrealtype s3f(const sunrealtype *udata, int nx) {
  return 1.0 / 6.0 * udata[-2 * nx] - 1.0 / 1.0 * udata[-1 * nx] +
         1.0 / 2.0 * udata[0] + 1.0 / 3.0 * udata[1 * nx];
}
inline sunrealtype s3b(sunrealtype *udata, int nx) {
  return -1.0 / 3.0 * udata[-1 * nx] - 1.0 / 2.0 * udata[0] + udata[1 * nx] -
         1.0 / 6.0 * udata[2 * nx];
}
inline sunrealtype s4f(const sunrealtype *udata, int nx) {
  return -1.0 / 12.0 * udata[-3 * nx] + 1.0 / 2.0 * udata[-2 * nx] -
         3.0 / 2.0 * udata[-1 * nx] + 5.0 / 6.0 * udata[0] +
         1.0 / 4.0 * udata[1 * nx];
}
inline sunrealtype s4c(const sunrealtype *udata, int nx) {
  return 1.0 / 12.0 * udata[-2 * nx] - 2.0 / 3.0 * udata[-1 * nx] + 0 +
         2.0 / 3.0 * udata[1 * nx] - 1.0 / 12.0 * udata[2 * nx];
}
inline sunrealtype s4b(const sunrealtype *udata, int nx) {
  return -1.0 / 4.0 * udata[-1 * nx] - 5.0 / 6.0 * udata[0] +
         3.0 / 2.0 * udata[1 * nx] - 1.0 / 2.0 * udata[2 * nx] +
         1.0 / 12.0 * udata[3 * nx];
}
inline sunrealtype s5f(const sunrealtype *udata, int nx) {
  return -1.0 / 30.0 * udata[-3 * nx] + 1.0 / 4.0 * udata[-2 * nx] -
         1.0 / 1.0 * udata[-1 * nx] + 1.0 / 3.0 * udata[0] +
         1.0 / 2.0 * udata[1 * nx] - 1.0 / 20.0 * udata[2 * nx];
}
inline sunrealtype s5b(sunrealtype *udata, int nx) {
  return 1.0 / 20.0 * udata[-2 * nx] - 1.0 / 2.0 * udata[-1 * nx] -
         1.0 / 3.0 * udata[0] + udata[1 * nx] - 1.0 / 4.0 * udata[2 * nx] +
         1.0 / 30.0 * udata[3 * nx];
}
inline sunrealtype s6f(const sunrealtype *udata, int nx) {
  return 1.0 / 60.0 * udata[-4 * nx] - 2.0 / 15.0 * udata[-3 * nx] +
         1.0 / 2.0 * udata[-2 * nx] - 4.0 / 3.0 * udata[-1 * nx] +
         7.0 / 12.0 * udata[0] + 2.0 / 5.0 * udata[1 * nx] -
         1.0 / 30.0 * udata[2 * nx];
}
inline sunrealtype s6c(const sunrealtype *udata, int nx) {
  return -1.0 / 60.0 * udata[-3 * nx] + 3.0 / 20.0 * udata[-2 * nx] -
         3.0 / 4.0 * udata[-1 * nx] + 0 + 3.0 / 4.0 * udata[1 * nx] -
         3.0 / 20.0 * udata[2 * nx] + 1.0 / 60.0 * udata[3 * nx];
}
inline sunrealtype s6b(const sunrealtype *udata, int nx) {
  return 1.0 / 30.0 * udata[-2 * nx] - 2.0 / 5.0 * udata[-1 * nx] -
         7.0 / 12.0 * udata[0] + 4.0 / 3.0 * udata[1 * nx] -
         1.0 / 2.0 * udata[2 * nx] + 2.0 / 15.0 * udata[3 * nx] -
         1.0 / 60.0 * udata[4 * nx];
}
inline sunrealtype s7f(const sunrealtype *udata, int nx) {
  return 1.0 / 140.0 * udata[-4 * nx] - 1.0 / 15.0 * udata[-3 * nx] +
         3.0 / 10.0 * udata[-2 * nx] - 1.0 / 1.0 * udata[-1 * nx] +
         1.0 / 4.0 * udata[0] + 3.0 / 5.0 * udata[1 * nx] -
         1.0 / 10.0 * udata[2 * nx] + 1.0 / 105.0 * udata[3 * nx];
}
inline sunrealtype s7b(sunrealtype *udata, int nx) {
  return -1.0 / 105.0 * udata[-3 * nx] + 1.0 / 10.0 * udata[-2 * nx] -
         3.0 / 5.0 * udata[-1 * nx] - 1.0 / 4.0 * udata[0] + udata[1 * nx] -
         3.0 / 10.0 * udata[2 * nx] + 1.0 / 15.0 * udata[3 * nx] -
         1.0 / 140.0 * udata[4 * nx];
}
inline sunrealtype s8f(const sunrealtype *udata, int nx) {
  return -1.0 / 280.0 * udata[-5 * nx] + 1.0 / 28.0 * udata[-4 * nx] -
         1.0 / 6.0 * udata[-3 * nx] + 1.0 / 2.0 * udata[-2 * nx] -
         5.0 / 4.0 * udata[-1 * nx] + 9.0 / 20.0 * udata[0] +
         1.0 / 2.0 * udata[1 * nx] - 1.0 / 14.0 * udata[2 * nx] +
         1.0 / 168.0 * udata[3 * nx];
}
inline sunrealtype s8c(const sunrealtype *udata, int nx) {
  return 1.0 / 280.0 * udata[-4 * nx] - 4.0 / 105.0 * udata[-3 * nx] +
         1.0 / 5.0 * udata[-2 * nx] - 4.0 / 5.0 * udata[-1 * nx] + 0 +
         4.0 / 5.0 * udata[1 * nx] - 1.0 / 5.0 * udata[2 * nx] +
         4.0 / 105.0 * udata[3 * nx] - 1.0 / 280.0 * udata[4 * nx];
}
inline sunrealtype s8b(const sunrealtype *udata, int nx) {
  return -1.0 / 168.0 * udata[-3 * nx] + 1.0 / 14.0 * udata[-2 * nx] -
         1.0 / 2.0 * udata[-1 * nx] - 9.0 / 20.0 * udata[0] +
         5.0 / 4.0 * udata[1 * nx] - 1.0 / 2.0 * udata[2 * nx] +
         1.0 / 6.0 * udata[3 * nx] - 1.0 / 28.0 * udata[4 * nx] +
         1.0 / 280.0 * udata[5 * nx];
}
inline sunrealtype s9f(const sunrealtype *udata, int nx) {
  return -1.0 / 630.0 * udata[-5 * nx] + 1.0 / 56.0 * udata[-4 * nx] -
         2.0 / 21.0 * udata[-3 * nx] + 1.0 / 3.0 * udata[-2 * nx] -
         1.0 / 1.0 * udata[-1 * nx] + 1.0 / 5.0 * udata[0] +
         2.0 / 3.0 * udata[1 * nx] - 1.0 / 7.0 * udata[2 * nx] +
         1.0 / 42.0 * udata[3 * nx] - 1.0 / 504.0 * udata[4 * nx];
}
inline sunrealtype s9b(sunrealtype *udata, int nx) {
  return 1.0 / 504.0 * udata[-4 * nx] - 1.0 / 42.0 * udata[-3 * nx] +
         1.0 / 7.0 * udata[-2 * nx] - 2.0 / 3.0 * udata[-1 * nx] -
         1.0 / 5.0 * udata[0] + udata[1 * nx] - 1.0 / 3.0 * udata[2 * nx] +
         2.0 / 21.0 * udata[3 * nx] - 1.0 / 56.0 * udata[4 * nx] +
         1.0 / 630.0 * udata[5 * nx];
}
inline sunrealtype s10f(const sunrealtype *udata, int nx) {
  return 1.0 / 1260.0 * udata[-6 * nx] - 1.0 / 105.0 * udata[-5 * nx] +
         3.0 / 56.0 * udata[-4 * nx] - 4.0 / 21.0 * udata[-3 * nx] +
         1.0 / 2.0 * udata[-2 * nx] - 6.0 / 5.0 * udata[-1 * nx] +
         11.0 / 30.0 * udata[0] + 4.0 / 7.0 * udata[1 * nx] -
         3.0 / 28.0 * udata[2 * nx] + 1.0 / 63.0 * udata[3 * nx] -
         1.0 / 840.0 * udata[4 * nx];
}
inline sunrealtype s10c(const sunrealtype *udata, int nx) {
  return -1.0 / 1260.0 * udata[-5 * nx] + 5.0 / 504.0 * udata[-4 * nx] -
         5.0 / 84.0 * udata[-3 * nx] + 5.0 / 21.0 * udata[-2 * nx] -
         5.0 / 6.0 * udata[-1 * nx] + 0 + 5.0 / 6.0 * udata[1 * nx] -
         5.0 / 21.0 * udata[2 * nx] + 5.0 / 84.0 * udata[3 * nx] -
         5.0 / 504.0 * udata[4 * nx] + 1.0 / 1260.0 * udata[5 * nx];
}
inline sunrealtype s10b(const sunrealtype *udata, int nx) {
  return 1.0 / 840.0 * udata[-4 * nx] - 1.0 / 63.0 * udata[-3 * nx] +
         3.0 / 28.0 * udata[-2 * nx] - 4.0 / 7.0 * udata[-1 * nx] -
         11.0 / 30.0 * udata[0] + 6.0 / 5.0 * udata[1 * nx] -
         1.0 / 2.0 * udata[2 * nx] + 4.0 / 21.0 * udata[3 * nx] -
         3.0 / 56.0 * udata[4 * nx] + 1.0 / 105.0 * udata[5 * nx] -
         1.0 / 1260.0 * udata[6 * nx];
}
inline sunrealtype s11f(const sunrealtype *udata, int nx) {
  return 1.0 / 2772.0 * udata[-6 * nx] - 1.0 / 210.0 * udata[-5 * nx] +
         5.0 / 168.0 * udata[-4 * nx] - 5.0 / 42.0 * udata[-3 * nx] +
         5.0 / 14.0 * udata[-2 * nx] - 1.0 / 1.0 * udata[-1 * nx] +
         1.0 / 6.0 * udata[0] + 5.0 / 7.0 * udata[1 * nx] -
         5.0 / 28.0 * udata[2 * nx] + 5.0 / 126.0 * udata[3 * nx] -
         1.0 / 168.0 * udata[4 * nx] + 1.0 / 2310.0 * udata[5 * nx];
}
inline sunrealtype s11b(sunrealtype *udata, int nx) {
  return -1.0 / 2310.0 * udata[-5 * nx] + 1.0 / 168.0 * udata[-4 * nx] -
         5.0 / 126.0 * udata[-3 * nx] + 5.0 / 28.0 * udata[-2 * nx] -
         5.0 / 7.0 * udata[-1 * nx] - 1.0 / 6.0 * udata[0] + udata[1 * nx] -
         5.0 / 14.0 * udata[2 * nx] + 5.0 / 42.0 * udata[3 * nx] -
         5.0 / 168.0 * udata[4 * nx] + 1.0 / 210.0 * udata[5 * nx] -
         1.0 / 2772.0 * udata[6 * nx];
}
inline sunrealtype s12f(const sunrealtype *udata, int nx) {
  return -1.0 / 5544.0 * udata[-7 * nx] + 1.0 / 396.0 * udata[-6 * nx] -
         1.0 / 60.0 * udata[-5 * nx] + 5.0 / 72.0 * udata[-4 * nx] -
         5.0 / 24.0 * udata[-3 * nx] + 1.0 / 2.0 * udata[-2 * nx] -
         7.0 / 6.0 * udata[-1 * nx] + 13.0 / 42.0 * udata[0] +
         5.0 / 8.0 * udata[1 * nx] - 5.0 / 36.0 * udata[2 * nx] +
         1.0 / 36.0 * udata[3 * nx] - 1.0 / 264.0 * udata[4 * nx] +
         1.0 / 3960.0 * udata[5 * nx];
}
inline sunrealtype s12c(const sunrealtype *udata, int nx) {
  return 1.0 / 5544.0 * udata[-6 * nx] - 1.0 / 385.0 * udata[-5 * nx] +
         1.0 / 56.0 * udata[-4 * nx] - 5.0 / 63.0 * udata[-3 * nx] +
         15.0 / 56.0 * udata[-2 * nx] - 6.0 / 7.0 * udata[-1 * nx] + 0 +
         6.0 / 7.0 * udata[1 * nx] - 15.0 / 56.0 * udata[2 * nx] +
         5.0 / 63.0 * udata[3 * nx] - 1.0 / 56.0 * udata[4 * nx] +
         1.0 / 385.0 * udata[5 * nx] - 1.0 / 5544.0 * udata[6 * nx];
}
inline sunrealtype s12b(const sunrealtype *udata, int nx) {
  return -1.0 / 3960.0 * udata[-5 * nx] + 1.0 / 264.0 * udata[-4 * nx] -
         1.0 / 36.0 * udata[-3 * nx] + 5.0 / 36.0 * udata[-2 * nx] -
         5.0 / 8.0 * udata[-1 * nx] - 13.0 / 42.0 * udata[0] +
         7.0 / 6.0 * udata[1 * nx] - 1.0 / 2.0 * udata[2 * nx] +
         5.0 / 24.0 * udata[3 * nx] - 5.0 / 72.0 * udata[4 * nx] +
         1.0 / 60.0 * udata[5 * nx] - 1.0 / 396.0 * udata[6 * nx] +
         1.0 / 5544.0 * udata[7 * nx];
}
inline sunrealtype s13f(const sunrealtype *udata, int nx) {
  return -1.0 / 12012.0 * udata[-7 * nx] + 1.0 / 792.0 * udata[-6 * nx] -
         1.0 / 110.0 * udata[-5 * nx] + 1.0 / 24.0 * udata[-4 * nx] -
         5.0 / 36.0 * udata[-3 * nx] + 3.0 / 8.0 * udata[-2 * nx] -
         1.0 / 1.0 * udata[-1 * nx] + 1.0 / 7.0 * udata[0] +
         3.0 / 4.0 * udata[1 * nx] - 5.0 / 24.0 * udata[2 * nx] +
         1.0 / 18.0 * udata[3 * nx] - 1.0 / 88.0 * udata[4 * nx] +
         1.0 / 660.0 * udata[5 * nx] - 1.0 / 10296.0 * udata[6 * nx];
}
inline sunrealtype s13b(sunrealtype *udata, int nx) {
  return 1.0 / 10296.0 * udata[-6 * nx] - 1.0 / 660.0 * udata[-5 * nx] +
         1.0 / 88.0 * udata[-4 * nx] - 1.0 / 18.0 * udata[-3 * nx] +
         5.0 / 24.0 * udata[-2 * nx] - 3.0 / 4.0 * udata[-1 * nx] -
         1.0 / 7.0 * udata[0] + udata[1 * nx] - 3.0 / 8.0 * udata[2 * nx] +
         5.0 / 36.0 * udata[3 * nx] - 1.0 / 24.0 * udata[4 * nx] +
         1.0 / 110.0 * udata[5 * nx] - 1.0 / 792.0 * udata[6 * nx] +
         1.0 / 12012.0 * udata[7 * nx];
}

////////////////////////////////
// Stencils with nx fixed to 6//
////////////////////////////////

inline sunrealtype s1f(sunrealtype *udata) { return s1f(udata, 6); }
inline sunrealtype s1b(sunrealtype *udata) { return s1b(udata, 6); }
inline sunrealtype s2f(sunrealtype *udata) { return s2f(udata, 6); }
inline sunrealtype s2c(sunrealtype *udata) { return s2c(udata, 6); }
inline sunrealtype s2b(sunrealtype *udata) { return s2b(udata, 6); }
inline sunrealtype s3f(sunrealtype *udata) { return s3f(udata, 6); }
inline sunrealtype s3b(sunrealtype *udata) { return s3b(udata, 6); }
inline sunrealtype s4f(sunrealtype *udata) { return s4f(udata, 6); }
inline sunrealtype s4c(sunrealtype *udata) { return s4c(udata, 6); }
inline sunrealtype s4b(sunrealtype *udata) { return s4b(udata, 6); }
inline sunrealtype s5f(sunrealtype *udata) { return s5f(udata, 6); }
inline sunrealtype s5b(sunrealtype *udata) { return s5b(udata, 6); }
inline sunrealtype s6f(sunrealtype *udata) { return s6f(udata, 6); }
inline sunrealtype s6c(sunrealtype *udata) { return s6c(udata, 6); }
inline sunrealtype s6b(sunrealtype *udata) { return s6b(udata, 6); }
inline sunrealtype s7f(sunrealtype *udata) { return s7f(udata, 6); }
inline sunrealtype s7b(sunrealtype *udata) { return s7b(udata, 6); }
inline sunrealtype s8f(sunrealtype *udata) { return s8f(udata, 6); }
inline sunrealtype s8c(sunrealtype *udata) { return s8c(udata, 6); }
inline sunrealtype s8b(sunrealtype *udata) { return s8b(udata, 6); }
inline sunrealtype s9f(sunrealtype *udata) { return s9f(udata, 6); }
inline sunrealtype s9b(sunrealtype *udata) { return s9b(udata, 6); }
inline sunrealtype s10f(sunrealtype *udata) { return s10f(udata, 6); }
inline sunrealtype s10c(sunrealtype *udata) { return s10c(udata, 6); }
inline sunrealtype s10b(sunrealtype *udata) { return s10b(udata, 6); }
inline sunrealtype s11f(sunrealtype *udata) { return s11f(udata, 6); }
inline sunrealtype s11b(sunrealtype *udata) { return s11b(udata, 6); }
inline sunrealtype s12f(sunrealtype *udata) { return s12f(udata, 6); }
inline sunrealtype s12c(sunrealtype *udata) { return s12c(udata, 6); }
inline sunrealtype s12b(sunrealtype *udata) { return s12b(udata, 6); }
inline sunrealtype s13f(sunrealtype *udata) { return s13f(udata, 6); }
inline sunrealtype s13b(sunrealtype *udata) { return s13b(udata, 6); }

