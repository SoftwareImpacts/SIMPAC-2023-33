//////////////////////////////////////////////////////////////////////////////////////
/// @file ICSetters.cpp
/// @brief Implementation of the plane wave and Gaussian wave packets in 1D, 2D,
/// 3D
//////////////////////////////////////////////////////////////////////////////////////

#include "ICSetters.h"

#include <math.h>

/** PlaneWave1D construction with */
PlaneWave1D::PlaneWave1D(vector<sunrealtype> k, vector<sunrealtype> p,
                         vector<sunrealtype> phi) {
  kx = k[0]; /** - wavevectors \f$ k_x \f$ */
  ky = k[1]; /** - \f$ k_y \f$ */
  kz = k[2]; /** - \f$ k_z \f$ normalized to \f$ 1/\lambda \f$ */
  // Amplitude bug: lower by factor 3
  px = p[0] / 3; /** - amplitude (polarization) in x-direction \f$ p_x \f$ */
  py = p[1] / 3; /** - amplitude (polarization) in y-direction \f$ p_y \f$ */
  pz = p[2] / 3; /** - amplitude (polarization) in z-direction \f$ p_z \f$ */
  phix = phi[0]; /** - phase shift in x-direction \f$ \phi_x \f$ */
  phiy = phi[1]; /** - phase shift in y-direction \f$ \phi_y \f$ */
  phiz = phi[2]; /** - phase shift in z-direction \f$ \phi_z \f$ */
}

/** PlaneWave1D implementation in space */
void PlaneWave1D::addToSpace(const sunrealtype x, const sunrealtype y, const sunrealtype z,
                             sunrealtype *pTo6Space) const {
  const sunrealtype wavelength =
      sqrt(kx * kx + ky * ky + kz * kz); /* \f$ 1/\lambda \f$ */
  const sunrealtype kScalarX = (kx * x + ky * y + kz * z) * 2 *
                         numbers::pi; /* \f$ 2\pi \ \vec{k} \cdot \vec{x} \f$ */
  // Plane wave definition
  const array<sunrealtype, 3> E{{                             /* E-field vector */
                           px * cos(kScalarX - phix),   /* \f$ E_x \f$ */
                           py * cos(kScalarX - phiy),   /* \f$ E_y \f$ */
                           pz * cos(kScalarX - phiz)}}; /* \f$ E_z \f$ */
  // Put E-field into space
  pTo6Space[0] += E[0];
  pTo6Space[1] += E[1];
  pTo6Space[2] += E[2];
  // and B-field
  pTo6Space[3] += (ky * E[2] - kz * E[1]) / wavelength;
  pTo6Space[4] += (kz * E[0] - kx * E[2]) / wavelength;
  pTo6Space[5] += (kx * E[1] - ky * E[0]) / wavelength;
}

/** PlaneWave2D construction with */
PlaneWave2D::PlaneWave2D(vector<sunrealtype> k, vector<sunrealtype> p,
                         vector<sunrealtype> phi) {
  kx = k[0]; /** - wavevectors \f$ k_x \f$ */
  ky = k[1]; /** - \f$ k_y \f$ */
  kz = k[2]; /** - \f$ k_z \f$ normalized to \f$ 1/\lambda \f$*/
  // Amplitude bug: lower by factor 9
  px = p[0] / 9; /** - amplitude (polarization) in x-direction \f$ p_x \f$ */
  py = p[1] / 9; /** - amplitude (polarization) in y-direction \f$ p_y \f$ */
  pz = p[2] / 9; /** - amplitude (polarization) in z-direction \f$ p_z \f$ */
  phix = phi[0]; /** - phase shift in x-direction \f$ \phi_x \f$ */
  phiy = phi[1]; /** - phase shift in y-direction \f$ \phi_y \f$ */
  phiz = phi[2]; /** - phase shift in z-direction \f$ \phi_z \f$ */
}

/** PlaneWave2D implementation in space */
void PlaneWave2D::addToSpace(const sunrealtype x, const sunrealtype y, const sunrealtype z,
                             sunrealtype *pTo6Space) const {
  const sunrealtype wavelength =
      sqrt(kx * kx + ky * ky + kz * kz); /* \f$ 1/\lambda \f$ */
  const sunrealtype kScalarX = (kx * x + ky * y + kz * z) * 2 *
                         numbers::pi; /* \f$ 2\pi \ \vec{k} \cdot \vec{x} \f$ */
  // Plane wave definition
  const array<sunrealtype, 3> E{{                             /* E-field vector */
                           px * cos(kScalarX - phix),   /* \f$ E_x \f$ */
                           py * cos(kScalarX - phiy),   /* \f$ E_y \f$ */
                           pz * cos(kScalarX - phiz)}}; /* \f$ E_z \f$ */
  // Put E-field into space
  pTo6Space[0] += E[0];
  pTo6Space[1] += E[1];
  pTo6Space[2] += E[2];
  // and B-field
  pTo6Space[3] += (ky * E[2] - kz * E[1]) / wavelength;
  pTo6Space[4] += (kz * E[0] - kx * E[2]) / wavelength;
  pTo6Space[5] += (kx * E[1] - ky * E[0]) / wavelength;
}

/** PlaneWave3D construction with */
PlaneWave3D::PlaneWave3D(vector<sunrealtype> k, vector<sunrealtype> p,
                         vector<sunrealtype> phi) {
  kx = k[0];     /** - wavevectors \f$ k_x \f$ */
  ky = k[1];     /** - \f$ k_y \f$ */
  kz = k[2];     /** - \f$ k_z \f$ normalized to \f$ 1/\lambda \f$ */
  px = p[0];     /** - amplitude (polarization) in x-direction \f$ p_x \f$ */
  py = p[1];     /** - amplitude (polarization) in y-direction \f$ p_y \f$ */
  pz = p[2];     /** - amplitude (polarization) in z-direction \f$ p_z \f$ */
  phix = phi[0]; /** - phase shift in x-direction \f$ \phi_x \f$ */
  phiy = phi[1]; /** - phase shift in y-direction \f$ \phi_y \f$ */
  phiz = phi[2]; /** - phase shift in z-direction \f$ \phi_z \f$ */
}

/** PlaneWave3D implementation in space */
void PlaneWave3D::addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                             sunrealtype *pTo6Space) const {
  const sunrealtype wavelength =
      sqrt(kx * kx + ky * ky + kz * kz); /* \f$ 1/\lambda \f$ */
  const sunrealtype kScalarX = (kx * x + ky * y + kz * z) * 2 *
                         numbers::pi; /* \f$ 2\pi \ \vec{k} \cdot \vec{x} \f$ */
  // Plane wave definition
  const array<sunrealtype, 3> E{{/* E-field vector \f$ \vec{E}\f$*/
                           px * cos(kScalarX - phix),   /* \f$ E_x \f$ */
                           py * cos(kScalarX - phiy),   /* \f$ E_y \f$ */
                           pz * cos(kScalarX - phiz)}}; /* \f$ E_z \f$ */
  // Put E-field into space
  pTo6Space[0] += E[0];
  pTo6Space[1] += E[1];
  pTo6Space[2] += E[2];
  // and B-field
  pTo6Space[3] += (ky * E[2] - kz * E[1]) / wavelength;
  pTo6Space[4] += (kz * E[0] - kx * E[2]) / wavelength;
  pTo6Space[5] += (kx * E[1] - ky * E[0]) / wavelength;
}

/** Gauss1D construction with */
Gauss1D::Gauss1D(vector<sunrealtype> k, vector<sunrealtype> p,
                 vector<sunrealtype> xo, sunrealtype phig_,
                 vector<sunrealtype> phi) {
  kx = k[0];     /** - wavevectors \f$ k_x \f$ */
  ky = k[1];     /** - \f$ k_y \f$ */
  kz = k[2];     /** - \f$ k_z \f$ normalized to \f$ 1/\lambda \f$*/
  px = p[0];     /** - amplitude (polarization) in x-direction */
  py = p[1];     /** - amplitude (polarization) in y-direction */
  pz = p[2];     /** - amplitude (polarization) in z-direction */
  phix = phi[0]; /** - phase shift in x-direction */
  phiy = phi[1]; /** - phase shift in y-direction */
  phiz = phi[2]; /** - phase shift in z-direction */
  phig = phig_;  /** - width */
  x0x = xo[0];   /** - shift from origin in x-direction*/
  x0y = xo[1];   /** - shift from origin in y-direction*/
  x0z = xo[2];   /** - shift from origin in z-direction*/
}

/** Gauss1D implementation in space */
void Gauss1D::addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                         sunrealtype *pTo6Space) const {
  const sunrealtype wavelength =
      sqrt(kx * kx + ky * ky + kz * kz); /* \f$ 1/\lambda \f$ */
  x = x - x0x; /* x-coordinate minus shift from origin */
  y = y - x0y; /* y-coordinate minus shift from origin */
  z = z - x0z; /* z-coordinate minus shift from origin */
  const sunrealtype kScalarX = (kx * x + ky * y + kz * z) * 2 *
                         numbers::pi; /* \f$ 2\pi \ \vec{k} \cdot \vec{x} \f$ */
  const sunrealtype envelopeAmp =
      exp(-(x * x + y * y + z * z) / phig / phig); /* enveloping Gauss shape */
  // Gaussian wave definition
  const array<sunrealtype, 3> E{
      {                                           /* E-field vector */
       px * cos(kScalarX - phix) * envelopeAmp,   /* \f$ E_x \f$ */
       py * cos(kScalarX - phiy) * envelopeAmp,   /* \f$ E_y \f$ */
       pz * cos(kScalarX - phiz) * envelopeAmp}}; /* \f$ E_z \f$ */
  // Put E-field into space
  pTo6Space[0] += E[0];
  pTo6Space[1] += E[1];
  pTo6Space[2] += E[2];
  // and B-field
  pTo6Space[3] += (ky * E[2] - kz * E[1]) / wavelength;
  pTo6Space[4] += (kz * E[0] - kx * E[2]) / wavelength;
  pTo6Space[5] += (kx * E[1] - ky * E[0]) / wavelength;
}

/** Gauss2D construction with */
Gauss2D::Gauss2D(vector<sunrealtype> dis_, vector<sunrealtype> axis_,
                 sunrealtype Amp_, sunrealtype phip_, sunrealtype w0_,
                 sunrealtype zr_, sunrealtype Ph0_, sunrealtype PhA_) {
  dis = dis_;           /** - center it approaches */
  axis = axis_;         /** - direction form where it comes */
  Amp = Amp_;           /** - amplitude */
  phip = phip_;         /** - polarization rotation from TE-mode */
  w0 = w0_;             /** - taille */
  zr = zr_;             /** - Rayleigh length */
  Ph0 = Ph0_;           /** - beam center */
  PhA = PhA_;           /** - beam length */
  A1 = Amp * cos(phip); // amplitude in z-direction
  A2 = Amp * sin(phip); // amplitude on xy-plane
  lambda = numbers::pi * w0 * w0 / zr; // formula for wavelength
}

void Gauss2D::addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                         sunrealtype *pTo6Space) const {
  //\f$ \vec{x} = \vec{x}_0-\vec{dis} \f$ // coordinates minus distance to
  //origin
  x -= dis[0];
  y -= dis[1];
  // z-=dis[2];
  z = NAN;
  //  \f$ z_g = \vec{x}\cdot\vec{e}_g \f$ projection on propagation axis
  const sunrealtype zg =
      x * axis[0] + y * axis[1]; //+z*axis[2];  // =z-z0 -> propagation
                                 //direction, minus origin
  // \f$ r = \sqrt{\vec{x}^2 -z_g^2} \f$ -> pythagoras of radius minus
  // projection on prop axis
  const sunrealtype r = sqrt((x * x + y * y /*+z*z*/) -
                       zg * zg); // radial distance to propagation axis
  // \f$  w(z) = w0\sqrt{1+(z_g/z_R)^2} \f$
  const sunrealtype wz = w0 * sqrt(1 + (zg * zg / zr / zr)); // waist at position z
  // \f$ g(z) = atan(z_g/z_r) \f$
  const sunrealtype gz = atan(zg / zr); // Gouy phase
  // \f$ R(z) = z_g*(1+(z_r/z_g)^2) \f$
  sunrealtype Rz = NAN; // beam curvature
  if (zg != 0)
    Rz = zg * (1 + (zr * zr / zg / zg));
  else
    Rz = 1e308;
  // wavenumber \f$ k = 2\pi/\lambda \f$
  const sunrealtype k = 2 * numbers::pi / lambda;
  // \f$ \Phi_F = kr^2/(2*R(z))+g(z)-kz_g \f$
  const sunrealtype PhF =
      -k * r * r / (2 * Rz) + gz - k * zg; // to be inserted into cosine
  // \f$ G = \sqrt{w_0/w_z}\e^{-(r/w(z))^2}\e^{(zg-Ph0)^2/PhA^2}\cos(PhF) \f$ –
  // CVode is a diva, no chance to remove the square in the second exponential
  // -> h too small
  const sunrealtype G2D = sqrt(w0 / wz) * exp(-r * r / wz / wz) *
                    exp(-(zg - Ph0) * (zg - Ph0) / PhA / PhA) *
                    cos(PhF); // gauss shape
  // \f$ c_\alpha =\vec{e}_x\cdot\vec{axis} \f$
  // projection components; do like this for CVode convergence -> otherwise
  // results in machine error values for non-existant field components if
  // axis[0] and axis[1] are given
  const sunrealtype ca =
      axis[0]; // x-component of propagation axis which is given as parameter
  const sunrealtype sa = sqrt(1 - ca * ca); // no z-component for 2D propagation
  // E-field to space: polarization in xy-plane (A2) is projection of
  // z-polarization (A1) on x- and y-directions
  pTo6Space[0] += sa * (G2D * A2);
  pTo6Space[1] += -ca * (G2D * A2);
  pTo6Space[2] += G2D * A1;
  // B-field -> negative derivative wrt polarization shift of E-field
  pTo6Space[3] += -sa * (G2D * A1);
  pTo6Space[4] += ca * (G2D * A1);
  pTo6Space[5] += G2D * A2;
}

/** Gauss3D construction with */
Gauss3D::Gauss3D(vector<sunrealtype> dis_, vector<sunrealtype> axis_,
                 sunrealtype Amp_,
                 // vector<sunrealtype> pol_,
                 sunrealtype phip_, sunrealtype w0_, sunrealtype zr_,
                 sunrealtype Ph0_, sunrealtype PhA_) {
  dis = dis_;   /** - center it approaches */
  axis = axis_; /** - direction from where it comes */
  Amp = Amp_;   /** - amplitude */
  // pol=pol_;
  phip = phip_; /** - polarization rotation form TE-mode */
  w0 = w0_;     /** - taille */
  zr = zr_;     /** - Rayleigh length */
  Ph0 = Ph0_;   /** - beam center */
  PhA = PhA_;   /** - beam length */
  lambda = numbers::pi * w0 * w0 / zr;
  A1 = Amp * cos(phip);
  A2 = Amp * sin(phip);
}

/** Gauss3D implementation in space */
void Gauss3D::addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                         sunrealtype *pTo6Space) const {
  x -= dis[0];
  y -= dis[1];
  z -= dis[2];
  const sunrealtype zg = x * axis[0] + y * axis[1] + z * axis[2];
  const sunrealtype r = sqrt((x * x + y * y + z * z) - zg * zg);
  const sunrealtype wz = w0 * sqrt(1 + (zg * zg / zr / zr));
  const sunrealtype gz = atan(zg / zr);
  sunrealtype Rz = NAN;
  if (zg != 0)
    Rz = zg * (1 + (zr * zr / zg / zg));
  else
    Rz = 1e308;
  const sunrealtype k = 2 * numbers::pi / lambda;
  const sunrealtype PhF = -k * r * r / (2 * Rz) + gz - k * zg;
  const sunrealtype G3D = (w0 / wz) * exp(-r * r / wz / wz) *
                    exp(-(zg - Ph0) * (zg - Ph0) / PhA / PhA) * cos(PhF);
  const sunrealtype ca = axis[0];
  const sunrealtype sa = sqrt(1 - ca * ca);
  pTo6Space[0] += sa * (G3D * A2);
  pTo6Space[1] += -ca * (G3D * A2);
  pTo6Space[2] += G3D * A1;
  pTo6Space[3] += -sa * (G3D * A1);
  pTo6Space[4] += ca * (G3D * A1);
  pTo6Space[5] += G3D * A2;
}

/** Evaluate lattice point values to zero and add field values */
void ICSetter::eval(sunrealtype x, sunrealtype y, sunrealtype z,
                    sunrealtype *pTo6Space) {
  pTo6Space[0] = 0;
  pTo6Space[1] = 0;
  pTo6Space[2] = 0;
  pTo6Space[3] = 0;
  pTo6Space[4] = 0;
  pTo6Space[5] = 0;
  add(x, y, z, pTo6Space);
}

/** Add all initial field values to the lattice space */
void ICSetter::add(sunrealtype x, sunrealtype y, sunrealtype z,
                   sunrealtype *pTo6Space) {
  for (const auto wave : planeWaves1D)
    wave.addToSpace(x, y, z, pTo6Space);
  for (const auto wave : planeWaves2D)
    wave.addToSpace(x, y, z, pTo6Space);
  for (const auto wave : planeWaves3D)
    wave.addToSpace(x, y, z, pTo6Space);
  for (const auto wave : gauss1Ds)
    wave.addToSpace(x, y, z, pTo6Space);
  for (const auto wave : gauss2Ds)
    wave.addToSpace(x, y, z, pTo6Space);
  for (const auto wave : gauss3Ds)
    wave.addToSpace(x, y, z, pTo6Space);
}

/** Add plane waves in 1D to their container vector */
void ICSetter::addPlaneWave1D(vector<sunrealtype> k, vector<sunrealtype> p,
                              vector<sunrealtype> phi) {
  planeWaves1D.emplace_back(PlaneWave1D(k, p, phi));
}

/** Add plane waves in 2D to their container vector */
void ICSetter::addPlaneWave2D(vector<sunrealtype> k, vector<sunrealtype> p,
                              vector<sunrealtype> phi) {
  planeWaves2D.emplace_back(PlaneWave2D(k, p, phi));
}

/** Add plane waves in 3D to their container vector */
void ICSetter::addPlaneWave3D(vector<sunrealtype> k, vector<sunrealtype> p,
                              vector<sunrealtype> phi) {
  planeWaves3D.emplace_back(PlaneWave3D(k, p, phi));
}

/** Add Gaussian waves in 1D to their container vector */
void ICSetter::addGauss1D(vector<sunrealtype> k, vector<sunrealtype> p,
                          vector<sunrealtype> xo, sunrealtype phig_,
                          vector<sunrealtype> phi) {
  gauss1Ds.emplace_back(Gauss1D(k, p, xo, phig_, phi));
}

/** Add Gaussian waves in 2D to their container vector */
void ICSetter::addGauss2D(vector<sunrealtype> dis_, vector<sunrealtype> axis_,
                          sunrealtype Amp_, sunrealtype phip_, sunrealtype w0_,
                          sunrealtype zr_, sunrealtype Ph0_, sunrealtype PhA_) {
  gauss2Ds.emplace_back(
      Gauss2D(dis_, axis_, Amp_, phip_, w0_, zr_, Ph0_, PhA_));
}

/** Add Gaussian waves in 3D to their container vector */
void ICSetter::addGauss3D(vector<sunrealtype> dis_, vector<sunrealtype> axis_,
                          sunrealtype Amp_, sunrealtype phip_, sunrealtype w0_,
                          sunrealtype zr_, sunrealtype Ph0_, sunrealtype PhA_) {
  gauss3Ds.emplace_back(
      Gauss3D(dis_, axis_, Amp_, phip_, w0_, zr_, Ph0_, PhA_));
}
