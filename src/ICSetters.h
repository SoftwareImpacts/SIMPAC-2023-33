/////////////////////////////////////////////////////////////////////////////////////
/// @file ICSetters.h
/// @brief Declaration of the plane wave and Gaussian wave packets in 1D, 2D, 3D
/////////////////////////////////////////////////////////////////////////////////////

// Include Guard
#ifndef ICSETTERS
#define ICSETTERS

// math, constants, vector, and array
#include <cmath>
//#include <mathimf.h>
#include <numbers>
#include <array>
#include <vector>

#include <sundials/sundials_types.h> /* definition of type sunrealtype */

using namespace std;

/** @brief super-class for plane waves
 *
 * They are given in the form \f$ \vec{E} = \vec{E}_0 \ \cos \left( \vec{k}
 * \cdot \vec{x} - \vec{\phi} \right) \f$ */
class PlaneWave {
protected:
  /// wavenumber \f$ k_x \f$
  sunrealtype kx;
  /// wavenumber \f$ k_y \f$
  sunrealtype ky;
  /// wavenumber \f$ k_z \f$
  sunrealtype kz;
  /// polarization & amplitude in x-direction, \f$ p_x \f$
  sunrealtype px;
  /// polarization & amplitude in y-direction, \f$ p_y \f$
  sunrealtype py;
  /// polarization & amplitude in z-direction, \f$ p_z \f$
  sunrealtype pz;
  /// phase shift in x-direction, \f$ \phi_x \f$
  sunrealtype phix;
  /// phase shift in y-direction, \f$ \phi_y \f$
  sunrealtype phiy;
  /// phase shift in z-direction, \f$ \phi_z \f$
  sunrealtype phiz;
};

/** @brief class for plane waves in 1D */
class PlaneWave1D : public PlaneWave {
public:
  /// construction with default parameters
  PlaneWave1D(vector<sunrealtype> k = {1, 0, 0},
              vector<sunrealtype> p = {0, 0, 1},
              vector<sunrealtype> phi = {0, 0, 0});
  /// function for the actual implementation in the lattice
  void addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                  sunrealtype *pTo6Space) const;
};

/** @brief class for plane waves in 2D */
class PlaneWave2D : public PlaneWave {
public:
  /// construction with default parameters
  PlaneWave2D(vector<sunrealtype> k = {1, 0, 0},
              vector<sunrealtype> p = {0, 0, 1},
              vector<sunrealtype> phi = {0, 0, 0});
  /// function for the actual implementation in the lattice
  void addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                  sunrealtype *pTo6Space) const;
};

/** @brief class for plane waves in 3D */
class PlaneWave3D : public PlaneWave {
public:
  /// construction with default parameters
  PlaneWave3D(vector<sunrealtype> k = {1, 0, 0},
              vector<sunrealtype> p = {0, 0, 1},
              vector<sunrealtype> phi = {0, 0, 0});
  /// function for the actual implementation in space
  void addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                  sunrealtype *pTo6Space) const;
};

/** @brief class for Gaussian waves in 1D
 *
 * They are given in the form \f$ \vec{E}=\vec{p} \, \exp \left(
 * -(\vec{x}-\vec{x}_0)^2 / \Phi_g^2 \right) \, \cos(\vec{k} \cdot \vec{x}) \f$
 */
class Gauss1D {
private:
  /// wavenumber \f$ k_x \f$
  sunrealtype kx;
  /// wavenumber \f$ k_y \f$
  sunrealtype ky;
  /// wavenumber \f$ k_z \f$
  sunrealtype kz;
  /// polarization & amplitude in x-direction, \f$ p_x \f$
  sunrealtype px;
  /// polarization & amplitude in y-direction, \f$ p_y \f$
  sunrealtype py;
  /// polarization & amplitude in z-direction, \f$ p_z \f$
  sunrealtype pz;
  /// phase shift in x-direction, \f$ \phi_x \f$
  sunrealtype phix;
  /// phase shift in y-direction, \f$ \phi_y \f$
  sunrealtype phiy;
  /// phase shift in z-direction, \f$ \phi_z \f$
  sunrealtype phiz;
  /// center of pulse in x-direction, \f$ x_0 \f$
  sunrealtype x0x;
  /// center of pulse in y-direction, \f$ y_0 \f$
  sunrealtype x0y;
  /// center of pulse in z-direction, \f$ z_0 \f$
  sunrealtype x0z;
  /// pulse width \f$ \Phi_g \f$
  sunrealtype phig;

public:
  /// construction with default parameters
  Gauss1D(vector<sunrealtype> k = {1, 0, 0}, vector<sunrealtype> p = {0, 0, 1},
          vector<sunrealtype> xo = {0, 0, 0}, sunrealtype phig_ = 1.0l,
          vector<sunrealtype> phi = {0, 0, 0});
  /// function for the actual implementation in space
  void addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                  sunrealtype *pTo6Space) const;

public:
};

/** @brief class for Gaussian waves in 2D
 *
 * They are given in the form
 * \f$ \vec{E}= A \, \vec{\epsilon} \ \sqrt{\frac{\omega_0}{\omega(z)}} \, \exp
 * \left(-r/\omega(z) \right)^2 \, \exp \left(-((z_g-\Phi_0)/\Phi_A)^2 \right)
 * \, \cos \left( \frac{k \, r^2}{2R(z)} + g(z) - k\, z_g \right)  \f$ with
 * - propagation direction (subtracted distance to origin) \f$ z_g \f$
 * - radial distance to propagation axis \f$ r = \sqrt{\vec{x}^2 -z_g^2} \f$
 * - \f$ k = 2\pi / \lambda \f$
 * - waist at position z, \f$ \omega(z) = w_0 \, \sqrt{1+(z_g/z_R)^2} \f$
 * - Gouy phase \f$ g(z) = \tan^{-1}(z_g/z_r) \f$
 * - beam curvature \f$ R(z) = z_g \, (1+(z_r/z_g)^2) \f$
 * obtained via the chosen parameters */
class Gauss2D {
private:
  /// distance maximum to origin
  vector<sunrealtype> dis;
  /// normalized propagation axis
  vector<sunrealtype> axis;
  /// amplitude \f$ A\f$
  sunrealtype Amp;
  /// polarization rotation from TE-mode around propagation direction
  // that determines \f$ \vec{\epsilon}\f$ above
  sunrealtype phip;
  /// taille \f$ \omega_0 \f$
  sunrealtype w0;
  /// Rayleigh length \f$ z_R = \pi \omega_0^2 / \lambda \f$
  sunrealtype zr;
  /// center of beam \f$ \Phi_0 \f$
  sunrealtype Ph0;
  /// length of beam \f$ \Phi_A \f$
  sunrealtype PhA;
  /// amplitude projection on TE-mode
  sunrealtype A1;
  /// amplitude projection on xy-plane
  sunrealtype A2;
  /// wavelength \f$ \lambda \f$
  sunrealtype lambda;

public:
  /// construction with default parameters
  Gauss2D(vector<sunrealtype> dis_ = {0, 0, 0},
          vector<sunrealtype> axis_ = {1, 0, 0}, sunrealtype Amp_ = 1.0l,
          sunrealtype phip_ = 0, sunrealtype w0_ = 1e-5, sunrealtype zr_ = 4e-5,
          sunrealtype Ph0_ = 2e-5, sunrealtype PhA_ = 0.45e-5);
  /// function for the actual implementation in space
  void addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                  sunrealtype *pTo6Space) const;

public:
};

/** @brief class for Gaussian waves in 3D
 *
 * They are given in the form
 * \f$ \vec{E}= A \, \vec{\epsilon} \ \frac{\omega_0}{\omega(z)} \, \exp
 * \left(-r/\omega(z) \right)^2 \, \exp \left(-((z_g-\Phi_0)/\Phi_A)^2 \right)
 * \, \cos \left( \frac{k \, r^2}{2R(z)} + g(z) - k\, z_g \right)  \f$ with
 * - propagation direction (subtracted distance to origin) \f$ z_g \f$
 * - radial distance to propagation axis \f$ r = \sqrt{\vec{x}^2 -z_g^2} \f$
 * - \f$ k = 2\pi / \lambda \f$
 * - waist at position z, \f$ \omega(z) = w_0 \, \sqrt{1+(z_g/z_R)^2} \f$
 * - Gouy phase \f$ g(z) = \tan^{-1}(z_g/z_r) \f$
 * - beam curvature \f$ R(z) = z_g \, (1+(z_r/z_g)^2) \f$
 * obtained via the chosen parameters */
class Gauss3D {
private:
  /// distance maximum to origin
  vector<sunrealtype> dis;
  /// normalized propagation axis
  vector<sunrealtype> axis;
  /// amplitude \f$ A\f$
  sunrealtype Amp;
  /// polarization rotation from TE-mode around propagation direction
  // that determines \f$ \vec{\epsilon}\f$ above
  sunrealtype phip;
  // polarization
  // vector<sunrealtype> pol;
  /// taille \f$ \omega_0 \f$
  sunrealtype w0;
  /// Rayleigh length \f$ z_R = \pi \omega_0^2 / \lambda \f$
  sunrealtype zr;
  /// center of beam \f$ \Phi_0 \f$
  sunrealtype Ph0;
  /// length of beam \f$ \Phi_A \f$
  sunrealtype PhA;
  /// amplitude projection on TE-mode (z-axis)
  sunrealtype A1;
  /// amplitude projection on xy-plane
  sunrealtype A2;
  /// wavelength \f$ \lambda \f$
  sunrealtype lambda;

public:
  /// construction with default parameters
  Gauss3D(vector<sunrealtype> dis_ = {0, 0, 0},
          vector<sunrealtype> axis_ = {1, 0, 0}, sunrealtype Amp_ = 1.0l,
          sunrealtype phip_ = 0,
          // sunrealtype pol_={0,0,1},
          sunrealtype w0_ = 1e-5, sunrealtype zr_ = 4e-5,
          sunrealtype Ph0_ = 2e-5, sunrealtype PhA_ = 0.45e-5);
  /// function for the actual implementation in space
  void addToSpace(sunrealtype x, sunrealtype y, sunrealtype z,
                  sunrealtype *pTo6Space) const;

public:
};

/** @brief ICSetter class to initialize wave types with default parameters */
class ICSetter {
private:
  /// container vector for plane waves in 1D
  vector<PlaneWave1D> planeWaves1D;
  /// container vector for plane waves in 2D
  vector<PlaneWave2D> planeWaves2D;
  /// container vector for plane waves in 3D
  vector<PlaneWave3D> planeWaves3D;
  /// container vector for Gaussian waves in 1D
  vector<Gauss1D> gauss1Ds;
  /// container vector for Gaussian waves in 2D
  vector<Gauss2D> gauss2Ds;
  /// container vector for Gaussian waves in 3D
  vector<Gauss3D> gauss3Ds;

public:
  /// function to set all coordinates to zero and then `add` the field values
  void eval(sunrealtype x, sunrealtype y, sunrealtype z,
            sunrealtype *pTo6Space);
  /// function to fill the lattice space with initial field values
  // of all field vector containers
  void add(sunrealtype x, sunrealtype y, sunrealtype z, sunrealtype *pTo6Space);
  /// function to add plane waves in 1D to their container vector
  void addPlaneWave1D(vector<sunrealtype> k = {1, 0, 0},
                      vector<sunrealtype> p = {0, 0, 1},
                      vector<sunrealtype> phi = {0, 0, 0});
  /// function to add plane waves in 2D to their container vector
  void addPlaneWave2D(vector<sunrealtype> k = {1, 0, 0},
                      vector<sunrealtype> p = {0, 0, 1},
                      vector<sunrealtype> phi = {0, 0, 0});
  /// function to add plane waves in 3D to their container vector
  void addPlaneWave3D(vector<sunrealtype> k = {1, 0, 0},
                      vector<sunrealtype> p = {0, 0, 1},
                      vector<sunrealtype> phi = {0, 0, 0});
  /// function to add Gaussian waves in 1D to their container vector
  void addGauss1D(vector<sunrealtype> k = {1, 0, 0},
                  vector<sunrealtype> p = {0, 0, 1},
                  vector<sunrealtype> xo = {0, 0, 0}, sunrealtype phig_ = 1.0l,
                  vector<sunrealtype> phi = {0, 0, 0});
  /// function to add Gaussian waves in 2D to their container vector
  void addGauss2D(vector<sunrealtype> dis_ = {0, 0, 0},
                  vector<sunrealtype> axis_ = {1, 0, 0},
                  sunrealtype Amp_ = 1.0l, sunrealtype phip_ = 0,
                  sunrealtype w0_ = 1e-5, sunrealtype zr_ = 4e-5,
                  sunrealtype Ph0_ = 2e-5, sunrealtype PhA_ = 0.45e-5);
  /// function to add Gaussian waves in 3D to their container vector
  void addGauss3D(vector<sunrealtype> dis_ = {0, 0, 0},
                  vector<sunrealtype> axis_ = {1, 0, 0},
                  sunrealtype Amp_ = 1.0l, sunrealtype phip_ = 0,
                  sunrealtype w0_ = 1e-5, sunrealtype zr_ = 4e-5,
                  sunrealtype Ph0_ = 2e-5, sunrealtype PhA_ = 0.45e-5);
};

// End of Includeguard
#endif
