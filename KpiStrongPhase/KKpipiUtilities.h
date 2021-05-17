// Martin Duy Tat 26th March 2021
/**
 * KKpipiUtilities is a namespace with a few useful general purpose functions
 */

#ifndef KKPIPIUTILITIES
#define KKPIPIUTILITIES

// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace KKpipiUtilities {
  /**
   * Helper function to get four-momentum from shower information
   * @param Energy Energy of photon
   * @param Theta Polar angle of photon
   * @param Phi Azimuthal angle of photon
   */
  CLHEP::HepLorentzVector GetPhoton4Vector(double Energy, double Theta, double Phi);
  /**
   * Helper function to calculate the angular separation between the EMC shower and the nearest charged track
   * @param EMCPosition Position of the EMC shower
   * @param Angle Output, the angular separation between the shower and the nearest charged track
   * @param Theta Output, the polar angle separation between the shower and the nearest charged track
   * @param Phi Output, the azimuthal angle separation between the shower and the nearest charged track
   * @return Returns true if the calculation was successful
   */
  bool GetPhotonAngularSeparation(const CLHEP::Hep3Vector &EMCPosition, double &Angle, double &Theta, double &Phi);
  /**
   * Function that calculates the missing momentum
   * The \f$D\f$ momentum is calculated more accurately using the beam energy and the \f$D\f$ mass
   * @param P_D Four-momentum of the reconstructed \f$D\f$ meson
   * @param P_X Four-momentum of all reconstructed particles on the other side
   * @param BeamE Beam energy
   * @return The missing momentum
   */
  CLHEP::HepLorentzVector GetMissingMomentum(CLHEP::HepLorentzVector P_D, CLHEP::HepLorentzVector P_X, double BeamE);
}

#endif
