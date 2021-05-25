// Martin Duy Tat 24th May 2021
/**
 * FindKLKKTagInfo is class for finding \f$K_L^0K^+K^-\f$ tags
 * It reconstructs the kaons on the other side of the signal and looks for missing energy and mass
 */

#ifndef FINDKLKKTAGINFO
#define FINDKLKKTAGINFO
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Boss
#include "DTagTool/DTagTool.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>

class FindKLKKTagInfo {
  public: 
    /**
     * Default constructor, initializes everything to zero
     */
    FindKLKKTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKLKKTagInfo();
    /**
     * Start looking for \f$K_L\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param DTTool DTagTool object with all the tag information
     */
    StatusCode findKLKKTagInfo(DTagToolIterator DTTool_iter, DTagTool DTTool);
    /**
     * Get the ith component of the \f$K^+\f$ four-momentum vector on the other side
     */
    double GetKPlusP(int i) const;
    /**
     * Get the ith component of the \f$K^-\f$ four-momentum vector on the other side
     */
    double GetKMinusP(int i) const;
    /**
     * Get the ith component of the $K_L^0\f$ four-momentum vector, from missing momentum
     */
    double GetKLongP(int i) const;
    /**
     * Get missing mass squared
     */
    double GetMMiss2() const;
    /**
     * Get flag of Kalman fit success
     */
    int GetKalmanFitSuccess() const;
    /**
     * Get Kalman fit \f$chi^2\f$
     */
    double GetKalmanFitChi2() const;
    /**
     * Get the ith component of the \f$K^+\f$ four-momentum vector on the other side, after Kalman fit
     */
    double GetKPlusPKalmanFit(int i) const;
    /**
     * Get the ith component of the \f$K^-\f$ four-momentum vector on the other side, after Kalman fit
     */
    double GetKMinusPKalmanFit(int i) const;
    /**
     * Get the ith component of the $K_L^0\f$ four-momentum vector, from missing momentum, after Kalman fit
     */
    double GetKLongPKalmanFit(int i) const;
    /**
     * Get the energy of the nearest shower
     */
    double GetNearestShowerEnergy() const;
    /**
     * Get the (cosine) angle of the nearest shower
     */
    double GetNearestShowerCosAngle() const;
    /**
     * Get vector of track IDs of \f$K^+K^-\f$ pair
     */
    std::vector<int> GetDaughterTrackID() const;
  private:
    /**
     * The \f$\pi^+\f$ four-momentum vector on the other side
     */
    CLHEP::HepLorentzVector m_KPlusP;
    /**
     * The \f$\pi^-\f$ four-momentum vector on the other side
     */
    CLHEP::HepLorentzVector m_KMinusP;
    /**
     * The $K_L^0\f$ four-momentum vector, from missing momentum
     */
    CLHEP::HepLorentzVector m_KLongP;
    /**
     * Flag of Kalman fit success
     */
    int m_KalmanFitSuccess;
    /**
     * Kalman fit \f$chi^2\f$
     */
    double m_KalmanFitChi2;
    /**
     * The \f$\pi^+\f$ four-momentum vector on the other side, after Kalman fit
     */
    CLHEP::HepLorentzVector m_KPlusPKalmanFit;
    /**
     * The \f$\pi^-\f$ four-momentum vector on the other side, after Kalman fit
     */
    CLHEP::HepLorentzVector m_KMinusPKalmanFit;
    /**
     * The $K_L^0\f$ four-momentum vector, from missing momentum, after Kalman fit
     */
    CLHEP::HepLorentzVector m_KLongPKalmanFit;
    /**
     * Shower energy of shower nearest \f$K_L^0\f$
     */
    double m_NearestShowerEnergy;
    /**
     * Cosine of angle between \f$K_L^0\f$ and shower
     */
    double m_NearestShowerCosAngle;
    /**
     * Vector of \f$K^+K^-\f$ daughter track IDs, in the order \f$K^+\f$ \f$K^-\f$
     */
    std::vector<int> m_DaughterTrackID;
};

#endif
