// Martin Duy Tat 12th February 2021, based on code by Yu Zhang
/**
 * FindPi0Eta is class for finding \f$\pi^0\f$ or \f$\eta\f$ candidates and returning the unconstrained and kinematically constrained four-momenta of the photons
 * By default one \f$\pi^0\f$ is assumed, unless otherwise is specified in the constructor and the getters
 */

#ifndef FINDPI0ETA
#define FINDPI0ETA
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include <vector>
#include <string>

class FindPi0Eta {
  public: 
    /**
     * Default constructor, initializes all properties to zero and the number of \f$\pi^0\f$ or \f$\eta\f$ to one
     * @param npi0eta Number of \f$\pi^0\f$ or \f$\eta\f$ to look for
     * @param Particle "pi0" for \f$\pi^0\f$ and "eta" for \f$\eta\f$
     */
    FindPi0Eta(int npi0eta = 1, std::string Particle = "pi0");
    /**
     * Trivial destructor
     */
    ~FindPi0Eta();
    /**
     * Start looking for \f$\pi^0\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param PiTrackIndex List of length 2 with track indices to the two pions in the event
     */
    StatusCode findPi0Eta(DTagToolIterator &DTTool_iter, DTagTool &DTTool);
    /** 
     * Get the high energy photon unconstrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetHighEPhotonP(int i, int pi0eta_index = 0) const;
    /** 
     * Get the high energy photon unconstrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetLowEPhotonP(int i, int pi0eta_index = 0) const;
    /**
     * Get the \f$\gamma\gamma\f$ invariant mass
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetMgammagamma(int pi0eta_index = 0) const;
    /** 
     * Get the high energy photon kinematically constrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetHighEPhotonPConstrained(int i, int pi0eta_index = 0) const;
    /** 
     * Get the low energy photon kinematically constrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetLowEPhotonPConstrained(int i, int pi0eta_index = 0) const;
    /**
     * Get the \f$\chi^2\f$ of the kinematic fit
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetChi2Fit(int pi0eta_index = 0) const;
    /**
     * Get the track ID of the high energy photon
     * @param pi0eta_index 0 for first \f$pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    int GetHighEPhotonTrackID(int pi0eta_index = 0) const;
    /**
     * Get the track ID of the low energy photon
     * @param pi0eta_index 0 for first \f$pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    int GetLowEPhotonTrackID(int pi0eta_index = 0) const;
  private:
    /**
     * The high energy photon unconstrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_HighEPhotonP;
    /**
     * The low energy photon unconstrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_LowEPhotonP;
    /**
     * The high energy photon kinematically constrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_HighEPhotonPConstrained;
    /**
     * The high energy photon kinematically constrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_LowEPhotonPConstrained;
    /**
     * Kinematic fit \f$\chi^2\f$
     */
    std::vector<double> m_Chi2Fit;
    /**
     * Track ID of high energy photon
     */
    std::vector<int> m_HighEPhotonTrackID;
    /**
     * Track ID of low energy photon
     */
    std::vector<int> m_LowEPhotonTrackID;
    /**
     * Number of \f$\pi^0\f$ to look for
     */
    int m_npi0eta;
    /**
     * Particle to look for, "pi0" for \f$\pi^0\f$ and "eta" for \f$eta\f$
     */
    std::string m_Particle;
};

#endif
