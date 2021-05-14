// Martin Duy Tat 5th March 2021
/**
 * FindKSKKTag is a class for extracting all the variables of a KSKK tag
 */

#ifndef FINDKSKKTAG
#define FINDKSKKTAG

//KpiStrongPhase
#include "KpiStrongPhase/FindKS.h"
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindKSKKTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKSKKTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKSKKTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Get the daughter track IDs, in the order (pi+ pi-)KS pi+ pi-
     */
    std::vector<int> GetDaughterTrackID() const;
    /**
     * Enumeration to label daughter particles in the order K+ K- KS
     */
    enum DaughterParticle {KPLUS, KMINUS, KSHORT};
    /**
     * Get the \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double GetChi2VeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double GetKSMassVeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ decay length, from fit
     */
    double GetDecayLengthFit() const;
    /**
     * Get the \f$K_S^0\f$ decay length error, from fit
     */
    double GetDecayLengthErrorFit() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetChi2Fit() const;
    /**
     * Get the \f$\pi^+\f$ daughter momentum component from the MDC track
     */
    double GetKSPiPlusP(int i) const;
    /**
     * Get the \f$\pi^-\f$ daughter momentum component from the MDC track
     */
    double GetKSPiMinusP(int i) const;
    /**
     * Get \f$K_S^0\f$ momentum component after vertex fit
     * @param i Component
     */
    double GetKShortP(int i) const;
    /**
     * Get \f$K^+\f$ momentum component
     * @param i Component
     */
    double GetKPlusP(int i) const;
    /**
     * Get \f$K^-\f$ momentum component
     * @param i Component
     */
    double GetKMinusP(int i) const;
    /**
     * Get flag of Kalman fit success
     */
    int GetKalmanFitSuccess() const;
    /**
     * Get Kalman fit \f$chi^2\f$
     */
    double GetKalmanFitChi2() const;
    /**
     * Get \f$K_S^0\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetKShortPKalmanFit(int i) const;
    /**
     * Get \f$K^+\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetKPlusPKalmanFit(int i) const;
    /**
     * Get \f$K^-\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetKMinusPKalmanFit(int i) const;
  private:
    /**
     * Daughter track IDs, in the order (pi+ pi-)KS pi+ pi-
     */
    std::vector<int> m_DaughterTrackID;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double m_DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    double m_DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    double m_DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_Chi2Fit;
    /**
     * The \f$\pi^+\f$ daughter momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiPlusP;
    /**
     * The \f$\pi^-\f$ daughter momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiMinusP;
    /**
     * \f$K_S^0\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KShortP;
    /**
     * \f$K^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KPlusP;
    /**
     * \f$K^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KMinusP;
    /**
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of tracks
     */
    int m_KalmanFitSuccess;
    /**
     * \f$\chi^2\f$ of Kalman fit
     */
    double m_KalmanFitChi2;
    /**
     * Kalman fitted \f$K_S^0\f$ four- momentum
     */
    CLHEP::HepLorentzVector m_KShortPKalmanFit;
    /**
     * Kalman fitted \f$K^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KPlusPKalmanFit;
    /**
     * Kalman fitted \f$K^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KMinusPKalmanFit;
};

#endif
