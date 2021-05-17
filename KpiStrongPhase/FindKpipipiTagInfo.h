// Martin Duy Tat 25th March 2021
/**
 * FindKpipipiTag is a class for extracting all the variables of a Kpipipi tag
 */

#ifndef FINDKPIPIPITAG
#define FINDKPIPIPITAG

// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindKpipipiTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKpipipiTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKpipipiTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Enumeration to label daughter particles in the order K+ K- pi+ pi-
     */
    enum DaughterParticle {KAON, PION1, PION2, PION3};
    /**
     * Get the daughter track IDs
     */
    std::vector<int> GetDaughterTrackID() const;
    /**
     * Get \f$K\f$ momentum component
     * @param i Component
     */
    double GetKP(int i) const;
    /**
     * Get \f$\pi\f$ momentum component
     * @param i Component
     */
    double GetPi1P(int i) const;
    /**
     * Get \f$\pi\f$ momentum component
     * @param i Component
     */
    double GetPi2P(int i) const;
    /**
     * Get \f$\pi\f$ momentum component
     * @param i Component
     */
    double GetPi3P(int i) const;
    /**
     * Get \f$K\f$ charge
     */
    int GetKCharge() const;
    /**
     * Get \f$\pi\f$ charge
     */
    int GetPi1Charge() const;
    /**
     * Get \f$\pi\f$ charge
     */
    int GetPi2Charge() const;
    /**
     * Get \f$\pi\f$ charge
     */
    int GetPi3Charge() const;
    /**
     * Get flag of \f$K_S^0\f$ fit success of tracks
     */
    int GetKSFitSuccess12() const;
    /**
     * Get the \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex12() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double GetChi2VeeVertex12() const;
    /**
     * Get the \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double GetKSMassVeeVertex12() const;
    /**
     * Get the \f$K_S^0\f$ decay length, from fit
     */
    double GetDecayLengthFit12() const;
    /**
     * Get the \f$K_S^0\f$ decay length error, from fit
     */
    double GetDecayLengthErrorFit12() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetChi2Fit12() const;
    /**
     * Get flag of \f$K_S^0\f$ fit success of tracks
     */
    int GetKSFitSuccess13() const;
    /**
     * Get the \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex13() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double GetChi2VeeVertex13() const;
    /**
     * Get the \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double GetKSMassVeeVertex13() const;
    /**
     * Get the \f$K_S^0\f$ decay length, from fit
     */
    double GetDecayLengthFit13() const;
    /**
     * Get the \f$K_S^0\f$ decay length error, from fit
     */
    double GetDecayLengthErrorFit13() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetChi2Fit13() const;
  private:
    /**
     * Daughter track IDs
     */
    std::vector<int> m_DaughterTrackID;
    /**
     * \f$K\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KP;
    /**
     * \f$\pi\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_Pi1P;
    /**
     * \f$\pi\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_Pi2P;
    /**
     * \f$\pi\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_Pi3P;
    /**
     * \f$K\f$ charge
     */
    int m_KCharge;
    /**
     * \f$\pi\f$ charge
     */
    int m_Pi1Charge;
    /**
     * \f$\pi\f$ charge
     */
    int m_Pi2Charge;
    /**
     * \f$\pi\f$ charge
     */
    int m_Pi3Charge;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    int m_12KSFitSuccess;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double m_12DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_12Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_12KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    double m_12DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    double m_12DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_12Chi2Fit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    int m_13KSFitSuccess;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double m_13DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_13Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_13KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    double m_13DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    double m_13DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_13Chi2Fit;
};

#endif
