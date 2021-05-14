// Martin Duy Tat 5th February 2021, based on code by Yu Zhang
/**
 * FindKS is class for finding \f$K_S^0\f$ candidates and returning the flight distance, the flight distance error, the fit \f$\chi^2\f$ and the \f$K_S^0\f$ mass from VeeVertexAlg, plus the daughter four-momenta
 * In addition, the secondary vertex is fitted again using the \f$\pi^+\pi^-\f$ tracks and a primary vertex is fitted using the vertex parameters from the secondary vertex fit
 */

#ifndef FINDKS
#define FINDKS
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "VertexFit/WTrackParameter.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>

class FindKS {
  public: 
    /**
     * Default constructor, initializes everything to zero
     * @param KSTag Set to true if we're looking for an actual \f$K_S^0\f$, set to false if we're just checking whether or not a pair of $\f$\pi^+\pi^-\f$ could be a \f$K_S^0\f$ that decayed
     * @param VetoKSIDs Vector of \f$K_S^0\f$ IDs (from DTagTool) that should be skipped
     */
    FindKS(bool KSTag, const std::vector<int> &VetoKSIDs = std::vector<int>());
    /**
     * Trivial destructor
     */
    ~FindKS();
    /**
     * Start looking for \f$K_S\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param DTTool DTagTool object with all the tag information
     * @param PiTrackIndex List of length 2 with track indices to the two pions in the event, please use this when checking if a \f$\pi^+\pi^-\f$ pair could be a \f$K_S\f$ in disguise
     */
    StatusCode findKS(DTagToolIterator DTTool_iter, DTagTool DTTool, std::vector<int> PiTrackIDs = std::vector<int>());
    /** 
     * Get decay length from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex() const;
    /** 
     * Get \f$\chi^2\f$ from VeeVertexAlg
     */
    double GetChi2VeeVertex() const;
    /** 
     * Get \f$K_S^0\f$ mass from VeeVertexAlg
     */
    double GetKSMassVeeVertex() const;
    /** 
     * Get decay length from fit
     */
    double GetDecayLengthFit() const;
    /** 
     * Get decay length error from fit
     */
    double GetDecayLengthErrorFit() const;
    /** 
     * Get \f$\chi^2\f$ from fit of primary vertex
     */
    double GetChi2Fit() const;
    /**
     * Get the \f$\pi^+\f$ daughter four-momentum from the MDC track
     * @param i Momentum component
     */
    double GetKSPiPlusP(int i) const;
    /**
     * Get the \f$\pi^-\f$ daughter four-momentum from the MDC track
     * @param i Momentum component
     */
    double GetKSPiMinusP(int i) const;
    /**
     * Get the daughter track IDs
     */
    std::vector<int> GetDaughterTrackIDs() const;
    /**
     * Get the \f$K_S^0\f$ momentum after vertex fit
     */
    CLHEP::HepLorentzVector GetKShortPFit() const;
    /**
     * Get the \f$K_S^0\f$ WTrackParameter
     */
    WTrackParameter GetWTrackParameter() const;
  private:
    /**
     * The decay length, from VeeVertexAlg
     */
    double m_DecayLengthVeeVertex;
    /**
     * The \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_KSMassVeeVertex;
    /**
     * The decay length, from fit
     */
    double m_DecayLengthFit;
    /**
     * The decay length error, from fit
     */
    double m_DecayLengthErrorFit;
    /**
     * The \f$\chi^2\f$, from fit of primary vertex
     */
    double m_Chi2Fit;
    /**
     * The \f$\pi^+\f$ daughter four-momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiPlusP;
    /**
     * The \f$\pi^-\f$ daughter four-momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiMinusP;
    /**
     * The \f$K_S^0\f$ momentum after vertex fit
     */
    CLHEP::HepLorentzVector m_KShortPFit;
    /**
     * If true, we are looking for an actual \f$K_S^0\f$ and its ID will be checked against the DTagTool::ksId, otherwise we are simply checking if a pair of \f$\pi^+\pi^-\f$ are actually \f$K_S^0\f$
     */
    bool m_KSTag;
    /**
     * List of \f$K_S^0\f$ IDs that are vetoed, if any of the \f$K_S\f$ are on this list, skip the event
     * Use this if there are other \f$K_S^0\f$ in the event that have already been identified
     */
    std::vector<int> m_VetoKSIDs;
    /**
     * List of the track IDs of the daughter \f$\pi\f$ tracks
     */
    std::vector<int> m_DaughterTrackIDs;
    /**
     * WTrackParameter of \f$K_S^0\f$
     */
    WTrackParameter m_WTrackParameter;
};

#endif
