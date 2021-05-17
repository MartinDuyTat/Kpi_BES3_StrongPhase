// Martin Duy Tat 5th March 2021
/**
 * KSKKVersusKpiDoubleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a double \f$D\to K_S^0K^+K^-\pi^+\pi^-\f$ vs \f$D^0\to K\pi\f$ tag
 * It also runs a fit for the decay \f$K_S^0\to\pi^+\pi^-\f$ by refitting the primary and secondary vertex in the class FindKS, from this the flight significance is used to eliminate peaking background
 * A kinematic fit of $K_S^0K^+K^-$ is done to get more accurate momenta for binning
 */

#ifndef KSKKPIPIVERSUSKPIDOUBLETAG
#define KSKKPIPIVERSUSKPIDOUBLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KSKKVersusKpiDoubleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KSKKVersusKpiDoubleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KSKKVersusKpiDoubleTag();
    /**
     * This function runs when algorithm is initialized
     */
    StatusCode initialize();
    /**
     * Execution of the algorithm
     */
    StatusCode execute();
    /**
     * This function runs when algorithm is finalized
     */
    StatusCode finalize();
    /**
     * Helper function to fill in information about the tag mode
     * @param DTTool_Signal_iter Iterator pointing to signal candidate
     * @param DTTool_Signal_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the event information
     */
    StatusCode FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool);
  private:
    /**
     * Dummy variable, placeholder for more important properties to be added later
     */
    int m_dummy;
    /**
     * The NTuple itself that is filled
     */
    NTuple::Tuple *m_tuple;
    /**
     * Run number (negative for MC event)
     */
    NTuple::Item<int> m_RunNumber;
    /**
     * Event number
     */
    NTuple::Item<int> m_EventNumber;
    /**
     * Number of particles in the decay chain
     * For example, \f$D^0\to K^-\pi^+\f$ counts as \f$3\f$ particles
     */
    NTuple::Item<int> m_NumberParticles;
    /**
     * Array of particle IDs of all particles in the decay chain
     */
    NTuple::Array<int> m_pdgID;
    /**
     * Array of indices referring to the particle mother
     * For example, an index 0 means that ```m_pdgID[0]``` is the mother
     */
    NTuple::Array<int> m_MotherIndex;
    /**
     * Number of particles in the decay chain without resonances
     */
    NTuple::Item<int> m_NumberParticlesStripped;
    /**
     * Array of particle IDs of all particles in the decay chain without resonances
     */
    NTuple::Array<int> m_pdgIDStripped;
    /**
     * Array of indices referring to the particle mother without resonances
     */
    NTuple::Array<int> m_MotherIndexStripped;
    /**
     * Generator label of the decay mode
     */
    NTuple::Item<int> m_MCmode;
    /**
     * Array of true x momenta
     */
    NTuple::Array<double> m_TruePx;
    /**
     * Array of true y momenta
     */
    NTuple::Array<double> m_TruePy;
    /**
     * Array of true z momenta
     */
    NTuple::Array<double> m_TruePz;
    /**
     * Array of true energies
     */
    NTuple::Array<double> m_TrueEnergy;
    /**
     * Invariant mass of signal \f$D\f$ meson
     */
    NTuple::Item<double> m_SignalDMass;
    /**
     * Beam constrained mass of signal
     */
    NTuple::Item<double> m_SignalMBC;
    /**
     * \f$E_D - E_{\text{beam}}\f$ of signal
     */
    NTuple::Item<double> m_SignalDeltaE;
    /**
     * Beam energy of signal
     */
    NTuple::Item<double> m_SignalBeamE;
    /**
     * Invariant mass of tag \f$D\f$ meson
     */
    NTuple::Item<double> m_TagDMass;
    /**
     * Beam constrained mass of tag
     */
    NTuple::Item<double> m_TagMBC;
    /**
     * \f$E_D - E_{\text{beam}}\f$ of tag
     */
    NTuple::Item<double> m_TagDeltaE;
    /**
     * Beam energy of tag
     */
    NTuple::Item<double> m_TagBeamE;
    /**
     * Signal \f$D\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalDpx;
    /**
     * Signal \f$D\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalDpy;
    /**
     * Signal \f$D\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalDpz;
    /**
     * Signal \f$D\f$ energy
     */
    NTuple::Item<double> m_SignalDenergy;
    /**
     * Tag \f$D\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagDpx;
    /**
     * Tag \f$D\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagDpy;
    /**
     * Tag \f$D\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagDpz;
    /**
     * Tag \f$D\f$ energy
     */
    NTuple::Item<double> m_TagDenergy;
    /**
     * The signal \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_SignalDecayLengthVeeVertex;
    /**
     * The signal \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_SignalChi2VeeVertex;
    /**
     * The signal \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_SignalKSMassVeeVertex;
    /**
     * The signal \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_SignalDecayLengthFit;
    /**
     * The signal \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_SignalDecayLengthErrorFit;
    /**
     * The signal \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_SignalChi2Fit;
    /**
     * Signal \f$K_S^0\f$ momentum along \f$x\f$ after vertex fit
     */
    NTuple::Item<double> m_SignalKSpx;
    /**
     * Signal \f$K_S^0\f$ momentum along \f$y\f$ after vertex fit
     */
    NTuple::Item<double> m_SignalKSpy;
    /**
     * Signal \f$K_S^0\f$ momentum along \f$z\f$ after vertex fit
     */
    NTuple::Item<double> m_SignalKSpz;
    /**
     * Signal \f$K_S^0\f$ energy after vertex fit
     */
    NTuple::Item<double> m_SignalKSenergy;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiPluspx;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$y\f$ from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiPluspy;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$z\f$ from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiPluspz;
    /**
     * The \f$\pi^+\f$ daughter energy from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiPlusenergy;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiMinuspx;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$y\f$ from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiMinuspy;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$z\f$ from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiMinuspz;
    /**
     * The \f$\pi^-\f$ daughter energy from the MDC track
     */
    NTuple::Item<double> m_SignalKSPiMinusenergy;
    /**
     * Signal \f$K^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKPluspx;
    /**
     * Signal \f$K^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKPluspy;
    /**
     * Signal \f$K^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKPluspz;
    /**
     * Signal \f$K^+\f$ energy
     */
    NTuple::Item<double> m_SignalKPlusenergy;
    /**
     * Signal \f$K^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKMinuspx;
    /**
     * Signal \f$K^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKMinuspy;
    /**
     * Signal \f$K^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKMinuspz;
    /**
     * Signal \f$K^-\f$ energy
     */
    NTuple::Item<double> m_SignalKMinusenergy;
    /**
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of signal tracks
     */
    NTuple::Item<int> m_SignalKalmanFitSuccess;
    /**
     * Signal \f$\chi^2\f$ of Kalman fit
     */
    NTuple::Item<double> m_SignalKalmanFitChi2;
    /**
     * Signal Kalman fitted \f$K^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKPluspxKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKPluspyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKPluspzKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ energy
     */
    NTuple::Item<double> m_SignalKPlusenergyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKMinuspxKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKMinuspyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKMinuspzKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ energy
     */
    NTuple::Item<double> m_SignalKMinusenergyKalmanFit;
    /**
     * Signal Kalman fitted \f$K_S^0\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKSpxKalmanFit;
    /**
     * Signal Kalman fitted \f$K_S^0\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKSpyKalmanFit;
    /**
     * Signal Kalman fitted \f$K_S^0\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKSpzKalmanFit;
    /**
     * Signal Kalman fitted \f$K_S^0\f$ energy
     */
    NTuple::Item<double> m_SignalKSenergyKalmanFit;
    /**
     * Equal to 1 if the signal daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_SignalIsSameDMother;
    /**
     * Equal to 1 if the signal daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_SignalPIDTrue;
    /**
     * The signal \f$\pi^+\f$ from \f$K_S^0\f$ true PID
     */
    NTuple::Item<int> m_SignalKSPiPlusTrueID;
    /**
     * The signal \f$\pi^-\f$ from \f$K_S^0\f$ true PID
     */
    NTuple::Item<int> m_SignalKSPiMinusTrueID;
    /**
     * The signal \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_SignalKPlusTrueID;
    /**
     * The signal \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_SignalKMinusTrueID;
    /**
     * Signal \f$\pi^+\f$ from \f$K_S^0\f$ true mother PID
     */
    NTuple::Item<int> m_SignalKSPiPlusMotherTrueID;
    /**
     * Signal \f$\pi^-\f$ from \f$K_S^0\f$ true mother PID
     */
    NTuple::Item<int> m_SignalKSPiMinusMotherTrueID;
    /**
     * Tag \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagPipx;
    /**
     * Tag \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagPipy;
    /**
     * Tag \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagPipz;
    /**
     * Tag \f$\pi\f$ energy
     */
    NTuple::Item<double> m_TagPienergy;
    /**
     * Tag \f$\pi\f$ charge
     */
    NTuple::Item<int> m_TagPiCharge;
    /**
     * Tag \f$K\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagKpx;
    /**
     * Tag \f$K\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagKpy;
    /**
     * Tag \f$K\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagKpz;
    /**
     * Tag \f$K\f$ energy
     */
    NTuple::Item<double> m_TagKenergy;
    /**
     * Tag \f$K\f$ charge
     */
    NTuple::Item<int> m_TagKCharge;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagHighEPi0px;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagHighEPi0py;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagHighEPi0pz;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained energy
     */
    NTuple::Item<double> m_TagHighEPi0energy;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagLowEPi0px;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagLowEPi0py;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagLowEPi0pz;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained energy
     */
    NTuple::Item<double> m_TagLowEPi0energy;
    /**
     * The tag \f$\gamma\gamma\f$ invariant mass
     */
    NTuple::Item<double> m_TagMgammagamma;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedpx;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedpy;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedpz;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedenergy;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedpx;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedpy;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedpz;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedenergy;
    /**
     * Tag \f$\pi^0\f$ kinematic fit \f$\chi^2\f$
     */
    NTuple::Item<double> m_Pi0Chi2Fit;
    /**
     * Equal to 1 if the tag daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_TagIsSameDMother;
    /**
     * Equal to 1 if all tag daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_TagIsSameDMotherAll;
    /**
     * Equal to 1 if the tag daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_TagPIDTrue;
    /**
     * The tag \f$K\f$ true PID
     */
    NTuple::Item<int> m_TagKTrueID;
    /**
     * The tag \f$\pi\f$ true PID
     */
    NTuple::Item<int> m_TagPiTrueID;
    /**
     * Tag high energy photon from \f$\pi^0\f$ true PID
     */
    NTuple::Item<int> m_TagHighEPi0PhotonTrueID;
    /**
     * Tag low energy photon from \f$\pi^0\f$ true PID
     */
    NTuple::Item<int> m_TagLowEPi0PhotonTrueID;
    /**
     * Tag high energy photon from \f$\pi^0\f$ true mother PID
     */
    NTuple::Item<int> m_TagHighEPi0PhotonMotherTrueID;
    /**
     * Tag low energy photon from \f$\pi^0\f$ true mother PID
     */
    NTuple::Item<int> m_TagLowEPi0PhotonMotherTrueID;
};

#endif
