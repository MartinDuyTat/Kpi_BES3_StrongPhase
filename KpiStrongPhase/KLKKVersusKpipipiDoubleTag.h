// Martin Duy Tat 25th May 2021
/**
 * KLKKVersusKpipipiDoubleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a double \f$D\to K_L^0K^+K^-\f$ vs \f$D^0\to K\pi\pi\pi\f$ tag
 * A kinematic fit of $K_L^0K^+K^-$ is done to get more accurate momenta for binning
 * The tag side is constructed first because the signal mode has missing momentum
 */

#ifndef KLKKPIPIVERSUSKPIPIPIDOUBLETAG
#define KLKKPIPIVERSUSKPIPIPIDOUBLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KLKKVersusKpipipiDoubleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KLKKVersusKpipipiDoubleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KLKKVersusKpipipiDoubleTag();
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
    StatusCode FillTuple(DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool);
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
    NTuple::Item<double> m_TagDMass;
    /**
     * Beam constrained mass of signal
     */
    NTuple::Item<double> m_TagMBC;
    /**
     * \f$E_D - E_{\text{beam}}\f$ of signal
     */
    NTuple::Item<double> m_TagDeltaE;
    /**
     * Beam energy of signal
     */
    NTuple::Item<double> m_TagBeamE;
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
     * Signal \f$K_L^0\f$ momentum along \f$x\f$ from missing momentum
     */
    NTuple::Item<double> m_SignalKLpx;
    /**
     * Signal \f$K_L^0\f$ momentum along \f$y\f$ from missing momentum
     */
    NTuple::Item<double> m_SignalKLpy;
    /**
     * Signal \f$K_L^0\f$ momentum along \f$z\f$ from missing momentum
     */
    NTuple::Item<double> m_SignalKLpz;
    /**
     * Signal \f$K_L^0\f$ energy from missing momentum
     */
    NTuple::Item<double> m_SignalKLenergy;
    /**
     * Signal missing mass squared
     */
    NTuple::Item<double> m_SignalMMiss2;
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
     * Signal Kalman fitted \f$K_L^0\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKLpxKalmanFit;
    /**
     * Signal Kalman fitted \f$K_L^0\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKLpyKalmanFit;
    /**
     * Signal Kalman fitted \f$K_L^0\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKLpzKalmanFit;
    /**
     * Signal Kalman fitted \f$K_L^0\f$ energy
     */
    NTuple::Item<double> m_SignalKLenergyKalmanFit;
    /**
     * Energy of shower closest to missing momentum
     */
    NTuple::Item<double> m_SignalNearestShowerEnergy;
    /**
     * Angular separation of shower closest to missing momentum (cosine)
     */
    NTuple::Item<double> m_SignalNearestShowerCosAngle;
    /**
     * Equal to 1 if the signal daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_SignalIsSameDMother;
    /**
     * Equal to 1 if the signal daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_SignalPIDTrue;
    /**
     * The signal \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_SignalKPlusTrueID;
    /**
     * The signal \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_SignalKMinusTrueID;
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
     * Tag \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagPi1px;
    /**
     * Tag \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagPi1py;
    /**
     * Tag \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagPi1pz;
    /**
     * Tag \f$\pi\f$ energy
     */
    NTuple::Item<double> m_TagPi1energy;
    /**
     * Tag \f$\pi\f$ charge
     */
    NTuple::Item<int> m_TagPi1Charge;
    /**
     * Tag \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagPi2px;
    /**
     * Tag \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagPi2py;
    /**
     * Tag \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagPi2pz;
    /**
     * Tag \f$\pi\f$ energy
     */
    NTuple::Item<double> m_TagPi2energy;
    /**
     * Tag \f$\pi\f$ charge
     */
    NTuple::Item<int> m_TagPi2Charge;
    /**
     * Tag \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagPi3px;
    /**
     * Tag \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagPi3py;
    /**
     * Tag \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagPi3pz;
    /**
     * Tag \f$\pi\f$ energy
     */
    NTuple::Item<double> m_TagPi3energy;
    /**
     * Tag \f$\pi\f$ charge
     */
    NTuple::Item<int> m_TagPi3Charge;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tag tracks
     */
    NTuple::Item<int> m_Tag12KSFitSuccess;
    /**
     * The tag \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_Tag12DecayLengthVeeVertex;
    /**
     * The tag \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_Tag12Chi2VeeVertex;
    /**
     * The tag \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_Tag12KSMassVeeVertex;
    /**
     * The tag \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_Tag12DecayLengthFit;
    /**
     * The tag \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_Tag12DecayLengthErrorFit;
    /**
     * The tag \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_Tag12Chi2Fit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tag tracks
     */
    NTuple::Item<int> m_Tag13KSFitSuccess;
    /**
     * The tag \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_Tag13DecayLengthVeeVertex;
    /**
     * The tag \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_Tag13Chi2VeeVertex;
    /**
     * The tag \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_Tag13KSMassVeeVertex;
    /**
     * The tag \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_Tag13DecayLengthFit;
    /**
     * The tag \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_Tag13DecayLengthErrorFit;
    /**
     * The tag \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_Tag13Chi2Fit;
    /**
     * Equal to 1 if the tag daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_TagIsSameDMother;
    /**
     * Equal to 1 if the tag daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_TagPIDTrue;
    /**
     * The tag \f$K^+\f$ true PID
     */
    NTuple::Item<int> m_TagKTrueID;
    /**
     * The tag \f$K^-\f$ true PID
     */
    NTuple::Item<int> m_TagPi1TrueID;
    /**
     * The tag \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_TagPi2TrueID;
    /**
     * The tag \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_TagPi3TrueID;
};

#endif
