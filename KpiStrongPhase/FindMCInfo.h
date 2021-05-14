// Martin Duy Tat 11th February 2021, based on code by Yu Zhang
/**
 * FindMCInfo is a class for finding the MC truth info of an event
 */

#ifndef FINDMCINFO
#define FINDMCINFO

// Gaudi
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"
// Boss
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>

class FindMCInfo {
  public: 
    /**
     * Default constructor which initializes all variables to zero
     */
    FindMCInfo();
    /**
     * Trivial destructor
     */
    ~FindMCInfo();
    /**
     * Helper function to calculate the MC truth information
     * @param KeepID ID of a particular resonance to keep in the MC truth list
     */
    StatusCode CalculateMCInfo(SmartDataPtr<Event::McParticleCol> MCParticleCol, IMcDecayModeSvc *IMcDecayModeService, int KeepID = 0);
    /**
     * Get the number of particles
     */
    int GetNumberParticles();
    /**
     * Get PDG ID
     * @param i Vector component
     */
    int GetpdgID(int i);
    /**
     * Get mother index
     * @param i Vector component
     */
    int GetMotherIndex(int i);
    /**
     * Get the number of particles without resonances
     */
    int GetNumberParticlesStripped();
    /**
     * Get PDG ID without resonances
     * @param i Vector component
     */
    int GetpdgIDStripped(int i);
    /**
     * Get mother index without resonances
     * @param i Vector component
     */
    int GetMotherIndexStripped(int i);
    /**
     * Get MC mode label
     */
    int GetMCmode();
    /**
     * Get true momentum in the \f$x\f$ direction
     * @param i Vector component
     */
    double GetTruePx(int i);
    /**
     * Get true momentum in the \f$y\f$ direction
     * @param i Vector component
     */
    double GetTruePy(int i);
    /**
     * Get true momentum in the \f$y\f$ direction
     * @param i Vector component
     */
    double GetTruePz(int i);
    /**
     * Get true energy
     * @param i Vector component
     */
    double GetTrueEnergy(int i);
  private:
    /**
     * Number of particles in the decay chain
     * For example, \f$D^0\to K^-\pi^+\f$ counts as \f$3\f$ particles
     */
    int m_NumberParticles;
    /**
     * Array of particle IDs of all particles in the decay chain
     */
    std::vector<int> m_pdgID;
    /**
     * Vector of indices referring to the particle mother
     * For example, an index 0 means that ```m_pdgID[0]``` is the mother
     */
    std::vector<int> m_MotherIndex;
    /**
     * Number of particles in the decay chain without resonances
     */
    int m_NumberParticlesStripped;
    /**
     * Array of particle IDs of all particles in the decay chain without resonances
     */
    std::vector<int> m_pdgIDStripped;
    /**
     * Vector of indices referring to the particle mother without resonances
     */
    std::vector<int> m_MotherIndexStripped;
    /**
     * Generator label of the decay mode
     */
    int m_MCmode;
    /**
     * True x momentum of the \f$D\f$ meson
     */
    std::vector<double> m_TruePx;
    /**
     * True y momentum of the \f$D\f$ meson
     */
    std::vector<double> m_TruePy;
    /**
     * True z momentum of the \f$D\f$ meson
     */
    std::vector<double> m_TruePz;
    /**
     * True total energy of the \f$D\f$ meson
     */
    std::vector<double> m_TrueEnergy;
};

#endif
