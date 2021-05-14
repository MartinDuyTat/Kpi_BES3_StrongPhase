// Martin Duy Tat 23rd March 2021
/**
 * PIDTruth is a class used to check the true PID values of reconstructed tracks
 */

// Gaudi
#include "GaudiKernel/Algorithm.h"
// STL
#include <vector>

#ifndef PIDTRUTH
#define PIDTRUTH

class PIDTruth {
  public:
    /**
     * Constructor that takes in a vector of track IDs of reconstructed tracks and an Algorithm object and saves these
     * The charged track IDs must come first in the input vector!
     * @param TrackID Vector of track IDs
     * @param NumberCharged Number of charged tracks
     * @param algorithm Pointer to the Algorithm object that created this class, need this to access the MC info!
     */
    PIDTruth(const std::vector<int> &TrackID, int NumberCharged, const Algorithm *algorithm);
    /**
     * Function from Alex Gilman that matches charged reconstructed particles with generator particles by comparing the hits of the reconstructed object and the truth level trajectory
     * @param tkID Track ID of particle of interest
     * @param mcPDG Set to -1 to find PID of the particle of interest, set to the reconstructed PID to check if the PID has been correctly assigned, set to 0 if we're not interested in this
     * @param mcParPDG Set to -1 to find PID of the parent, set to the reconstructed PID to check if the PID has been correctly assigned to the parent, set to 0 if we're not interested in this
     * @param GParPDG Set to the PID of any parent or grand parent if we want to check if the particle of interest originated from this particle through a (potentially long) decay chain, set to 0 if we're not interested in this
     * @return Returns 0 for false, 1 for true or the PDG PID number where appropriate
     */
    int MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const;
    /**
     * Function from Alex Gilman that matches neutral reconstructed particles with generator particles by comparing the hits of the reconstructed object and the truth level trajectory
     * @param tkID Track ID of particle of interest
     * @param mcPDG Set to -1 to find PID of the particle of interest, set to the reconstructed PID to check if the PID has been correctly assigned, set to 0 if we're not interested in this
     * @param mcParPDG Set to -1 to find PID of the parent, set to the reconstructed PID to check if the PID has been correctly assigned to the parent, set to 0 if we're not interested in this
     * @param GParPDG Set to the PID of any parent or grand parent if we want to check if the particle of interest originated from this particle through a (potentially long) decay chain, set to 0 if we're not interested in this
     * @return Returns 0 for false, 1 for true or the PDG PID number where appropriate
     */
    int MCSHPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const;
    /**
     * Function that determines which $D$ meson that a particle originated from
     * @param TrackID Track ID of the particle of interest
     * @param Charged True if particle is charged, false if particle is neutral
     * @return 421 for \f$D^0\f$, -421 for \f$\bar{D^0}\f$ and 0 for unknown
     */
    int FindDOrigin(int TrackID, bool Charged) const;
    /**
     * Function that checks if all the particles originate from the same \f$D\f$ meson
     * @param Vector of track IDs that should not be included in this check
     */
    bool SameDMother(const std::vector<int> &IgnoreTrackID = std::vector<int>()) const;
    /**
     * Function that determines the true PID of the particles of interest
     * If some particle is to not be included in the truth matching, put 0 as the input
     * @param ParticleID Input vector of reconstructed PIDs, the function will replace these with the MC truth
     * @return Returns true if the reconstructed PIDs match perfectly with the MC truth
     */
    bool FindTrueID(std::vector<int> &ParticleID) const;
    /**
     * Function for obtaining the true mother PID
     * @param TrackID Track ID of particle of interest
     * @param Charge True if particle is charged, false if it's neutral
     * @return Returns the ID of the mother particle
     */
    int GetTrueMotherID(int TrackID, bool Charged) const;
  private:
    /**
     * Vector with track IDs of reconstructed daughters
     */
    std::vector<int> m_TrackID;
    /**
     * Number of charged particles
     */
    int m_NumberCharged;
    /**
     * Pointer to the Algorithm object that created this class, need this to access the MC info!
     */
    const Algorithm *m_algorithm;
};

#endif
