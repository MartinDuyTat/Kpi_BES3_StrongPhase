// Martin Duy Tat 23rd March 2021
/**
 * StripResonances is a class used to strip away all intermediate resonances from an array of MC truth particle IDs and their mother indices
 * The only particles left in the array are \f$\psi(3770)\f$, \f$D^0\f$, \f$\bar{D^0}\f$, \f$K^\pm\f$, \f$K_S\f$, \f$K_L\f$, \f$\pi^\pm$, \f$\pi^0\f$, \f$e^\pm\f$, \f$\mu^\pm\f$, \f$\nu\f$, \f$p\f$, \f$\bar{p}\f$ and \f$\gamma\f$
 * One can also add an optional list of resonances to not remove in the constructor
 */

#include <vector>

#ifndef STRIPRESONANCES
#define STRIPRESONANCES

class StripResonances {
  public:
    /**
     * Constructor that sets up the list of particles that are not removed
     * @param KeepList List of particle IDs that should be kept, useful if one wants to study \f$\phi\f$ resonances for example, which ordinarily would be removed
     */
    StripResonances(const std::vector<int> &KeepList = std::vector<int>());
    /**
     * Function that return the index of the first resonance it can find
     * @param pidID Vector of particle IDs
     * @return Returns the index of the resonance, return -1 if none is found
     */
    int FindResonanceIndex(const std::vector<int> &pidID) const;
    /**
     * Function that removes the resonance at a specific index
     * The input array is passed by reference because the resonance is removed directly
     * @param Index of resonance to be removed
     * @param pidID Vector of particle IDs
     * @param MotherIndex Vector of the mother indices
     */
    void RemoveResonance(int ResonanceIndex, std::vector<int> &pidID, std::vector<int> &MotherIndex) const;
    /**
     * Function that removes all intermediate resonances
     * @param pidID Vector of particle IDs
     * @param MotherIndex Vector of the mother indices
     */
    void RemoveIntermediateResonances(std::vector<int> &pidID, std::vector<int> &MotherIndex) const;
  private:
    /**
     * List of particles that should not be removed
     */
    std::vector<int> m_KeepList;
};

#endif



