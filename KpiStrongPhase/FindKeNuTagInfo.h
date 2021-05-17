// Martin Duy Tat 25th March 2021, based on code by Yu Zhang
/**
 * FindKeNuTagInfo is class for finding \f$Ke\nu\f$ tag info
 * It reconstructs the electron and kaon on the other side of the signal and looks for missing energy and mass
 */

#ifndef FINDKENUTAGINFO
#define FINDKENUTAGINFO
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

class FindKeNuTagInfo {
  public: 
    /**
     * Default constructor, initializes everything to zero
     */
    FindKeNuTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKeNuTagInfo();
    /**
     * Start looking for \f$Ke\nu\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param DTTool DTagTool object with all the tag information
     */
    StatusCode findKeNuTagInfo(DTagToolIterator DTTool_iter, DTagTool DTTool);
    /**
     * Get the ith component of \f$e^\pm\f$ four-momentum, including FSR
     */
    double GetElectronP(int i) const;
    /**
     * Get the \f$e^\pm\f$ charge
     */
    int GetElectronCharge() const;
    /**
     * Get the ith component of FSR four-momentum
     */
    double GetFSRP(int i) const;
    /**
     * Get the ith component of the \f$K^\pm\f$ four-momentum
     */
    double GetKaonP(int i) const;
    /**
     * Get the \f$K^\pm\f$ charge
     */
    int GetKaonCharge() const;
    /**
     * Get the missing four-momentum
     */
    double GetMissP(int i) const;
    /**
     * Get \f$U_\text{miss} = E_\text{miss} - |p_\text{miss}|\f$
     */
    double GetUMiss() const;
    /**
     * Get the energy of the jth photon shower that is not FSR
     */
    double GetExtraShowerEnergy(int j) const;
    /**
     * Get the number of photons that are not FSR
     */
    int GetNumberGamma() const;
    /**
     * Get the daughter track ID, in the order \f$K\f$ \f$e\f$
     */
    std::vector<int> GetDaughterTrackID() const;
  private:
    /**
     * \f$e^\pm\f$ four-momentum, including FSR
     */
    CLHEP::HepLorentzVector m_ElectronP;
    /**
     * \f$e^\pm\f$ charge
     */
    int m_ElectronCharge;
    /**
     * FSR four-momentum
     */
    CLHEP::HepLorentzVector m_FSRP;
    /**
     * \f$K^\pm\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KaonP;
    /**
     * \f$K^\pm\f$ charge
     */
    int m_KaonCharge;
    /**
     * Missing four-momentum
     */
    CLHEP::HepLorentzVector m_MissP;
    /**
     * \f$U_\text{miss} = E_\text{miss} - |p_\text{miss}|\f$
     */
    double m_UMiss;
    /**
     * Vector of photon shower energies that is not FSR
     */
    std::vector<double> m_ExtraShowerEnergy;
    /**
     * Number of photons that are not FSR
     */
    int m_NumberGamma;
    /**
     * Daughter track ID, in the order \f$K\f$ \f$e\f$
     */
    std::vector<int> m_DaughterTrackID;
};

#endif
