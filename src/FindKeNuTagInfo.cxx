// Martin Duy Tat 25th March 2021, based on code by Yu Zhang

// KpiStrongPhase
#include "KpiStrongPhase/FindKeNuTagInfo.h"
#include "KpiStrongPhase/KKpipiUtilities.h"
#include "KpiStrongPhase/ParticleMasses.h"
// Gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Event information
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "ExtEvent/RecExtTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// STL
#include <cmath>
#include <vector>
// ROOT
#include "TMath.h"

FindKeNuTagInfo::FindKeNuTagInfo(): m_ElectronCharge(0), m_KaonCharge(0), m_UMiss(0.0), m_ExtraShowerEnergy(0.0), m_NumberGamma(0), m_DaughterTrackID(std::vector<int>(2)) {
}

FindKeNuTagInfo::~FindKeNuTagInfo() {
}

StatusCode FindKeNuTagInfo::findKeNuTagInfo(DTagToolIterator DTTool_iter, DTagTool DTTool) {
  // Prepare message service
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindKeNuTagInfo");
  // Get tracks on the other side of the reconstructed D meson
  SmartRefVector<EvtRecTrack> OtherTracks = (*DTTool_iter)->otherTracks();
  // Loop over all tracks on the other side to find pi+ pi-
  int NumberElectronTracks = 0, NumberKaonTracks = 0;
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = OtherTracks.begin(); Track_iter != OtherTracks.end(); Track_iter++) {
    // First check if track is valid
    if(!(*Track_iter)->isMdcTrackValid() || !(*Track_iter)->isMdcKalTrackValid()) {
      continue;
    }
    if(DTTool.isGoodTrack(*Track_iter)) {
      RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
      if(DTTool.isElectron(*Track_iter)) {
	MDCKalTrack->setPidType(RecMdcKalTrack::electron);
	m_ElectronP = MDCKalTrack->p4(MASS::ELECTRON_MASS);
	m_ElectronCharge = MDCKalTrack->charge();
	m_DaughterTrackID[1] = (*Track_iter)->trackId();
	NumberElectronTracks++;
      } else if(DTTool.isKaon(*Track_iter)) {
	MDCKalTrack->setPidType(RecMdcKalTrack::kaon);
	m_KaonP = MDCKalTrack->p4(MASS::K_MASS);
	m_KaonCharge = MDCKalTrack->charge();
	m_DaughterTrackID[0] = (*Track_iter)->trackId();
	NumberKaonTracks++;
      } else {
	return StatusCode::FAILURE;
      }
    }
  }
  // If more than one electron or kaon is found, or if their charges are not opposite, skip event
  if(NumberElectronTracks != 1 || NumberKaonTracks != 1 || m_ElectronCharge*m_KaonCharge != -1) {
    return StatusCode::FAILURE;
  }
  // Get showers on the other side of the reconstructed D meson
  SmartRefVector<EvtRecTrack> OtherShowers = (*DTTool_iter)->otherShowers();
  // Loop over all showers to find FSR photons
  for(SmartRefVector<EvtRecTrack>::iterator Shower_iter = OtherShowers.begin(); Shower_iter != OtherShowers.end(); Shower_iter++) {
    // Check if shower is valid
    if(!(*Shower_iter)->isEmcShowerValid()) {
      continue;
    }
    // Get reconstructed EMC shower
    RecEmcShower *EMCShower = (*Shower_iter)->emcShower();
    if(EMCShower->module() == 1 && EMCShower->energy() < 0.025) {
      // Shower in the barrel must have energy larger than 25 MeV
      continue;
    } else if(EMCShower->module() != 1 && EMCShower->energy() < 0.050) {
      // Shower in endcap must have energy larger than 50 MeV
      continue;
    }
    if(EMCShower->time() < 0 || EMCShower->time() > 14) {
      // EMC shower time requirement 0 <= T <= 14 (in units of 50 ns)
      continue;
    }
    // Get EMC position of shower
    CLHEP::Hep3Vector EMCPosition(EMCShower->x(), EMCShower->y(), EMCShower->z());
    // Initialize angles to their maximum
    double Theta;
    double Phi;
    double Angle;
    // If the nearest charged track is separated by less than 20 degrees from the shower, skip event
    if(!KKpipiUtilities::GetPhotonAngularSeparation(EMCPosition, Angle, Theta, Phi) || Angle < 20.0*TMath::Pi()/180.0) {
      continue;
    }
    // Get four-momentum of shower
    CLHEP::HepLorentzVector ShowerP = KKpipiUtilities::GetPhoton4Vector(EMCShower->energy(), EMCShower->theta(), EMCShower->phi());
    // Calculate angle between lepton and shower
    double PhotonLeptonAngle = m_ElectronP.vect().theta(ShowerP.vect());
    if(PhotonLeptonAngle < 5.0*TMath::Pi()/180.0) {
      // If photon and lepton are less than 5 degrees apart, this must be FSR
      m_FSRP += ShowerP;
    } else {
      // If not, save shower four-momentum for later background study
      m_ExtraShowerEnergy.push_back(ShowerP.e());
      m_NumberGamma++;
    }
  }
  m_ElectronP += m_FSRP;
  m_MissP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), m_ElectronP + m_KaonP, (*DTTool_iter)->beamE());
  m_UMiss = m_MissP.e() - m_MissP.vect().mag();
  return StatusCode::SUCCESS;
}

double FindKeNuTagInfo::GetElectronP(int i) const {
  return m_ElectronP[i];
}

int FindKeNuTagInfo::GetElectronCharge() const {
  return m_ElectronCharge;
}

double FindKeNuTagInfo::GetFSRP(int i) const {
  return m_FSRP[i];
}

double FindKeNuTagInfo::GetKaonP(int i) const {
  return m_KaonP[i];
}

int FindKeNuTagInfo::GetKaonCharge() const {
  return m_KaonCharge;
}

double FindKeNuTagInfo::GetMissP(int i) const {
  return m_MissP[i];
}

double FindKeNuTagInfo::GetUMiss() const {
  return m_UMiss;
}

double FindKeNuTagInfo::GetExtraShowerEnergy(int j) const {
  return m_ExtraShowerEnergy[j];
}

int FindKeNuTagInfo::GetNumberGamma() const {
  return m_NumberGamma;
}

std::vector<int> FindKeNuTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}
