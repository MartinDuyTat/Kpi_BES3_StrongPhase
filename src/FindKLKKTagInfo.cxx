// Martin Duy Tat 24th May 2021

// KpiStrongPhase
#include "KpiStrongPhase/FindKLKKTagInfo.h"
#include "KpiStrongPhase/ParticleMasses.h"
#include "KpiStrongPhase/KpiStrongPhaseUtilities.h"
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
#include "VertexFit/KalmanKinematicFit.h"
// STL
#include <cmath>
#include <vector>
// ROOT
#include "TMath.h"

FindKLKKTagInfo::FindKLKKTagInfo(): m_KalmanFitSuccess(0), m_KalmanFitChi2(0), m_NearestShowerEnergy(0.0), m_NearestShowerCosAngle(0.0) {
}

FindKLKKTagInfo::~FindKLKKTagInfo() {
}

StatusCode FindKLKKTagInfo::findKL(DTagToolIterator DTTool_iter, DTagTool DTTool) {
  // Prepare message service
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindKLKKTagInfo");
  // Prepare event data service
  IDataProviderSvc *EventDataService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventDataService);
  // Prepare pi0 service
  SmartDataPtr<EvtRecPi0Col> evtRecPi0Col(EventDataService, "/Event/EvtRec/EvtRecPi0Col");
  if(!evtRecPi0Col) {
    log << MSG::ERROR << "EvtRecPi0Col not found" << endreq;
  }
  // Prepare eta service
  SmartDataPtr<EvtRecEtaToGGCol> evtRecEtaToGGCol(EventDataService, "/Event/EvtRec/EvtRecEtaToGGCol");
  if(!evtRecEtaToGGCol) {
    log << MSG::ERROR << "EvtRecEtaToGGCol not found" << endreq;
  }
  // Get tracks on the other side of the reconstructed D meson
  SmartRefVector<EvtRecTrack> OtherTracks = (*DTTool_iter)->otherTracks();
  // Loop over all tracks on the other side to find pi+ pi-
  std::vector<RecMdcKalTrack*> KalmanTracks(2); //In the order K+ K-
  int NumberKPlusTracks = 0, NumberKMinusTracks = 0;
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = OtherTracks.begin(); Track_iter != OtherTracks.end(); Track_iter++) {
    // First check if track is valid
    if(!(*Track_iter)->isMdcTrackValid() || !(*Track_iter)->isMdcKalTrackValid()) {
      continue;
    }
    if(DTTool.isGoodTrack(*Track_iter)) {
      if(DTTool.isKaon(*Track_iter)) {
	RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
	MDCKalTrack->setPidType(RecMdcKalTrack::kaon);
	if(MDCKalTrack->charge() == +1) {
	  NumberKPlusTracks++;
	  m_KPlusP = MDCKalTrack->p4(MASS::K_MASS);
	  KalmanTracks[0] = MDCKalTrack;
	  m_DaughterTrackID[0] = (*Track_iter)->trackId();
	} else if(MDCKalTrack->charge() == -1) {
	  NumberKMinusTracks++;
	  m_KMinusP = MDCKalTrack->p4(MASS::K_MASS);
	  KalmanTracks[1] = MDCKalTrack;
	  m_DaughterTrackID[1] = (*Track_iter)->trackId();
	} else {
	  return StatusCode::FAILURE;
	}
      } else {
	return StatusCode::FAILURE;
      }
    }
  }
  // If K+ K- pair is not found, reject event
  if(NumberKPlusTracks != 1 || NumberKMinusTracks != 1) {
    return StatusCode::FAILURE;
  }
  // Get missing momentum, which is the KL momentum
  m_KLongP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), m_KPlusP + m_KMinusP, (*DTTool_iter)->beamE());
  // Check if there are any showers already present on the reconstructed side
  bool ShowersUsed;
  SmartRefVector<EvtRecTrack> Showers = (*DTTool_iter)->showers();
  if(Showers.size() == 0) {
    // If no showers, it's probably a Kpi or Kpipipi tag on the reconstructed side
    ShowersUsed = false;
  } else if(Showers.size() > 2) {
    // There shouldn't be more than 2 showers or just 1 shower
    return StatusCode::FAILURE;
  } else if(Showers.size() == 1) {
    log << MSG::ERROR << "Found a single shower when looking for KL" << endreq;
    return StatusCode::FAILURE;
  } else {
    // If 2 showers, it's probably a Kpipi0 tag on the reconstructed side
    ShowersUsed = true;
  }
  // Look for pi0
  int NumberPi0 = 0;
  for(EvtRecPi0Col::iterator Pi0_iter = evtRecPi0Col->begin(); Pi0_iter != evtRecPi0Col->end(); Pi0_iter++) {
    // Get photon tracks...?
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->hiEnGamma());
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->loEnGamma());
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = HighEnergyPhotonTrack->emcShower();
    RecEmcShower *LowEPhotonShower = LowEnergyPhotonTrack->emcShower();
    // Get photon track ID
    int HighEnergyPhotonTrackID = HighEnergyPhotonTrack->trackId();
    int LowEnergyPhotonTrackID = LowEnergyPhotonTrack->trackId();
    // If reconstructed side constains a shower, make sure this is not the same shower
    if(ShowersUsed) {
      if(HighEnergyPhotonTrackID == Showers[0]->trackId() || HighEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
      if(LowEnergyPhotonTrackID == Showers[0]->trackId() || LowEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
    }
    NumberPi0++;
  }
  // If there are no pi0, reject event
  if(NumberPi0 != 0) {
    return StatusCode::FAILURE;
  }
  // Look for eta
  int NumberEta = 0;
  for(EvtRecEtaToGGCol::iterator Eta_iter = evtRecEtaToGGCol->begin(); Eta_iter != evtRecEtaToGGCol->end(); Eta_iter++) {
    // Get photon tracks...?
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Eta_iter)->hiEnGamma());
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Eta_iter)->loEnGamma());
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = HighEnergyPhotonTrack->emcShower();
    RecEmcShower *LowEPhotonShower = LowEnergyPhotonTrack->emcShower();
    // Get photon track ID
    int HighEnergyPhotonTrackID = HighEnergyPhotonTrack->trackId();
    int LowEnergyPhotonTrackID = LowEnergyPhotonTrack->trackId();
    // If reconstructed side constains a shower, make sure this is not the same shower
    if(ShowersUsed) {
      if(HighEnergyPhotonTrackID == Showers[0]->trackId() || HighEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
      if(LowEnergyPhotonTrackID == Showers[0]->trackId() || LowEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
    }
    NumberEta++;
  }
  // Get showers on the other side of the reconstructed D meson
  SmartRefVector<EvtRecTrack> OtherShowers = (*DTTool_iter)->otherShowers();
  // Loop over all showers to find nearest photon
  double LargestCosAngle = -1.0;
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
    double PhotonP = KKpipiUtilities::GetPhoton4Vector(EMCShower->energy(), EMCShower->theta(), EMCShower->phi());
    double CosAngle = PhotonP.vect().unit().dot(m_KLong.vect().unit());
    if(CosAngle > LargestCosAngle) {
      LargestCosAngle = CosAngle;
      m_NearestShowerEnergy = PhotonP[3];
      m_NearestShowerCosAngle = CosAngle;
    }
  }
  // Do a Kalman fit of the KL K+ K- daughters
  WTrackParameter WTrackKplus(MASS::K_MASS, KalmanTracks[0]->getZHelix(), KalmanTracks[0]->getZError());
  WTrackParameter WTrackKminus(MASS::K_MASS, KalmanTracks[1]->getZHelix(), KalmanTracks[1]->getZError());
  KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
  KalmanFit->init();
  KalmanFit->AddTrack(0, WTrackKplus);
  KalmanFit->AddTrack(1, WTrackKminus);
  KalmanFit->AddMissTrack(2, MASS::KS_MASS, m_KLongP);
  KalmanFit->AddResonance(0, MASS::D_MASS, 0, 1, 2);
  m_KalmanFitSuccess = KalmanFit->Fit();
  if(m_KalmanFitSuccess) {
    m_KalmanFitChi2 = KalmanFit->chisq();
    m_KPlusPKalmanFit = KalmanFit->pfit(0);
    m_KMinusPKalmanFit = KalmanFit->pfit(1);
    m_KLongPKalmanFit = KalmanFit->pfit(2);
  }
  return StatusCode::SUCCESS;
}

double FindKLKKTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKLKKTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}

double FindKLKKTagInfo::GetKLongP(int i) const {
  return m_KLongP[i];
}

double FindKLKKTagInfo::GetMMiss2() const {
  return m_KLongP.M2();
}

double FindKLKKTagInfo::GetKPlusPKalmanFit(int i) const {
  return m_KPlusPKalmanFit[i];
}

double FindKLKKTagInfo::GetKMinusPKalmanFit(int i) const {
  return m_KMinusPKalmanFit[i];
}

double FindKLKKTagInfo::GetKLongPKalmanFit(int i) const {
  return m_KLongPKalmanFit[i];
}

int FindKLKKTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKLKKTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKLKKTagInfo::GetNearestShowerEnergy() const {
  return m_NearestShowerEnergy;
}

double FindKLKKTagInfo::GetNearestShowerCosAngle() const {
  return m_NearestShowerCosAngle;
}
