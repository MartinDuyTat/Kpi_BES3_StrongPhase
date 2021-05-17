// Martin Duy Tat 12th February 2021, based on code by Yu Zhang

// KKpipi file
#include "KpiStrongPhase/FindPi0Eta.h"
#include "KpiStrongPhase/KKpipiUtilities.h"
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
#include "EmcRecEventModel/RecEmcShower.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "EvtRecEvent/EvtRecEtaToGG.h"
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/VertexParameter.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/WTrackParameter.h"
// STL
#include<vector>
#include<string>
// ROOT
#include "TMath.h"
// Particle masses
#include "KpiStrongPhase/ParticleMasses.h"

FindPi0Eta::FindPi0Eta(int npi0eta, std::string Particle): m_Chi2Fit(0.0), m_npi0eta(npi0eta), m_Particle(Particle) {
  m_HighEPhotonP.resize(m_npi0eta);
  m_LowEPhotonP.resize(m_npi0eta);
  m_HighEPhotonPConstrained.resize(m_npi0eta);
  m_LowEPhotonPConstrained.resize(m_npi0eta);
  m_Chi2Fit.resize(m_npi0eta);
}

FindPi0Eta::~FindPi0Eta() {
}

StatusCode FindPi0Eta::findPi0Eta(DTagToolIterator &DTTool_iter, DTagTool &DTTool) {
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindPi0Eta");
  IDataProviderSvc *EventService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventService);
  SmartDataPtr<EvtRecPi0Col> RecPi0Col(EventService, "/Event/EvtRec/EvtRecPi0Col");
  SmartDataPtr<EvtRecEtaToGGCol> RecEtaCol(EventService, "/Event/EvtRec/EvtRecEtaToGGCol");
  if(m_Particle != "pi0" && m_Particle != "eta") {
    log << "Particle must be pi0 or eta" << endreq;
  } else if(!RecPi0Col && m_Particle == "pi0") {
    log << "Could not find EvtRecPi0Col" << endreq;
    return StatusCode::FAILURE;
  } else if(!RecEtaCol && m_Particle == "eta") {
    log << "Could not find EvtRecEtaCol" << endreq;
    return StatusCode::FAILURE;
  }
  if(m_Particle == "pi0") {
    int i = 0;
    // Get the track ID of the tagged pi0 candidates
    std::vector<int> Pi0TrackIDs = DTTool.pi0Id(DTTool_iter, m_npi0eta);
    // Loop over all pi0 candidates found by EvtRecPi0
    for(EvtRecPi0Col::iterator Pi0_iter = RecPi0Col->begin(); Pi0_iter != RecPi0Col->end(); Pi0_iter++) {
      // Check if pi0 candidate from the DTagTool is the same as the candidate found by EvtRecPi0
      if(Pi0_iter - RecPi0Col->begin() != Pi0TrackIDs[i]) {
	continue;
      }
      // Get EM shower four-momenta of photons
      RecEmcShower *HighEPhotonShower = const_cast<EvtRecTrack*>((*Pi0_iter)->hiEnGamma())->emcShower();
      RecEmcShower *LowEPhotonShower = const_cast<EvtRecTrack*>((*Pi0_iter)->loEnGamma())->emcShower();
      m_HighEPhotonP[i] = KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi());
      m_LowEPhotonP[i] = KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi());
      // Get kinematically constrained four-momenta of photons
      m_HighEPhotonPConstrained[i] = (*Pi0_iter)->hiPfit();
      m_LowEPhotonPConstrained[i] = (*Pi0_iter)->loPfit();
      m_Chi2Fit[i] = (*Pi0_iter)->chisq();
      // Get photon track IDs
      m_HighEPhotonTrackID.push_back((*Pi0_iter)->hiEnGamma()->trackId());
      m_LowEPhotonTrackID.push_back((*Pi0_iter)->loEnGamma()->trackId());
      i++;
      if(i == m_npi0eta) {
	return StatusCode::SUCCESS;
      }
    }
  } else if(m_Particle == "eta") {
    int i = 0;
    // Get the track ID of the tagged eta candidates
    std::vector<int> EtaTrackIDs = DTTool.etaId(DTTool_iter, m_npi0eta);
    // Loop over all eta candidates found by EvtRecEtaToGG
    for(EvtRecEtaToGGCol::iterator Eta_iter = RecEtaCol->begin(); Eta_iter != RecEtaCol->end(); Eta_iter++) {
      // Check if eta candidate from the DTagTool is the same as the candidate found by EvtRecEtaToGG
      if(Eta_iter - RecEtaCol->begin() != EtaTrackIDs[i]) {
	continue;
      }
      // Get EM shower four-momenta of photons
      RecEmcShower *HighEPhotonShower = const_cast<EvtRecTrack*>((*Eta_iter)->hiEnGamma())->emcShower();
      RecEmcShower *LowEPhotonShower = const_cast<EvtRecTrack*>((*Eta_iter)->loEnGamma())->emcShower();
      m_HighEPhotonP[i] = KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi());
      m_LowEPhotonP[i] = KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi());
      // Get kinematically constrained four-momenta of photons
      m_HighEPhotonPConstrained[i] = (*Eta_iter)->hiPfit();
      m_LowEPhotonPConstrained[i] = (*Eta_iter)->loPfit();
      m_Chi2Fit[i] = (*Eta_iter)->chisq();
      // Get photon track IDs
      m_HighEPhotonTrackID.push_back((*Eta_iter)->hiEnGamma()->trackId());
      m_LowEPhotonTrackID.push_back((*Eta_iter)->loEnGamma()->trackId());
      i++;
      if(i == m_npi0eta) {
	return StatusCode::SUCCESS;
      }
    }
  }
  log << "Could not find the requested number of " << m_Particle << " candidates" << endreq;
  return StatusCode::FAILURE;
}

double FindPi0Eta::GetHighEPhotonP(int i, int pi0eta_index) const {
  return m_HighEPhotonP[pi0eta_index][i];
}

double FindPi0Eta::GetLowEPhotonP(int i, int pi0eta_index) const {
  return m_LowEPhotonP[pi0eta_index][i];
}

double FindPi0Eta::GetMgammagamma(int pi0eta_index) const {
  return (m_HighEPhotonP[pi0eta_index] + m_LowEPhotonP[pi0eta_index]).m();
}

double FindPi0Eta::GetHighEPhotonPConstrained(int i, int pi0eta_index) const {
  return m_HighEPhotonPConstrained[pi0eta_index][i];
}

double FindPi0Eta::GetLowEPhotonPConstrained(int i, int pi0eta_index) const {
  return m_LowEPhotonPConstrained[pi0eta_index][i];
}

double FindPi0Eta::GetChi2Fit(int pi0eta_index) const {
  return m_Chi2Fit[pi0eta_index];
}

int FindPi0Eta::GetHighEPhotonTrackID(int pi0eta_index) const {
  return m_HighEPhotonTrackID[pi0eta_index];
}

int FindPi0Eta::GetLowEPhotonTrackID(int pi0eta_index) const {
  return m_LowEPhotonTrackID[pi0eta_index];
}
