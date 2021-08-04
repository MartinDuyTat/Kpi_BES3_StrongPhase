// Martin Duy Tat 4th August 2021

// KpiStrongPhase
#include "KpiStrongPhase/FindKSpipiTagInfo.h"
#include "KpiStrongPhase/FindKS.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// ROOT
#include "TMath.h"
// STL
#include<vector>
#include<string>
// Particle masses
#include "KpiStrongPhase/ParticleMasses.h"

FindKSpipiTagInfo::FindKSpipiTagInfo(): m_DaughterTrackID(std::vector<int>(4)), m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0), m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2Fit(0.0), m_KalmanFitSuccess(0), m_KalmanFitChi2(0.0) {
}

FindKSpipiTagInfo::~FindKSpipiTagInfo() {
}

StatusCode FindKSpipiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // Find the actual KS candidate
  FindKS findKS(true);
  StatusCode status = findKS.findKS(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  // Fill out the information about the KS candidate
  m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_DecayLengthFit = findKS.GetDecayLengthFit();
  m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
  m_Chi2Fit = findKS.GetChi2Fit();
  m_KSPiPlusP = CLHEP::HepLorentzVector(findKS.GetKSPiPlusP(0), findKS.GetKSPiPlusP(1), findKS.GetKSPiPlusP(2), findKS.GetKSPiPlusP(3));
  m_KSPiMinusP = CLHEP::HepLorentzVector(findKS.GetKSPiMinusP(0), findKS.GetKSPiMinusP(1), findKS.GetKSPiMinusP(2), findKS.GetKSPiMinusP(3));
  m_KShortP = findKS.GetKShortPFit();
  // Get the track ID of the KS daughters
  std::vector<int> KSDaughterTrackIDs = findKS.GetDaughterTrackIDs();
  m_DaughterTrackID[0] = KSDaughterTrackIDs[0];
  m_DaughterTrackID[1] = KSDaughterTrackIDs[1];
  // Get all tracks
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  std::vector<SmartRefVector<EvtRecTrack>::iterator> DaughterTrackIterators(2); // In the order pi+ pi-
  std::vector<RecMdcKalTrack*> KalmanTracks(2); //In the order pi+ pi-
  std::vector<RecMdcKalTrack*> KSDaughterKalmanTracks;
  // Loop over all tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    // If track is from KS daughters, skip
    if(std::find(KSDaughterTrackIDs.begin(), KSDaughterTrackIDs.end(), (*Track_iter)->trackId()) != KSDaughterTrackIDs.end()) {
      KSDaughterKalmanTracks.push_back(MDCKalTrack);
      continue;
    }
    // Fill out track information
    if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[PIPLUS] = Track_iter;
	KalmanTracks[PIPLUS] = MDCKalTrack;
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[2] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[PIMINUS] = Track_iter;
	KalmanTracks[PIMINUS] = MDCKalTrack;
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[3] = (*Track_iter)->trackId();
      }
    }
  }
  // Do a Kalman kinematic fit of the KS daughter tracks, and constrain the KS mass to its PDG value
  WTrackParameter WTrackKSPIplus(MASS::PI_MASS, KSDaughterKalmanTracks[0]->getZHelix(), KSDaughterKalmanTracks[0]->getZError());
  WTrackParameter WTrackKSPIminus(MASS::PI_MASS, KSDaughterKalmanTracks[1]->getZHelix(), KSDaughterKalmanTracks[1]->getZError());
  KalmanKinematicFit *KSKalmanFit = KalmanKinematicFit::instance();
  KSKalmanFit->init();
  KSKalmanFit->AddTrack(0, WTrackKSPIplus);
  KSKalmanFit->AddTrack(1, WTrackKSPIminus);
  KSKalmanFit->AddResonance(0, MASS::KS_MASS, 0, 1);
  bool KSKalmanFitSuccess = KSKalmanFit->Fit();
  // Do a Kalman kinematic fit of the D daughter tracks, and constrain the KS and D masses to their PDG values
  if(KSKalmanFitSuccess) {
    KSKalmanFit->BuildVirtualParticle(0);
    WTrackParameter WTrackKS = KSKalmanFit->wVirtualTrack(0);
    WTrackParameter WTrackPiplus(MASS::PI_MASS, KalmanTracks[PIPLUS]->getZHelix(), KalmanTracks[PIPLUS]->getZError());
    WTrackParameter WTrackPiminus(MASS::PI_MASS, KalmanTracks[PIMINUS]->getZHelix(), KalmanTracks[PIMINUS]->getZError());
    KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
    KalmanFit->init();
    KalmanFit->AddTrack(PIPLUS, WTrackPiplus);
    KalmanFit->AddTrack(PIMINUS, WTrackPiminus);
    KalmanFit->AddTrack(KSHORT, WTrackKS);
    KalmanFit->AddResonance(0, MASS::D_MASS, PIPLUS, PIMINUS, KSHORT);
    m_KalmanFitSuccess = KalmanFit->Fit();
    if(m_KalmanFitSuccess) {
      m_KalmanFitChi2 = KalmanFit->chisq();
      m_KShortPKalmanFit = KalmanFit->pfit(KSHORT);
      m_PiPlusPKalmanFit = KalmanFit->pfit(PIPLUS);
      m_PiMinusPKalmanFit = KalmanFit->pfit(PIMINUS);
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<int> FindKSpipiTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}

double FindKSpipiTagInfo::GetDecayLengthVeeVertex() const {
  return m_DecayLengthVeeVertex;
}

double FindKSpipiTagInfo::GetChi2VeeVertex() const {
  return m_Chi2VeeVertex;
}

double FindKSpipiTagInfo::GetKSMassVeeVertex() const {
  return m_KSMassVeeVertex;
}

double FindKSpipiTagInfo::GetDecayLengthFit() const {
  return m_DecayLengthFit;
}

double FindKSpipiTagInfo::GetDecayLengthErrorFit() const {
  return m_DecayLengthErrorFit;
}

double FindKSpipiTagInfo::GetChi2Fit() const {
  return m_Chi2Fit;
}

double FindKSpipiTagInfo::GetKSPiPlusP(int i) const {
  return m_KSPiPlusP[i];
}

double FindKSpipiTagInfo::GetKSPiMinusP(int i) const {
  return m_KSPiMinusP[i];
}

double FindKSpipiTagInfo::GetKShortP(int i) const {
  return m_KShortP[i];
}

double FindKSpipiTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKSpipiTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

int FindKSpipiTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKSpipiTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKSpipiTagInfo::GetKShortPKalmanFit(int i) const {
  return m_KShortPKalmanFit[i];
}

double FindKSpipiTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKSpipiTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}
