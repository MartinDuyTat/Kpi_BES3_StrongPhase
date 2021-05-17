// Martin Duy Tat 25th March 2021

// KKpipi
#include "KpiStrongPhase/FindKpipipiTagInfo.h"
#include "KpiStrongPhase/FindKS.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
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

FindKpipipiTagInfo::FindKpipipiTagInfo(): m_DaughterTrackID(std::vector<int>(4)), m_KCharge(0), m_Pi1Charge(0), m_Pi2Charge(0), m_Pi3Charge(0), m_12KSFitSuccess(0), m_12DecayLengthVeeVertex(0.0), m_12Chi2VeeVertex(0.0), m_12KSMassVeeVertex(0.0), m_12DecayLengthFit(0.0), m_12DecayLengthErrorFit(0.0), m_12Chi2Fit(0.0), m_13KSFitSuccess(0), m_13DecayLengthVeeVertex(0.0), m_13Chi2VeeVertex(0.0), m_13KSMassVeeVertex(0.0), m_13DecayLengthFit(0.0), m_13DecayLengthErrorFit(0.0), m_13Chi2Fit(0.0) {
}

FindKpipipiTagInfo::~FindKpipipiTagInfo() {
}

StatusCode FindKpipipiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  std::vector<SmartRefVector<EvtRecTrack>::iterator> DaughterTrackIterators(4); // In the order K pi pi pi
  // First find the K daughter
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    if(DTTool.isKaon(*Track_iter)) {
      RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
      DaughterTrackIterators[KAON] = Track_iter;
      m_KP = MDCKalTrack->p4(MASS::K_MASS);
      m_DaughterTrackID[KAON] = (*Track_iter)->trackId();
      m_KCharge = MDCKalTrack->charge();
      break;
    }
  }
  // Go through all tracks and find the pi daughter tracks
  int PionCounter = 0;
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == m_KCharge) {
	DaughterTrackIterators[PION1] = Track_iter;
	m_Pi1P = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[PION1] = (*Track_iter)->trackId();
	m_Pi1Charge = MDCKalTrack->charge();
      } else if(MDCKalTrack->charge() == -1*m_KCharge) {
	if(PionCounter == 0) {
	  DaughterTrackIterators[PION2] = Track_iter;
	  m_Pi2P = MDCKalTrack->p4(MASS::PI_MASS);
	  m_DaughterTrackID[PION2] = (*Track_iter)->trackId();
	  m_Pi2Charge = MDCKalTrack->charge();
	  PionCounter++;
	} else {
	  DaughterTrackIterators[PION3] = Track_iter;
	  m_Pi3P = MDCKalTrack->p4(MASS::PI_MASS);
	  m_DaughterTrackID[PION3] = (*Track_iter)->trackId();
	  m_Pi3Charge = MDCKalTrack->charge();
	}
      }
    }
  }
  // Check if the \f$\pi\pi\f$ pair is a \f$K_S\f$ in disguise
  double Mpipi12 = (m_Pi1P + m_Pi2P).m();
  m_12KSFitSuccess = 0;
  if(Mpipi12 - MASS::KS_MASS < 0.050 && Mpipi12 - MASS::KS_MASS > -0.060) {
    FindKS findKS(false);
    std::vector<int> PionTrackIDs;
    PionTrackIDs.push_back((*DaughterTrackIterators[PION1])->trackId());
    PionTrackIDs.push_back((*DaughterTrackIterators[PION2])->trackId());
    StatusCode statuscode = findKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    m_12DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
    m_12Chi2VeeVertex = findKS.GetChi2VeeVertex();
    m_12KSMassVeeVertex = findKS.GetKSMassVeeVertex();
    if(statuscode == StatusCode::SUCCESS) {
      m_12KSFitSuccess = 1;
      m_12DecayLengthFit = findKS.GetDecayLengthFit();
      m_12DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
      m_12Chi2Fit = findKS.GetChi2Fit();
    }
  }
  // Check if the \f$\pi\pi\f$ pair is a \f$K_S\f$ in disguise
  double Mpipi13 = (m_Pi1P + m_Pi3P).m();
  m_13KSFitSuccess = 0;
  if(Mpipi13 - MASS::KS_MASS < 0.050 && Mpipi13 - MASS::KS_MASS > -0.060) {
    FindKS findKS(false);
    std::vector<int> PionTrackIDs;
    PionTrackIDs.push_back((*DaughterTrackIterators[PION1])->trackId());
    PionTrackIDs.push_back((*DaughterTrackIterators[PION3])->trackId());
    StatusCode statuscode = findKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    m_13DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
    m_13Chi2VeeVertex = findKS.GetChi2VeeVertex();
    m_13KSMassVeeVertex = findKS.GetKSMassVeeVertex();
    if(statuscode == StatusCode::SUCCESS) {
      m_13KSFitSuccess = 1;
      m_13DecayLengthFit = findKS.GetDecayLengthFit();
      m_13DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
      m_13Chi2Fit = findKS.GetChi2Fit();
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<int> FindKpipipiTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}

double FindKpipipiTagInfo::GetKP(int i) const {
  return m_KP[i];
}

double FindKpipipiTagInfo::GetPi1P(int i) const {
  return m_Pi1P[i];
}

double FindKpipipiTagInfo::GetPi2P(int i) const {
  return m_Pi2P[i];
}

double FindKpipipiTagInfo::GetPi3P(int i) const {
  return m_Pi3P[i];
}

int FindKpipipiTagInfo::GetKCharge() const {
  return m_KCharge;
}

int FindKpipipiTagInfo::GetPi1Charge() const {
  return m_Pi1Charge;
}

int FindKpipipiTagInfo::GetPi2Charge() const {
  return m_Pi2Charge;
}

int FindKpipipiTagInfo::GetPi3Charge() const {
  return m_Pi3Charge;
}

int FindKpipipiTagInfo::GetKSFitSuccess12() const {
  return m_12KSFitSuccess;
}

double FindKpipipiTagInfo::GetDecayLengthVeeVertex12() const {
  return m_12DecayLengthVeeVertex;
}

double FindKpipipiTagInfo::GetChi2VeeVertex12() const {
  return m_12Chi2VeeVertex;
}

double FindKpipipiTagInfo::GetKSMassVeeVertex12() const {
  return m_12KSMassVeeVertex;
}

double FindKpipipiTagInfo::GetDecayLengthFit12() const {
  return m_12DecayLengthFit;
}

double FindKpipipiTagInfo::GetDecayLengthErrorFit12() const {
  return m_12DecayLengthErrorFit;
}

double FindKpipipiTagInfo::GetChi2Fit12() const {
  return m_12Chi2Fit;
}

int FindKpipipiTagInfo::GetKSFitSuccess13() const {
  return m_13KSFitSuccess;
}

double FindKpipipiTagInfo::GetDecayLengthVeeVertex13() const {
  return m_13DecayLengthVeeVertex;
}

double FindKpipipiTagInfo::GetChi2VeeVertex13() const {
  return m_13Chi2VeeVertex;
}

double FindKpipipiTagInfo::GetKSMassVeeVertex13() const {
  return m_13KSMassVeeVertex;
}

double FindKpipipiTagInfo::GetDecayLengthFit13() const {
  return m_13DecayLengthFit;
}

double FindKpipipiTagInfo::GetDecayLengthErrorFit13() const {
  return m_13DecayLengthErrorFit;
}

double FindKpipipiTagInfo::GetChi2Fit13() const {
  return m_13Chi2Fit;
}
