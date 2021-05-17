// Martin Duy Tat 17th May 2021

// KpiStrongPhase
#include "KpiStrongPhase/KSKKVersusKpipipiDoubleTag.h"
#include "KpiStrongPhase/FindKSKKTagInfo.h"
#include "KpiStrongPhase/FindKpipipiTagInfo.h"
#include "KpiStrongPhase/FindMCInfo.h"
#include "KpiStrongPhase/PIDTruth.h"
// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
// Event information
#include "EventModel/Event.h"
#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>
#include<string>

KSKKVersusKpipipiDoubleTag::KSKKVersusKpipipiDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSKKVersusKpipipiDoubleTag::~KSKKVersusKpipipiDoubleTag() {
}

StatusCode KSKKVersusKpipipiDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSKK vs Kpipipi Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KPISTRONGPHASE/KSKKVersusKpipipiDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KPISTRONGPHASE/KSKKVersusKpipipiDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KSKK vs D->Kpipipi events");
    if(m_tuple) {
      status = m_tuple->addItem("Run", m_RunNumber);
      status = m_tuple->addItem("Event", m_EventNumber);
      status = m_tuple->addItem("NumberOfParticles", m_NumberParticles, 0, 100);
      status = m_tuple->addIndexedItem("ParticleIDs", m_NumberParticles, m_pdgID);
      status = m_tuple->addIndexedItem("MotherIndex", m_NumberParticles, m_MotherIndex);
      status = m_tuple->addItem("NumberOfParticlesStripped", m_NumberParticlesStripped, 0, 100);
      status = m_tuple->addIndexedItem("ParticleIDsStripped", m_NumberParticlesStripped, m_pdgIDStripped);
      status = m_tuple->addIndexedItem("MotherIndexStripped", m_NumberParticlesStripped, m_MotherIndexStripped);
      status = m_tuple->addItem("MCmode", m_MCmode);
      status = m_tuple->addIndexedItem("True_Px", m_NumberParticles, m_TruePx);
      status = m_tuple->addIndexedItem("True_Py", m_NumberParticles, m_TruePy);
      status = m_tuple->addIndexedItem("True_Pz", m_NumberParticles, m_TruePz);
      status = m_tuple->addIndexedItem("True_Energy", m_NumberParticles, m_TrueEnergy);
      status = m_tuple->addItem("SignalDMass", m_SignalDMass);
      status = m_tuple->addItem("SignalMBC", m_SignalMBC);
      status = m_tuple->addItem("SignalDeltaE", m_SignalDeltaE);
      status = m_tuple->addItem("SignalBeamE", m_SignalBeamE);
      status = m_tuple->addItem("TagDMass", m_TagDMass);
      status = m_tuple->addItem("TagMBC", m_TagMBC);
      status = m_tuple->addItem("TagDeltaE", m_TagDeltaE);
      status = m_tuple->addItem("TagBeamE", m_TagBeamE);
      status = m_tuple->addItem("SignalDpx", m_SignalDpx);
      status = m_tuple->addItem("SignalDpy", m_SignalDpy);
      status = m_tuple->addItem("SignalDpz", m_SignalDpz);
      status = m_tuple->addItem("SignalDenergy", m_SignalDenergy);
      status = m_tuple->addItem("TagDpx", m_TagDpx);
      status = m_tuple->addItem("TagDpy", m_TagDpy);
      status = m_tuple->addItem("TagDpz", m_TagDpz);
      status = m_tuple->addItem("TagDenergy", m_TagDenergy);
      status = m_tuple->addItem("SignalKSDecayLengthVeeVertex", m_SignalDecayLengthVeeVertex);
      status = m_tuple->addItem("SignalKSChi2VeeVertex", m_SignalChi2VeeVertex);
      status = m_tuple->addItem("SignalKSMassVeeVertex", m_SignalKSMassVeeVertex);
      status = m_tuple->addItem("SignalKSDecayLengthFit", m_SignalDecayLengthFit);
      status = m_tuple->addItem("SignalKSDecayLengthErrorFit", m_SignalDecayLengthErrorFit);
      status = m_tuple->addItem("SignalKSChi2Fit", m_SignalChi2Fit);
      status = m_tuple->addItem("SignalKSpx", m_SignalKSpx);
      status = m_tuple->addItem("SignalKSpy", m_SignalKSpy);
      status = m_tuple->addItem("SignalKSpz", m_SignalKSpz);
      status = m_tuple->addItem("SignalKSenergy", m_SignalKSenergy);
      status = m_tuple->addItem("SignalKSPiPluspx", m_SignalKSPiPluspx);
      status = m_tuple->addItem("SignalKSPiPluspy", m_SignalKSPiPluspy);
      status = m_tuple->addItem("SignalKSPiPluspz", m_SignalKSPiPluspz);
      status = m_tuple->addItem("SignalKSPiPlusenergy", m_SignalKSPiPlusenergy);
      status = m_tuple->addItem("SignalKSPiMinuspx", m_SignalKSPiMinuspx);
      status = m_tuple->addItem("SignalKSPiMinuspy", m_SignalKSPiMinuspy);
      status = m_tuple->addItem("SignalKSPiMinuspz", m_SignalKSPiMinuspz);
      status = m_tuple->addItem("SignalKSPiMinusenergy", m_SignalKSPiMinusenergy);
      status = m_tuple->addItem("SignalKPluspx", m_SignalKPluspx);
      status = m_tuple->addItem("SignalKPluspy", m_SignalKPluspy);
      status = m_tuple->addItem("SignalKPluspz", m_SignalKPluspz);
      status = m_tuple->addItem("SignalKPlusenergy", m_SignalKPlusenergy);
      status = m_tuple->addItem("SignalKMinuspx", m_SignalKMinuspx);
      status = m_tuple->addItem("SignalKMinuspy", m_SignalKMinuspy);
      status = m_tuple->addItem("SignalKMinuspz", m_SignalKMinuspz);
      status = m_tuple->addItem("SignalKMinusenergy", m_SignalKMinusenergy);
      status = m_tuple->addItem("SignalKalmanFitSuccess", m_SignalKalmanFitSuccess);
      status = m_tuple->addItem("SignalKalmanFitChi2", m_SignalKalmanFitChi2);
      status = m_tuple->addItem("SignalKPluspxKalmanFit", m_SignalKPluspxKalmanFit);
      status = m_tuple->addItem("SignalKPluspyKalmanFit", m_SignalKPluspyKalmanFit);
      status = m_tuple->addItem("SignalKPluspzKalmanFit", m_SignalKPluspzKalmanFit);
      status = m_tuple->addItem("SignalKPlusenergyKalmanFit", m_SignalKPlusenergyKalmanFit);
      status = m_tuple->addItem("SignalKMinuspxKalmanFit", m_SignalKMinuspxKalmanFit);
      status = m_tuple->addItem("SignalKMinuspyKalmanFit", m_SignalKMinuspyKalmanFit);
      status = m_tuple->addItem("SignalKMinuspzKalmanFit", m_SignalKMinuspzKalmanFit);
      status = m_tuple->addItem("SignalKMinusenergyKalmanFit", m_SignalKMinusenergyKalmanFit);
      status = m_tuple->addItem("SignalKSpxKalmanFit", m_SignalKSpxKalmanFit);
      status = m_tuple->addItem("SignalKSpyKalmanFit", m_SignalKSpyKalmanFit);
      status = m_tuple->addItem("SignalKSpzKalmanFit", m_SignalKSpzKalmanFit);
      status = m_tuple->addItem("SignalKSenergyKalmanFit", m_SignalKSenergyKalmanFit);
      status = m_tuple->addItem("SignalIsSameDMother", m_SignalIsSameDMother);
      status = m_tuple->addItem("SignalPIDTrue", m_SignalPIDTrue);
      status = m_tuple->addItem("SignalKSPiPlusTrueID", m_SignalKSPiPlusTrueID);
      status = m_tuple->addItem("SignalKSPiMinusTrueID", m_SignalKSPiMinusTrueID);
      status = m_tuple->addItem("SignalKPlusTrueID", m_SignalKPlusTrueID);
      status = m_tuple->addItem("SignalKMinusTrueID", m_SignalKMinusTrueID);
      status = m_tuple->addItem("SignalKSPiPlusMotherTrueID", m_SignalKSPiPlusMotherTrueID);
      status = m_tuple->addItem("SignalKSPiMinusMotherTrueID", m_SignalKSPiMinusMotherTrueID);
      status = m_tuple->addItem("TagKpx", m_TagKpx);
      status = m_tuple->addItem("TagKpy", m_TagKpy);
      status = m_tuple->addItem("TagKpz", m_TagKpz);
      status = m_tuple->addItem("TagKenergy", m_TagKenergy);
      status = m_tuple->addItem("TagKCharge", m_TagKCharge);
      status = m_tuple->addItem("TagPi1px", m_TagPi1px);
      status = m_tuple->addItem("TagPi1py", m_TagPi1py);
      status = m_tuple->addItem("TagPi1pz", m_TagPi1pz);
      status = m_tuple->addItem("TagPi1energy", m_TagPi1energy);
      status = m_tuple->addItem("TagPi1Charge", m_TagPi1Charge);
      status = m_tuple->addItem("TagPi2px", m_TagPi2px);
      status = m_tuple->addItem("TagPi2py", m_TagPi2py);
      status = m_tuple->addItem("TagPi2pz", m_TagPi2pz);
      status = m_tuple->addItem("TagPi2energy", m_TagPi2energy);
      status = m_tuple->addItem("TagPi2Charge", m_TagPi2Charge);
      status = m_tuple->addItem("TagPi3px", m_TagPi3px);
      status = m_tuple->addItem("TagPi3py", m_TagPi3py);
      status = m_tuple->addItem("TagPi3pz", m_TagPi3pz);
      status = m_tuple->addItem("TagPi3energy", m_TagPi3energy);
      status = m_tuple->addItem("TagPi3Charge", m_TagPi3Charge);
      status = m_tuple->addItem("Tag12KSFitSuccess", m_Tag12KSFitSuccess);
      status = m_tuple->addItem("Tag12KSDecayLengthVeeVertex", m_Tag12DecayLengthVeeVertex);
      status = m_tuple->addItem("Tag12KSChi2VeeVertex", m_Tag12Chi2VeeVertex);
      status = m_tuple->addItem("Tag12KSMassVeeVertex", m_Tag12KSMassVeeVertex);
      status = m_tuple->addItem("Tag12KSDecayLengthFit", m_Tag12DecayLengthFit);
      status = m_tuple->addItem("Tag12KSDecayLengthErrorFit", m_Tag12DecayLengthErrorFit);
      status = m_tuple->addItem("Tag12KSChi2Fit", m_Tag12Chi2Fit);
      status = m_tuple->addItem("Tag13KSFitSuccess", m_Tag13KSFitSuccess);
      status = m_tuple->addItem("Tag13KSDecayLengthVeeVertex", m_Tag13DecayLengthVeeVertex);
      status = m_tuple->addItem("Tag13KSChi2VeeVertex", m_Tag13Chi2VeeVertex);
      status = m_tuple->addItem("Tag13KSMassVeeVertex", m_Tag13KSMassVeeVertex);
      status = m_tuple->addItem("Tag13KSDecayLengthFit", m_Tag13DecayLengthFit);
      status = m_tuple->addItem("Tag13KSDecayLengthErrorFit", m_Tag13DecayLengthErrorFit);
      status = m_tuple->addItem("Tag13KSChi2Fit", m_Tag13Chi2Fit);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKTrueID", m_TagKTrueID);
      status = m_tuple->addItem("TagPi1TrueID", m_TagPi1TrueID);
      status = m_tuple->addItem("TagPi2TrueID", m_TagPi2TrueID);
      status = m_tuple->addItem("TagPi3TrueID", m_TagPi3TrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSKK vs Kpipipi Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSKKVersusKpipipiDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSKK vs Kpipipi Double Tag Algorithm" << endreq;
  DTagTool DTTool;
  DTTool.setPID(true);
  if(DTTool.isDTagListEmpty()) {
    log << MSG::DEBUG << "No D candidates found" << endreq;
    return StatusCode::SUCCESS;
  }
  if(!DTTool.cosmicandleptonVeto()) {
    log << MSG::DEBUG << "Cosmic and lepton veto" << endreq;
    return StatusCode::SUCCESS;
  }
  if(DTTool.findDTag(EvtRecDTag::kD0toKsKK, EvtRecDTag::kD0toKPiPiPi)) {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
    DTagToolIterator DTTool_Signal_iter = DTTool.dtag1();
    DTagToolIterator DTTool_Tag_iter = DTTool.dtag2();
    StatusCode FillTupleStatus = FillTuple(DTTool_Signal_iter, DTTool_Tag_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      if(FillTupleStatus == StatusCode::RECOVERABLE) {
	log << MSG::WARNING << "Vertex fit of KS failed, skipping event" << endreq;
	return StatusCode::SUCCESS;
      }
      log << MSG::FATAL << "Assigning KSKK vs Kpipipi tuple info failed" << endreq;
      return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KSKKVersusKpipipiDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSKK vs Kpipipi Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSKKVersusKpipipiDoubleTag::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
  if(m_RunNumber < 0) {
    SmartDataPtr<Event::McParticleCol> MCParticleCol(eventSvc(), "/Event/MC/McParticleCol");
    if(!MCParticleCol) {
      return StatusCode::FAILURE;
    }
    IMcDecayModeSvc *IMcDecayModeService;
    StatusCode McDecayModeSVC_Status = service("McDecayModeSvc", IMcDecayModeService);
    if(McDecayModeSVC_Status.isFailure()) {
      return StatusCode::FAILURE;
    }
    FindMCInfo findMCInfo;
    StatusCode MCStatus = findMCInfo.CalculateMCInfo(MCParticleCol, IMcDecayModeService);
    if(MCStatus != StatusCode::SUCCESS) {
      return MCStatus;
    }
    m_NumberParticles = findMCInfo.GetNumberParticles();
    m_MCmode = findMCInfo.GetMCmode();
    for(int i = 0; i < m_NumberParticles; i++) {
      m_pdgID[i] = findMCInfo.GetpdgID(i);
      m_MotherIndex[i] = findMCInfo.GetMotherIndex(i);
      m_TruePx[i] = findMCInfo.GetTruePx(i);
      m_TruePy[i] = findMCInfo.GetTruePy(i);
      m_TruePz[i] = findMCInfo.GetTruePz(i);
      m_TrueEnergy[i] = findMCInfo.GetTrueEnergy(i);
    }
    m_NumberParticlesStripped = findMCInfo.GetNumberParticlesStripped();
    for(int i = 0; i < m_NumberParticlesStripped; i++) {
      m_pdgIDStripped[i] = findMCInfo.GetpdgIDStripped(i);
      m_MotherIndexStripped[i] = findMCInfo.GetMotherIndexStripped(i);
    }
  }
  m_SignalDMass = (*DTTool_Signal_iter)->mass();
  m_SignalMBC = (*DTTool_Signal_iter)->mBC();
  m_SignalDeltaE = (*DTTool_Signal_iter)->deltaE();
  m_SignalBeamE = (*DTTool_Signal_iter)->beamE();
  m_TagDMass = (*DTTool_Tag_iter)->mass();
  m_TagMBC = (*DTTool_Tag_iter)->mBC();
  m_TagDeltaE = (*DTTool_Tag_iter)->deltaE();
  m_TagBeamE = (*DTTool_Tag_iter)->beamE();
  m_SignalDpx = (*DTTool_Signal_iter)->p4().x();
  m_SignalDpy = (*DTTool_Signal_iter)->p4().y();
  m_SignalDpz = (*DTTool_Signal_iter)->p4().z();
  m_SignalDenergy = (*DTTool_Signal_iter)->p4().t();
  m_TagDpx = (*DTTool_Tag_iter)->p4().x();
  m_TagDpy = (*DTTool_Tag_iter)->p4().y();
  m_TagDpz = (*DTTool_Tag_iter)->p4().z();
  m_TagDenergy = (*DTTool_Tag_iter)->p4().t();
  FindKSKKTagInfo findKSKKTagInfo;
  StatusCode status = findKSKKTagInfo.CalculateTagInfo(DTTool_Signal_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_SignalDecayLengthVeeVertex = findKSKKTagInfo.GetDecayLengthVeeVertex();
  m_SignalChi2VeeVertex = findKSKKTagInfo.GetChi2VeeVertex();
  m_SignalKSMassVeeVertex = findKSKKTagInfo.GetKSMassVeeVertex();
  m_SignalDecayLengthFit = findKSKKTagInfo.GetDecayLengthFit();
  m_SignalDecayLengthErrorFit = findKSKKTagInfo.GetDecayLengthErrorFit();
  m_SignalChi2Fit = findKSKKTagInfo.GetChi2Fit();
  m_SignalKSPiPluspx = findKSKKTagInfo.GetKSPiPlusP(0);
  m_SignalKSPiPluspy = findKSKKTagInfo.GetKSPiPlusP(1);
  m_SignalKSPiPluspz = findKSKKTagInfo.GetKSPiPlusP(2);
  m_SignalKSPiPlusenergy = findKSKKTagInfo.GetKSPiPlusP(3);
  m_SignalKSPiMinuspx = findKSKKTagInfo.GetKSPiMinusP(0);
  m_SignalKSPiMinuspy = findKSKKTagInfo.GetKSPiMinusP(1);
  m_SignalKSPiMinuspz = findKSKKTagInfo.GetKSPiMinusP(2);
  m_SignalKSPiMinusenergy = findKSKKTagInfo.GetKSPiMinusP(3);
  m_SignalKSpx = findKSKKTagInfo.GetKShortP(0);
  m_SignalKSpy = findKSKKTagInfo.GetKShortP(1);
  m_SignalKSpz = findKSKKTagInfo.GetKShortP(2);
  m_SignalKSenergy = findKSKKTagInfo.GetKShortP(3);
  m_SignalKPluspx = findKSKKTagInfo.GetKPlusP(0);
  m_SignalKPluspy = findKSKKTagInfo.GetKPlusP(1);
  m_SignalKPluspz = findKSKKTagInfo.GetKPlusP(2);
  m_SignalKPlusenergy = findKSKKTagInfo.GetKPlusP(3);
  m_SignalKMinuspx = findKSKKTagInfo.GetKMinusP(0);
  m_SignalKMinuspy = findKSKKTagInfo.GetKMinusP(1);
  m_SignalKMinuspz = findKSKKTagInfo.GetKMinusP(2);
  m_SignalKMinusenergy = findKSKKTagInfo.GetKMinusP(3);
  m_SignalKalmanFitSuccess = findKSKKTagInfo.GetKalmanFitSuccess();
  m_SignalKalmanFitChi2 = findKSKKTagInfo.GetKalmanFitChi2();
  m_SignalKSpxKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(0);
  m_SignalKSpyKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(1);
  m_SignalKSpzKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(2);
  m_SignalKSenergyKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(3);
  m_SignalKPluspxKalmanFit = findKSKKTagInfo.GetKPlusPKalmanFit(0);
  m_SignalKPluspyKalmanFit = findKSKKTagInfo.GetKPlusPKalmanFit(1);
  m_SignalKPluspzKalmanFit = findKSKKTagInfo.GetKPlusPKalmanFit(2);
  m_SignalKPlusenergyKalmanFit = findKSKKTagInfo.GetKPlusPKalmanFit(3);
  m_SignalKMinuspxKalmanFit = findKSKKTagInfo.GetKMinusPKalmanFit(0);
  m_SignalKMinuspyKalmanFit = findKSKKTagInfo.GetKMinusPKalmanFit(1);
  m_SignalKMinuspzKalmanFit = findKSKKTagInfo.GetKMinusPKalmanFit(2);
  m_SignalKMinusenergyKalmanFit = findKSKKTagInfo.GetKMinusPKalmanFit(3);
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKSKKTagInfo.GetDaughterTrackID();
    PIDTruth PID_Truth(DaughterTrackIDs, 4, this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 321, -321};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKSPiPlusTrueID = ReconstructedPID[0];
    m_SignalKSPiMinusTrueID = ReconstructedPID[1];
    m_SignalKPlusTrueID = ReconstructedPID[2];
    m_SignalKMinusTrueID = ReconstructedPID[3];
    m_SignalKSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_SignalKSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
  }
  FindKpipipiTagInfo findKpipipiTagInfo;
  status = findKpipipiTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_TagKpx = findKpipipiTagInfo.GetKP(0);
  m_TagKpy = findKpipipiTagInfo.GetKP(1);
  m_TagKpz = findKpipipiTagInfo.GetKP(2);
  m_TagKenergy = findKpipipiTagInfo.GetKP(3);
  m_TagKCharge = findKpipipiTagInfo.GetKCharge();
  m_TagPi1px = findKpipipiTagInfo.GetPi1P(0);
  m_TagPi1py = findKpipipiTagInfo.GetPi1P(1);
  m_TagPi1pz = findKpipipiTagInfo.GetPi1P(2);
  m_TagPi1energy = findKpipipiTagInfo.GetPi1P(3);
  m_TagPi1Charge = findKpipipiTagInfo.GetPi1Charge();
  m_TagPi2px = findKpipipiTagInfo.GetPi2P(0);
  m_TagPi2py = findKpipipiTagInfo.GetPi2P(1);
  m_TagPi2pz = findKpipipiTagInfo.GetPi2P(2);
  m_TagPi2energy = findKpipipiTagInfo.GetPi2P(3);
  m_TagPi2Charge = findKpipipiTagInfo.GetPi2Charge();
  m_TagPi3px = findKpipipiTagInfo.GetPi3P(0);
  m_TagPi3py = findKpipipiTagInfo.GetPi3P(1);
  m_TagPi3pz = findKpipipiTagInfo.GetPi3P(2);
  m_TagPi3energy = findKpipipiTagInfo.GetPi3P(3);
  m_TagPi3Charge = findKpipipiTagInfo.GetPi3Charge();
  m_Tag12KSFitSuccess = findKpipipiTagInfo.GetKSFitSuccess12();
  m_Tag12DecayLengthVeeVertex = findKpipipiTagInfo.GetDecayLengthVeeVertex12();
  m_Tag12Chi2VeeVertex = findKpipipiTagInfo.GetChi2VeeVertex12();
  m_Tag12KSMassVeeVertex = findKpipipiTagInfo.GetKSMassVeeVertex12();
  m_Tag12DecayLengthFit = findKpipipiTagInfo.GetDecayLengthFit12();
  m_Tag12DecayLengthErrorFit = findKpipipiTagInfo.GetDecayLengthErrorFit12();
  m_Tag12Chi2Fit = findKpipipiTagInfo.GetChi2Fit12();
  m_Tag13KSFitSuccess = findKpipipiTagInfo.GetKSFitSuccess13();
  m_Tag13DecayLengthVeeVertex = findKpipipiTagInfo.GetDecayLengthVeeVertex13();
  m_Tag13Chi2VeeVertex = findKpipipiTagInfo.GetChi2VeeVertex13();
  m_Tag13KSMassVeeVertex = findKpipipiTagInfo.GetKSMassVeeVertex13();
  m_Tag13DecayLengthFit = findKpipipiTagInfo.GetDecayLengthFit13();
  m_Tag13DecayLengthErrorFit = findKpipipiTagInfo.GetDecayLengthErrorFit13();
  m_Tag13Chi2Fit = findKpipipiTagInfo.GetChi2Fit13();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKpipipiTagInfo.GetDaughterTrackID(), 4, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {321*m_TagKCharge, 211*m_TagPi1Charge, 211*m_TagPi2Charge, 211*m_TagPi3Charge};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKTrueID = ReconstructedPID[0];
    m_TagPi1TrueID = ReconstructedPID[1];
    m_TagPi2TrueID = ReconstructedPID[2];
    m_TagPi3TrueID = ReconstructedPID[3];
  }
  return StatusCode::SUCCESS;
}
