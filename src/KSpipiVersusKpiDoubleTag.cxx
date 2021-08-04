// Martin Duy Tat 4th August 2021

// KpiStrongPhase
#include "KpiStrongPhase/KSpipiVersusKpiDoubleTag.h"
#include "KpiStrongPhase/FindKSpipiTagInfo.h"
#include "KpiStrongPhase/FindKpiTagInfo.h"
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

KSpipiVersusKpiDoubleTag::KSpipiVersusKpiDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSpipiVersusKpiDoubleTag::~KSpipiVersusKpiDoubleTag() {
}

StatusCode KSpipiVersusKpiDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSpipi vs Kpi Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KPISTRONGPHASE/KSpipiVersusKpiDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KPISTRONGPHASE/KSpipiVersusKpiDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KSpipi vs D->Kpi events");
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
      status = m_tuple->addItem("SignalPiPluspx", m_SignalPiPluspx);
      status = m_tuple->addItem("SignalPiPluspy", m_SignalPiPluspy);
      status = m_tuple->addItem("SignalPiPluspz", m_SignalPiPluspz);
      status = m_tuple->addItem("SignalPiPlusenergy", m_SignalPiPlusenergy);
      status = m_tuple->addItem("SignalPiMinuspx", m_SignalPiMinuspx);
      status = m_tuple->addItem("SignalPiMinuspy", m_SignalPiMinuspy);
      status = m_tuple->addItem("SignalPiMinuspz", m_SignalPiMinuspz);
      status = m_tuple->addItem("SignalPiMinusenergy", m_SignalPiMinusenergy);
      status = m_tuple->addItem("SignalKalmanFitSuccess", m_SignalKalmanFitSuccess);
      status = m_tuple->addItem("SignalKalmanFitChi2", m_SignalKalmanFitChi2);
      status = m_tuple->addItem("SignalPiPluspxKalmanFit", m_SignalPiPluspxKalmanFit);
      status = m_tuple->addItem("SignalPiPluspyKalmanFit", m_SignalPiPluspyKalmanFit);
      status = m_tuple->addItem("SignalPiPluspzKalmanFit", m_SignalPiPluspzKalmanFit);
      status = m_tuple->addItem("SignalPiPlusenergyKalmanFit", m_SignalPiPlusenergyKalmanFit);
      status = m_tuple->addItem("SignalPiMinuspxKalmanFit", m_SignalPiMinuspxKalmanFit);
      status = m_tuple->addItem("SignalPiMinuspyKalmanFit", m_SignalPiMinuspyKalmanFit);
      status = m_tuple->addItem("SignalPiMinuspzKalmanFit", m_SignalPiMinuspzKalmanFit);
      status = m_tuple->addItem("SignalPiMinusenergyKalmanFit", m_SignalPiMinusenergyKalmanFit);
      status = m_tuple->addItem("SignalKSpxKalmanFit", m_SignalKSpxKalmanFit);
      status = m_tuple->addItem("SignalKSpyKalmanFit", m_SignalKSpyKalmanFit);
      status = m_tuple->addItem("SignalKSpzKalmanFit", m_SignalKSpzKalmanFit);
      status = m_tuple->addItem("SignalKSenergyKalmanFit", m_SignalKSenergyKalmanFit);
      status = m_tuple->addItem("SignalIsSameDMother", m_SignalIsSameDMother);
      status = m_tuple->addItem("SignalPIDTrue", m_SignalPIDTrue);
      status = m_tuple->addItem("SignalKSPiPlusTrueID", m_SignalKSPiPlusTrueID);
      status = m_tuple->addItem("SignalKSPiMinusTrueID", m_SignalKSPiMinusTrueID);
      status = m_tuple->addItem("SignalPiPlusTrueID", m_SignalPiPlusTrueID);
      status = m_tuple->addItem("SignalPiMinusTrueID", m_SignalPiMinusTrueID);
      status = m_tuple->addItem("SignalKSPiPlusMotherTrueID", m_SignalKSPiPlusMotherTrueID);
      status = m_tuple->addItem("SignalKSPiMinusMotherTrueID", m_SignalKSPiMinusMotherTrueID);
      status = m_tuple->addItem("TagPipx", m_TagPipx);
      status = m_tuple->addItem("TagPipy", m_TagPipy);
      status = m_tuple->addItem("TagPipz", m_TagPipz);
      status = m_tuple->addItem("TagPienergy", m_TagPienergy);
      status = m_tuple->addItem("TagPiCharge", m_TagPiCharge);
      status = m_tuple->addItem("TagKpx", m_TagKpx);
      status = m_tuple->addItem("TagKpy", m_TagKpy);
      status = m_tuple->addItem("TagKpz", m_TagKpz);
      status = m_tuple->addItem("TagKenergy", m_TagKenergy);
      status = m_tuple->addItem("TagKCharge", m_TagKCharge);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKTrueID", m_TagKTrueID);
      status = m_tuple->addItem("TagPiTrueID", m_TagPiTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSpipi vs Kpi Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSpipiVersusKpiDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSpipi vs Kpi Double Tag Algorithm" << endreq;
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
  if(DTTool.findDTag(EvtRecDTag::kD0toKsPiPi, EvtRecDTag::kD0toKPi)) {
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
      log << MSG::FATAL << "Assigning KSpipi vs Kpi tuple info failed" << endreq;
      return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KSpipiVersusKpiDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSpipi vs Kpi Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSpipiVersusKpiDoubleTag::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
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
  FindKSpipiTagInfo findKSpipiTagInfo;
  StatusCode status = findKSpipiTagInfo.CalculateTagInfo(DTTool_Signal_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_SignalDecayLengthVeeVertex = findKSpipiTagInfo.GetDecayLengthVeeVertex();
  m_SignalChi2VeeVertex = findKSpipiTagInfo.GetChi2VeeVertex();
  m_SignalKSMassVeeVertex = findKSpipiTagInfo.GetKSMassVeeVertex();
  m_SignalDecayLengthFit = findKSpipiTagInfo.GetDecayLengthFit();
  m_SignalDecayLengthErrorFit = findKSpipiTagInfo.GetDecayLengthErrorFit();
  m_SignalChi2Fit = findKSpipiTagInfo.GetChi2Fit();
  m_SignalKSPiPluspx = findKSpipiTagInfo.GetKSPiPlusP(0);
  m_SignalKSPiPluspy = findKSpipiTagInfo.GetKSPiPlusP(1);
  m_SignalKSPiPluspz = findKSpipiTagInfo.GetKSPiPlusP(2);
  m_SignalKSPiPlusenergy = findKSpipiTagInfo.GetKSPiPlusP(3);
  m_SignalKSPiMinuspx = findKSpipiTagInfo.GetKSPiMinusP(0);
  m_SignalKSPiMinuspy = findKSpipiTagInfo.GetKSPiMinusP(1);
  m_SignalKSPiMinuspz = findKSpipiTagInfo.GetKSPiMinusP(2);
  m_SignalKSPiMinusenergy = findKSpipiTagInfo.GetKSPiMinusP(3);
  m_SignalKSpx = findKSpipiTagInfo.GetKShortP(0);
  m_SignalKSpy = findKSpipiTagInfo.GetKShortP(1);
  m_SignalKSpz = findKSpipiTagInfo.GetKShortP(2);
  m_SignalKSenergy = findKSpipiTagInfo.GetKShortP(3);
  m_SignalPiPluspx = findKSpipiTagInfo.GetPiPlusP(0);
  m_SignalPiPluspy = findKSpipiTagInfo.GetPiPlusP(1);
  m_SignalPiPluspz = findKSpipiTagInfo.GetPiPlusP(2);
  m_SignalPiPlusenergy = findKSpipiTagInfo.GetPiPlusP(3);
  m_SignalPiMinuspx = findKSpipiTagInfo.GetPiMinusP(0);
  m_SignalPiMinuspy = findKSpipiTagInfo.GetPiMinusP(1);
  m_SignalPiMinuspz = findKSpipiTagInfo.GetPiMinusP(2);
  m_SignalPiMinusenergy = findKSpipiTagInfo.GetPiMinusP(3);
  m_SignalKalmanFitSuccess = findKSpipiTagInfo.GetKalmanFitSuccess();
  m_SignalKalmanFitChi2 = findKSpipiTagInfo.GetKalmanFitChi2();
  m_SignalKSpxKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(0);
  m_SignalKSpyKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(1);
  m_SignalKSpzKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(2);
  m_SignalKSenergyKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(3);
  m_SignalPiPluspxKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(0);
  m_SignalPiPluspyKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(1);
  m_SignalPiPluspzKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(2);
  m_SignalPiPlusenergyKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(3);
  m_SignalPiMinuspxKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(0);
  m_SignalPiMinuspyKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(1);
  m_SignalPiMinuspzKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(2);
  m_SignalPiMinusenergyKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(3);
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKSpipiTagInfo.GetDaughterTrackID();
    PIDTruth PID_Truth(DaughterTrackIDs, 4, this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKSPiPlusTrueID = ReconstructedPID[0];
    m_SignalKSPiMinusTrueID = ReconstructedPID[1];
    m_SignalPiPlusTrueID = ReconstructedPID[2];
    m_SignalPiMinusTrueID = ReconstructedPID[3];
    m_SignalKSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_SignalKSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
  }
  FindKpiTagInfo findKpiTagInfo;
  status = findKpiTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_TagPipx = findKpiTagInfo.GetPiP(0);
  m_TagPipy = findKpiTagInfo.GetPiP(1);
  m_TagPipz = findKpiTagInfo.GetPiP(2);
  m_TagPienergy = findKpiTagInfo.GetPiP(3);
  m_TagPiCharge = findKpiTagInfo.GetPiCharge();
  m_TagKpx = findKpiTagInfo.GetKP(0);
  m_TagKpy = findKpiTagInfo.GetKP(1);
  m_TagKpz = findKpiTagInfo.GetKP(2);
  m_TagKenergy = findKpiTagInfo.GetKP(3);
  m_TagKCharge = findKpiTagInfo.GetKCharge();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKpiTagInfo.GetDaughterTrackID(), 2, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[2] = {321*m_TagKCharge, 211*m_TagPiCharge};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 2);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKTrueID = ReconstructedPID[0];
    m_TagPiTrueID = ReconstructedPID[1];
  }
  return StatusCode::SUCCESS;
}
