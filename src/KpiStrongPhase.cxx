// Martin Duy Tat 12th February 2021

// KpiStrongPhase
#include "KpiStrongPhase/KpiStrongPhase.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/StatusCode.h"
// STL
#include<string>

KpiStrongPhase::KpiStrongPhase(const std::string& name, ISvcLocator* pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("KSKKVersusKpiDoubleTag", m_recKSKKVersusKpiTag = true);
  declareProperty("KSKKVersusKpipi0DoubleTag", m_recKSKKVersusKpipi0Tag = true);
  declareProperty("KSKKVersusKpipipiDoubleTag", m_recKSKKVersusKpipipiTag = true);
  declareProperty("KSKKVersusKeNuDoubleTag", m_recKSKKVersusKeNuTag = true);
  declareProperty("KLKKVersusKpiDoubleTag", m_recKLKKVersusKpiTag = true);
  declareProperty("KLKKVersusKpipi0DoubleTag", m_recKLKKVersusKpipi0Tag = true);
  declareProperty("KLKKVersusKpipipiDoubleTag", m_recKLKKVersusKpipipiTag = true);
  declareProperty("KSpipiVersusKpiDoubleTag", m_recKSpipiVersusKpiTag = true);
  declareProperty("KSpipiVersusKeNuDoubleTag", m_recKSpipiVersusKeNuTag = true);
}

KpiStrongPhase::~KpiStrongPhase() {
}

StatusCode KpiStrongPhase::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Creating KpiStrongPhase Double Tag Algorithm" << endreq;
  StatusCode sc;
  if(m_recKSKKVersusKpiTag) {
    sc = createSubAlgorithm("KSKKVersusKpiDoubleTag", "KSKKVersusKpiDoubleTag", m_KSKKVersusKpiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSKKVersusKpiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSKKVersusKpipi0Tag) {
    sc = createSubAlgorithm("KSKKVersusKpipi0DoubleTag", "KSKKVersusKpipi0DoubleTag", m_KSKKVersusKpipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSKKVersusKpipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSKKVersusKpiTag) {
    sc = createSubAlgorithm("KSKKVersusKpipipiDoubleTag", "KSKKVersusKpipipiDoubleTag", m_KSKKVersusKpipipiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSKKVersusKpipipiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSKKVersusKeNuTag) {
    sc = createSubAlgorithm("KSKKVersusKeNuDoubleTag", "KSKKVersusKeNuDoubleTag", m_KSKKVersusKeNuTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSKKVersusKeNuDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKLKKVersusKpiTag) {
    sc = createSubAlgorithm("KLKKVersusKpiDoubleTag", "KLKKVersusKpiDoubleTag", m_KLKKVersusKpiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KLKKVersusKpiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKLKKVersusKpipi0Tag) {
    sc = createSubAlgorithm("KLKKVersusKpipi0DoubleTag", "KLKKVersusKpipi0DoubleTag", m_KLKKVersusKpipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KLKKVersusKpipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKLKKVersusKpiTag) {
    sc = createSubAlgorithm("KLKKVersusKpipipiDoubleTag", "KLKKVersusKpipipiDoubleTag", m_KLKKVersusKpipipiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KLKKVersusKpipipiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpipiVersusKpiTag) {
    sc = createSubAlgorithm("KSpipiVersusKpiDoubleTag", "KSpipiVersusKpiDoubleTag", m_KSpipiVersusKpiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSpipiVersusKpiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpipiVersusKeNuTag) {
    sc = createSubAlgorithm("KSpipiVersusKeNuDoubleTag", "KSpipiVersusKeNuDoubleTag", m_KSpipiVersusKeNuTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSpipiVersusKeNuDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
}

StatusCode KpiStrongPhase::execute() {
  for(std::vector<Algorithm*>::const_iterator it = subAlgorithms()->begin(); it != subAlgorithms()->end(); it++) {
    StatusCode sc = (*it)->execute();
    if(sc.isFailure()) {
      return StatusCode::FAILURE;
    }
  }
}

StatusCode KpiStrongPhase::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "KpiStrongPhase Double Tag Algorithm finalized" << endmsg;
  return StatusCode::SUCCESS;
}
