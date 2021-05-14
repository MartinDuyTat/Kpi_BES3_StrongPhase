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
}

KpiStrongPhase::~KpiStrongPhase() {
}

StatusCode KpiStrongPhase::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Creating KpiStrongPhase Double Tag Algorithm" << endreq;
  StatusCode sc;
  if(m_recKSKKVersusKpiTag) {
    sc = createSubAlgorithm("KSKKVersusKpiDoubleTag", "KSKKVersusKpiDoubleTag", m_KKTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSKKVersusKpiDoubleTag" << endreq;
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
