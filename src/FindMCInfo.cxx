// Martin Duy Tat 11th February 2021, based on code by Yu Zhang

// KpiStrongPhase
#include "KpiStrongPhase/FindMCInfo.h"
#include "KpiStrongPhase/StripResonances.h"
// Gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
// Boss
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>

FindMCInfo::FindMCInfo(): m_NumberParticles(0), m_MCmode(0) {
}

FindMCInfo::~FindMCInfo() {
}

StatusCode FindMCInfo::CalculateMCInfo(SmartDataPtr<Event::McParticleCol> MCParticleCol, IMcDecayModeSvc *IMcDecayModeService, int KeepID) {
  std::vector<int> pdgID, MotherIndex;
  std::vector<double> TruePx, TruePy, TruePz, TrueEnergy;
  for(Event::McParticleCol::iterator MCParticleCol_iter = MCParticleCol->begin(); MCParticleCol_iter != MCParticleCol->end(); MCParticleCol_iter++) {
    if((*MCParticleCol_iter)->primaryParticle() || !(*MCParticleCol_iter)->decayFromGenerator()) {
      continue;
    }
    if((*MCParticleCol_iter)->particleProperty() == 30443) {      
      McDecayModeSvc *McDecayModeService = dynamic_cast<McDecayModeSvc*>(IMcDecayModeService);
      m_MCmode = McDecayModeService->extract(*MCParticleCol_iter, pdgID, MotherIndex, TruePx, TruePy, TruePz, TrueEnergy);
    }
  }
  m_NumberParticles = pdgID.size();
  for(int i = 0; i < m_NumberParticles; i++) {
    m_pdgID.push_back(pdgID[i]);
    m_MotherIndex.push_back(MotherIndex[i]);
    m_TruePx.push_back(TruePx[i]);
    m_TruePy.push_back(TruePy[i]);
    m_TruePz.push_back(TruePz[i]);
    m_TrueEnergy.push_back(TrueEnergy[i]);
  }
  StripResonances stripResonances;
  if(KeepID != 0) {
    std::vector<int> IDToKeep;
    IDToKeep.push_back(KeepID);
    stripResonances = StripResonances(IDToKeep);
  }
  m_pdgIDStripped = m_pdgID;
  m_MotherIndexStripped = m_MotherIndex;
  stripResonances.RemoveIntermediateResonances(m_pdgIDStripped, m_MotherIndexStripped);
  m_NumberParticlesStripped = m_pdgIDStripped.size();
  return StatusCode::SUCCESS;
}

int FindMCInfo::GetNumberParticles() {
  return m_NumberParticles;
}

int FindMCInfo::GetpdgID(int i) {
  return m_pdgID[i];
}

int FindMCInfo::GetMotherIndex(int i) {
  return m_MotherIndex[i];
}

int FindMCInfo::GetNumberParticlesStripped() {
  return m_NumberParticlesStripped;
}

int FindMCInfo::GetpdgIDStripped(int i) {
  return m_pdgIDStripped[i];
}

int FindMCInfo::GetMotherIndexStripped(int i) {
  return m_MotherIndexStripped[i];
}

int FindMCInfo::GetMCmode() {
  return m_MCmode;
}

double FindMCInfo::GetTruePx(int i) {
  return m_TruePx[i];
}

double FindMCInfo::GetTruePy(int i) {
  return m_TruePy[i];
}

double FindMCInfo::GetTruePz(int i) {
  return m_TruePz[i];
}

double FindMCInfo::GetTrueEnergy(int i) {
  return m_TrueEnergy[i];
}
