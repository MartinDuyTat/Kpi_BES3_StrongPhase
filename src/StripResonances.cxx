// Martin Duy Tat 23rd March 2021

// KpiStrongPhase
#include "KpiStrongPhase/StripResonances.h"
// ROOT
#include "TMath.h"
// STL
#include <vector>
#include <algorithm>

StripResonances::StripResonances(const std::vector<int> &KeepList) {
  int FinalParticles[15] = {30443, 421, 130, 310, 321, 111, 221, 331, 211, 11, 12, 13, 14, 22, 2212};
  m_KeepList = std::vector<int>(FinalParticles, FinalParticles + 15);
  m_KeepList.insert(m_KeepList.end(), KeepList.begin(), KeepList.end());
}

int StripResonances::FindResonanceIndex(const std::vector<int> &pidID) const {
  for(unsigned int i = 0; i < pidID.size(); i++) {
    if(std::find(m_KeepList.begin(), m_KeepList.end(), pidID[i]) == m_KeepList.end() && std::find(m_KeepList.begin(), m_KeepList.end(), -pidID[i]) == m_KeepList.end()) {
      return i;
    }
  }
  return -1;
}

void StripResonances::RemoveResonance(int ResonanceIndex, std::vector<int> &pidID, std::vector<int> &MotherIndex) const {
  int Mother = MotherIndex[ResonanceIndex];
  for(unsigned int i = 0; i < pidID.size(); i++) {
    if(MotherIndex[i] == ResonanceIndex) {
      MotherIndex[i] = Mother;
    } else if (MotherIndex[i] > ResonanceIndex) {
      MotherIndex[i]--;
    }
  }
  pidID.erase(pidID.begin() + ResonanceIndex);
  MotherIndex.erase(MotherIndex.begin() + ResonanceIndex);
}

void StripResonances::RemoveIntermediateResonances(std::vector<int> &pidID, std::vector<int> &MotherIndex) const {
  int ResonanceIndex = FindResonanceIndex(pidID);
  while(ResonanceIndex != -1) {
    RemoveResonance(ResonanceIndex, pidID, MotherIndex);
    ResonanceIndex = FindResonanceIndex(pidID);
  }
}
