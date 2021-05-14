// Martin Duy Tat 23rd March 2021

// KpiStrongPhase
#include "KpiStrongPhase/PIDTruth.h"
// Gaudi
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"
// BOSS
#include "EventNavigator/EventNavigator.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "McTruth/McParticle.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// STL
#include <vector>
#include <algorithm>

PIDTruth::PIDTruth(const std::vector<int> &TrackID, int NumberCharged, const Algorithm *algorithm): m_TrackID(TrackID), m_NumberCharged(NumberCharged), m_algorithm(algorithm) {
}

int PIDTruth::MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const {
  SmartDataPtr < EvtRecEvent > evtRecEvent(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(m_algorithm->eventSvc(), "/Event/Navigator");
  int ismatched = 0;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG == -1 && (GParPDG == 0 || GParPDG == thebestGParent)) {
        return thebestParent;
      }
      if (mcParPDG == -1 && GParPDG == 0) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }
      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }
  return ismatched;
}

int PIDTruth::MCSHPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const {
  SmartDataPtr < EvtRecEvent > evtRecEvent(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(m_algorithm->eventSvc(), "/Event/Navigator");
  int ismatched = -1;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }

      if (mcParPDG == -1) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }

      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
          thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }
  return ismatched;
}

int PIDTruth::FindDOrigin(int TrackID, bool Charged) const {
  if(Charged) {
    if(MCTKPIDCHG(TrackID, 0, 0, 421) == 1) {
      return 421;
    } else if(MCTKPIDCHG(TrackID, 0, 0, -421) == 1) {
      return -421;
    } else {
      return 0;
    }
  } else {
    if(MCSHPIDCHG(TrackID, 0, 0, 421) == 1) {
      return 421;
    } else if(MCSHPIDCHG(TrackID, 0, 0, -421) == 1) {
      return -421;
    } else {
      return 0;
    }
  }
}

bool PIDTruth::SameDMother(const std::vector<int> &IgnoreTrackID) const {
  std::vector<int> IsD0Mother, IsD0barMother;
  for(std::vector<int>::const_iterator iter = m_TrackID.begin(); iter != m_TrackID.end(); iter++) {
    // If track ID is on ignore list, skip
    if(IgnoreTrackID.size() != 0 && std::find(IgnoreTrackID.begin(), IgnoreTrackID.end(), *iter) != IgnoreTrackID.end()) {
      continue;
    }
    // Check if every track ID is a daughter originating from D0 or D0bar
    if(iter - m_TrackID.begin() < m_NumberCharged) {
      IsD0Mother.push_back(MCTKPIDCHG(*iter, 0, 0, 421));
      IsD0barMother.push_back(MCTKPIDCHG(*iter, 0, 0, -421));
    } else {
      IsD0Mother.push_back(MCSHPIDCHG(*iter, 0, 0, 421));
      IsD0barMother.push_back(MCSHPIDCHG(*iter, 0, 0, -421));
    }
  }
  bool isD0Mother = true, isD0barMother = true;
  // If any daughter doesn't originate from D0, set isD0Mother to false
  for(std::vector<int>::iterator iter = IsD0Mother.begin(); iter != IsD0Mother.end(); iter++) {
    if(*iter != 1) {
      isD0Mother = false;
      break;
    }
  }
  // If any daughter doesn't originate from D0bar, set isD0barMother to false
  for(std::vector<int>::iterator iter = IsD0barMother.begin(); iter != IsD0barMother.end(); iter++) {
    if(*iter != 1) {
      isD0barMother = false;
      break;
    }
  }
  // Return true of all daughters originate from D0, or if all daughters originate from D0bar, otherwise return false
  return isD0Mother || isD0barMother;
}

bool PIDTruth::FindTrueID(std::vector<int> &ParticleID) const {
  std::vector<int> TrueID;
  bool PIDMatch = true;
  for(unsigned int i = 0; i < m_TrackID.size(); i++) {
    if(i < m_NumberCharged) {
      TrueID.push_back(MCTKPIDCHG(m_TrackID[i], -1, 0, 0));
    } else {
      TrueID.push_back(MCSHPIDCHG(m_TrackID[i], -1, 0, 0));
    }
    if(PIDMatch && ParticleID[i] != TrueID[i] && ParticleID[i] != 0) {
      PIDMatch = false;
    }
  }
  std::swap(ParticleID, TrueID);
  return PIDMatch;
}

int PIDTruth::GetTrueMotherID(int TrackID, bool Charged) const {
  if(Charged) {
    return MCTKPIDCHG(TrackID, 0, -1, 0);
  } else {
    return MCSHPIDCHG(TrackID, 0, -1, 0);
  }
}
