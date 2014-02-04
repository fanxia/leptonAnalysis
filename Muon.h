#include "math.h"
#include "TMath.h"

#include <vector>

#include "../src/SusyEvent.h"

using namespace std;

bool isTightMuon(susy::Muon mu, vector<susy::Track> tracks, double d0, double dz) {

  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  bool hasTracks = (int)mu.trackIndex < (int)tracks.size() && 
    (int)mu.standAloneTrackIndex < (int)tracks.size() && 
    (int)mu.combinedTrackIndex < (int)tracks.size() && 
    (int)mu.bestTrackIndex() < (int)tracks.size() && 
    (int)mu.bestTrackIndex() >= 0;

  if(!hasTracks) return false;
  
  bool passes = mu.isGlobalMuon() && 
    mu.isPFMuon() && 
    tracks[mu.combinedTrackIndex].normChi2() < 10. && 
    mu.nValidMuonHits > 0 && 
    mu.nMatchedStations > 1 &&
    fabs(d0) < 0.2 &&
    fabs(dz) < 0.5 &&
    tracks[mu.trackIndex].numberOfValidPixelHits > 0 && 
    (mu.nPixelLayersWithMeasurement + mu.nStripLayersWithMeasurement) > 5 && 
    mu.momentum.Pt() > 30. && // STUDY THIS ONE (ttH(bb) for now) -- ttH(gg) uses 20
    fabs(mu.momentum.Eta()) < 2.1 && // STUDY THIS ONE (ttH(bb) for now) -- ttH(gg) uses 2.4
    mu_iso / mu.momentum.Pt() < 0.12; // STUDY THIS ONE (ttH(bb) for now) -- ttH(gg) uses 0.2
  
  return passes;
  
}

bool isVetoMuon(susy::Muon mu) {

  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  bool passes = (mu.isPFMuon() &&
		 (mu.isGlobalMuon() || mu.isTrackerMuon()) &&
		 mu.momentum.Pt() > 10. && // (ttH(bb) for now)
		 fabs(mu.momentum.Eta()) < 2.5 &&
		 mu_iso / mu.momentum.Pt() < 0.2);

  return passes;
}

bool isJetVetoMuon(susy::Muon mu, vector<susy::Track> tracks, double d0, double dz) {

  bool hasTracks = (int)mu.trackIndex < (int)tracks.size() &&
    (int)mu.combinedTrackIndex < (int)tracks.size()
    && (int)mu.combinedTrackIndex >= 0;
  if(!hasTracks) return false;

  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  if(mu.isGlobalMuon() &&
     mu.isPFMuon() &&
     fabs(mu.momentum.Eta()) < 2.6 &&
     mu.momentum.Pt() >= 15.0 &&
     tracks[mu.combinedTrackIndex].normChi2() < 10. &&
     mu.nValidMuonHits > 0 &&
     fabs(d0) < 0.2 &&
     fabs(dz) < 0.5 &&
     tracks[mu.trackIndex].numberOfValidPixelHits > 0 && 
     (mu.nPixelLayersWithMeasurement + mu.nStripLayersWithMeasurement) > 5 &&
     mu_iso / mu.momentum.Pt() < 0.2)
    return true;
  
  return false;

}
     
     
