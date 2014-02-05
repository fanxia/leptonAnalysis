#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TArrayI.h>

#include <map>
#include <set>
#include <fstream>

#include "../src/SusyEvent.h"
#include "Muon.h"
#include "Electron.h"
#include "Photon.h"
#include "Jet.h"
#include "ScaleFactorInfo.h"
#include "BtagInfo.h"
#include "Metadata.h"

using namespace std;

template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}

bool CorrPtGreater(TLorentzVector p1, TLorentzVector p2) {
  return (p1.Pt() > p2.Pt());
}

class HLTInfo {
 public :
  HLTInfo() {
    paths.clear();
    types.clear();
  }
  HLTInfo(vector<TString> pathList, vector<int> eventType) {
    SetHLTPaths(pathList);
    SetEventTypes(eventType);
  }
  virtual ~HLTInfo() { ; }

  vector<TString> GetHLTNames() { return paths; }
  vector<int> GetEventTypes() { return types; }

  void SetHLTPaths(vector<TString> p) { 
    paths.clear();
    for(unsigned int i = 0; i < p.size(); i++) paths.push_back(p[i]);
  }
  void SetEventTypes(vector<int> t) { 
    types.clear();
    for(unsigned int i = 0; i < t.size(); i++) types.push_back(t[i]);
  }
  void AddEventTypes(int s) { types.push_back(s); }

 private :
  vector<TString> paths;
  vector<int> types;
};

class SusyEventAnalyzer {
 public:
  SusyEventAnalyzer(TTree&);
  virtual ~SusyEventAnalyzer();

  virtual void Data();
  virtual void Acceptance();
  virtual void ttggStudy();
  virtual void CalculateBtagEfficiency();
  virtual void PileupWeights(TString puFile);
  virtual void SignalContent_gg();
  virtual void PhotonInfo();

  // utility functions
  float deltaR(TLorentzVector& p1, TLorentzVector& p2);
  float deltaR(TLorentzVector& p1, TVector3& p2);
  float deltaR(TVector3& p1, TVector3& p2);
  TVector3 FindJetVertex(susy::PFJet jet, vector<susy::Track>& tracks);
  double d0correction(TVector3& beamSpot, susy::Track& track) const;
  double dZcorrection(TVector3& beamSpot, susy::Track& track) const;
  void AddHlt(vector<TString> v, vector<int> eventType) {
    hltInfos.push_back(HLTInfo(v, eventType));
  }

  // metadata and btag stuff
  void SetUseSyncFile(bool v) { useSyncFile = v; }
  void IncludeSyncFile(char* file);
  void SetCheckSingleEvent(bool v) { singleEvent = v; }
  void AddCheckSingleEvent(Int_t run, Int_t lumi, ULong_t eventnum) {
    single_run = run;
    single_lumi = lumi;
    single_event = eventnum;
  }
  void SetScanName(TString v) { scan = v; }
  void SetIsMC(bool v)        { isMC = v; }
  void SetIsFastSim(bool v) {
    if(v) {
      isMC = v;
      isFastSim = v;
    }
    else isFastSim = v;
  }
  void SetRejectFakeElectrons(bool v)  { rejectFakeElectrons = v;}
  void SetDoBtagScaling(bool v)        { doBtagScaling = v; }
  void SetBtagTechnicalStop(TString t) { btagTechnicalStop = t; }
  void SetDoPileupReweighting(bool v)  { doPileupReweighting = v; }
  void SetBtagger(TString t) { btagger = t; }
  void AddValidTagger(TString t) { validTags.push_back(t); }
  pair<unsigned int, Float_t> tagEncoder(TString tagger) {
    if(tagger == "CSVL") return make_pair(5, 0.244);
    if(tagger == "CSVM") return make_pair(5, 0.679);
    if(tagger == "CSVT") return make_pair(5, 0.898);
    return make_pair(0, 1.e6);
  }
  void SetUseTrigger(bool v) { useTrigger = v; }
  void SetUseJson(bool v)    { useJson = v; }

  void SetUseDPhiCut(bool v) { useDPhiCut = v; }

  void IncludeAJson(TString const&);
  void SetOutput(TString const& v) { outputName = v; }
  void SetPrintInterval(int v) { printInterval = v; }
  void SetPrintLevel(int v) { printLevel = v; }
  void SetProcessNEvents(int v) { processNEvents = v; }
  void CopyEvents(bool v) { copyEvents = v; }

  // major analysis logic
  void findPhotons(susy::Event& ev, vector<susy::Photon*>& candidates,
		   vector<susy::Muon*> isoMuons, vector<susy::Muon*> looseMuons,
		   vector<susy::Electron*> isoEles, vector<susy::Electron*> looseEles,
		   vector<susy::PFJet*>& pfJets, float& HT)
  void findPhotons_prioritizeEt(susy::Event& ev, vector<susy::Photon*>& candidates, int& event_type, bool doDPhiCut);
  void findPhotons_simple(susy::Event& ev, vector<susy::Photon*>& candidates, int& event_type, int wp, bool doDPhiCut);
  void findPhotons_fakesWithSeeds(susy::Event& ev, vector<susy::Photon*>& candidates, int& event_type, int wp, bool doDPhiCut);
  void findJets(susy::Event& ev, vector<susy::Photon*> candidates,
		vector<susy::Muon*> isoMuons, vector<susy::Muon*> looseMuons,
		vector<susy::Electron*> isoEles, vector<susy::Electron*> looseEles,
		vector<susy::PFJet*>& pfJets, vector<susy::PFJet*>& btags, 
		ScaleFactorInfo sf,
		vector<BtagInfo>& tagInfos, vector<float>& csvValues, 
		vector<TLorentzVector>& pfJets_corrP4, vector<TLorentzVector>& btags_corrP4, 
		float& HT, TLorentzVector& hadronicSystem,
		TH2F*& h_DR_jet_gg);
  void findMuons(susy::Event& ev, vector<susy::Photon*> candidates, vector<susy::Muon*>& isoMuons, vector<susy::Muon*>& looseMuons, float& HT);
  void findElectrons(susy::Event& ev, vector<susy::Photon*> candidates, vector<susy::Electron*>& isoEles, vector<susy::Electron*>& looseElese, float& HT);
  bool GetDiJetPt(susy::Event& ev, vector<susy::Photon*> candidates, float& diJetPt, float& leadpt, float& trailpt);
  bool PhotonMatchesElectron(susy::Event& ev, vector<susy::Photon*> candidates, int& bothMatchCounter);
  int FigureTTbarDecayMode(susy::Event& ev);

  // lazy junk
  void FillMetFilter2D(susy::Event& ev, TH2F*& h);

 protected:
  bool IsGoodLumi(UInt_t, UInt_t) const;
  bool PassTrigger(TString path) const;
  bool PassTriggers(int eventType);

  susy::Event event;
  TTree *fTree;
  TString outputName;
  int printLevel;
  unsigned printInterval;
  int processNEvents;
  bool copyEvents;
  map<unsigned, set<unsigned> > goodLumiList;
  mutable pair<unsigned, unsigned> currentLumi;
  mutable bool currentLumiIsGood;

 private:
  vector<HLTInfo> hltInfos;

  vector<UInt_t> syncRuns;
  vector<UInt_t> syncLumi;
  vector<ULong_t> syncEvents;
  UInt_t single_run;
  UInt_t single_lumi;
  ULong_t single_event;
  bool useSyncFile;
  bool singleEvent;
  
  vector<TString> validTags;
  TString btagger;

  bool isMC;
  bool isFastSim;
  bool rejectFakeElectrons;
  bool doBtagScaling;
  TString btagTechnicalStop;
  bool doPileupReweighting;
  bool useTrigger;
  bool useJson;

  bool useDPhiCut;

  TString scan;

};

SusyEventAnalyzer::SusyEventAnalyzer(TTree& tree) :
  event(),
  fTree(&tree),
  outputName("analysis"),
  printLevel(0),
  printInterval(1000),
  processNEvents(-1),
  copyEvents(false),
  goodLumiList(),
  currentLumi(0, 0),
  currentLumiIsGood(true),
  rejectFakeElectrons(false)  
{
  event.setInput(tree);
}

SusyEventAnalyzer::~SusyEventAnalyzer() {}

void SusyEventAnalyzer::IncludeAJson(TString const& _fileName) {
  if(_fileName == "") return;

  ifstream inputFile(_fileName);
  if(!inputFile.is_open()){
    cerr << "Cannot open JSON file " << _fileName << endl;
    return;
  }

  string line;
  TString jsonText;
  while(true){
    getline(inputFile, line);
    if(!inputFile.good()) break;
    jsonText += line;
  }
  inputFile.close();

  TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
  TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

  TArrayI positions(2);
  positions[1] = 0;
  while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
    TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
    TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

    unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
    set<unsigned>& lumis(goodLumiList[run]);

    TArrayI lumiPos(2);
    lumiPos[1] = 0;
    while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
      TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
      int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
      int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
      for(int lumi(begin); lumi <= end; ++lumi)
        lumis.insert(lumi);
    }
  }
}

bool SusyEventAnalyzer::IsGoodLumi(UInt_t run, UInt_t lumi) const {
  if(goodLumiList.size() == 0) return true;
  if(run == currentLumi.first && lumi == currentLumi.second) return currentLumiIsGood;
  currentLumi.first = run;
  currentLumi.second = lumi;
  currentLumiIsGood = false;

  map<unsigned, set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
  if(rItr != goodLumiList.end()){
    set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
    if(lItr != rItr->second.end()) currentLumiIsGood = true;
  }

  return currentLumiIsGood;
}

bool SusyEventAnalyzer::PassTrigger(TString path) const {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event.hltMap.begin(); it != event.hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second))) {
      pass = true;
      break;
    }
  }
  return pass;
}

bool SusyEventAnalyzer::PassTriggers(int eventType) {

  // Types: gg = 1
  //	    eg = 2
  // 	    ee = 3
  // 	    ff = 4
  //        gf = 5

  bool pass = false;

  for(unsigned int i = 0; i < hltInfos.size(); i++) {

    bool correctType = false;
    for(unsigned int j = 0; j < ((hltInfos[i]).GetEventTypes()).size(); j++) {
      if(eventType == ((hltInfos[i]).GetEventTypes())[j]) {
	correctType = true;
	break;
      }
    }
    if(!correctType) continue;

    for(unsigned int j = 0; j < ((hltInfos[i]).GetHLTNames()).size(); j++) {
      if(PassTrigger( ((hltInfos[i]).GetHLTNames())[j]) ) {
	pass = true;
	break;
      }
    }

    if(pass == true) break;
  }

  return pass;
}

float SusyEventAnalyzer::deltaR(TLorentzVector& p1, TLorentzVector& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

float SusyEventAnalyzer::deltaR(TLorentzVector& p1, TVector3& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

float SusyEventAnalyzer::deltaR(TVector3& p1, TVector3& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

TVector3 SusyEventAnalyzer::FindJetVertex(susy::PFJet jet, vector<susy::Track>& tracks) {

  vector<TVector3> vertices;
  vector<float> sumPt;
  vector<int> nTracks;

  float epsilon = 1.e-6;

  for(unsigned int i = 0; i < jet.tracklist.size(); i++) {
    TVector3 vtx = tracks[jet.tracklist[i]].vertex;
    bool newVertex = true;
    unsigned int oldVertexIndex = -1;
    for(unsigned int j = 0; j < vertices.size(); j++) {
      if(abs(vtx.X() - vertices[j].X()) < epsilon &&
	 abs(vtx.Y() - vertices[j].Y()) < epsilon &&
	 abs(vtx.Z() - vertices[j].Z()) < epsilon) {
	newVertex = false;
	oldVertexIndex = j;
	break;
      }
    }

    if(newVertex) {
      vertices.push_back(vtx);
      sumPt.push_back(tracks[jet.tracklist[i]].momentum.Pt());
      nTracks.push_back(1);
    }
    else {
      sumPt[oldVertexIndex] += tracks[jet.tracklist[i]].momentum.Pt();
      nTracks[oldVertexIndex]++;
    }

  }

  int maxNtracks = -1;
  float maxSumPt = -1.;
  unsigned int matchByN = 0;
  unsigned int matchBySumPt = 0;

  for(unsigned int i = 0; i < vertices.size(); i++) {
    if(nTracks[i] > maxNtracks) {
      matchByN = i;
      maxNtracks = nTracks[i];
    }
    if(sumPt[i] > maxSumPt) {
      matchBySumPt = i;
      maxSumPt = sumPt[i];
    }
  }

  if(matchByN != matchBySumPt) cout << "Disagreeing nTracks and sumPt vertex-finding!" << endl;

  if(vertices.size() > 0) return vertices[matchByN];

  return TVector3(0.);

}

double SusyEventAnalyzer::d0correction(TVector3& beamspot, susy::Track& track) const {
  double d0 = track.d0() - beamspot.X()*sin(track.phi()) + beamspot.Y()*cos(track.phi());
  return d0;
}

double SusyEventAnalyzer::dZcorrection(TVector3& beamSpot, susy::Track& track) const {

  if(&(beamSpot) == 0x0 || &(track) == 0x0) {
    cout << endl << endl << "Something majorly fucked up in dZcorrection!!!" <<	endl <<	endl;
    return 1.e6;
  }

  if(track.momentum.Pt() == 0.) return 1.e6;

  double dz = (track.vertex.Z() - beamSpot.Z()) - ((track.vertex.X() - beamSpot.X())*track.momentum.Px() + (track.vertex.Y() - beamSpot.Y())*track.momentum.Py()) / track.momentum.Pt() * (track.momentum.Pz() / track.momentum.Pt());
  return dz;
}

void SusyEventAnalyzer::findPhotons(susy::Event& ev, vector<susy::Photon*>& candidates,
				    vector<susy::Muon*> isoMuons, vector<susy::Muon*> looseMuons,
				    vector<susy::Electron*> isoEles, vector<susy::Electron*> looseEles,
				    vector<susy::PFJet*>& pfJets, float& HT) {
  
  map<TString, vector<susy::Photon> >::iterator phoMap = ev.photons.find("photons");
  if(phoMap != event.photons.end()) {
    for(vector<susy::Photon>::iterator it = phoMap->second.begin();
	it != phoMap->second.end(); it++) {
      
      if(is_egf(*it, event.rho25)) {

	bool overlap = false;
	
	for(unsigned int i = 0; i < isoMuons.size(); i++) {
	  if(deltaR(isoMuons[i]->momentum, it->caloPosition) <= 0.5) {
	    overlap = true;
	    break;
	  }
	}

	for(unsigned int i = 0; i < looseMuons.size(); i++) {
	  if(deltaR(looseMuons[i]->momentum, it->caloPosition) <= 0.5) {
	    overlap = true;
	    break;
	  }
	}

	for(unsigned int i = 0; i < isoEles.size(); i++) {
	  if(deltaR(isoEles[i]->momentum, it->caloPosition) <= 0.5) {
	    overlap = true;
	    break;
	  }
	}

	for(unsigned int i = 0; i < looseEles.size(); i++) {
	  if(deltaR(looseEles[i]->momentum, it->caloPosition) <= 0.5) {
	    overlap = true;
	    break;
	  }
	}

	for(unsigned int i = 0; i < pfJets.size(); i++) {
	  if(deltaR(pfJets[i]->momentum, it->caloPosition) <= 0.5) {
	    overlap = true;
	    break;
	  }
	}

	if(overlap) continue;

	photons.push_back(&*it);
	HT += it->momentum.Pt();
	
      }
      
    } // for photon
  } // if

  sortphotons.begin(), photons.end(), EtGreater<susy::Photon>);
    
  return;

}

void SusyEventAnalyzer::findJets(susy::Event& ev, 
				 vector<susy::Muon*> isoMuons, vector<susy::Muon*> looseMuons,
				 vector<susy::Electron*> isoEles, vector<susy::Electron*> looseEles,
				 vector<susy::PFJet*>& pfJets, vector<susy::PFJet*>& btags, 
				 ScaleFactorInfo sf,
				 vector<BtagInfo>& tagInfos, vector<float>& csvValues, 
				 vector<TLorentzVector>& pfJets_corrP4, vector<TLorentzVector>& btags_corrP4, 
				 float& HT, TLorentzVector& hadronicSystem) {

  map<TString, susy::PFJetCollection>::iterator pfJets_it = ev.pfJets.find("ak5");
  if(pfJets_it != ev.pfJets.end()) {
    susy::PFJetCollection& jetColl = pfJets_it->second;
      
    for(vector<susy::PFJet>::iterator it = jetColl.begin();
	it != jetColl.end(); it++) {
	
      map<TString, Float_t>::iterator s_it = it->jecScaleFactors.find("L1FastL2L3");
      float scale = s_it->second;
	  
      TLorentzVector corrP4 = scale * it->momentum;

      if(!isGoodJet(*it, corrP4)) continue;
      if(JetOverlapsElectron(corrP4, isoEles, looseEles)) continue;
      if(JetOverlapsMuon(corrP4, isoMuons, looseMuons)) continue;

      pfJets.push_back(&*it);
      csvValues.push_back(it->bTagDiscriminators[susy::kCSV]);
      HT += corrP4.Pt();
      hadronicSystem += corrP4;
      pfJets_corrP4.push_back(corrP4);
	
      if(fabs(corrP4.Eta()) < 2.4) {
	if(isMC) {
	  BtagInfo info((*it), corrP4, btagger, 1., isMC, isFastSim, sf, btagTechnicalStop, (15906. + 77789.) / 328124., 117192. / 328124., 117237. / 328124.);
	  tagInfos.push_back(info);
	}

	if((btagger == "CSVL" && it->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && it->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && it->bTagDiscriminators[susy::kCSV] > 0.898)) {
	  btags.push_back(&*it);
	  btags_corrP4.push_back(corrP4);
	}
	 
      }
   
    } // loop over jet coll
  } // if the jet coll exists
  sort(pfJets.begin(), pfJets.end(), EtGreater<susy::PFJet>);
  sort(btags.begin(), btags.end(), EtGreater<susy::PFJet>);
  sort(csvValues.begin(), csvValues.end(), greater<float>());
  sort(pfJets_corrP4.begin(), pfJets_corrP4.end(), CorrPtGreater);
  sort(btags_corrP4.begin(), btags_corrP4.end(), CorrPtGreater);

}

void SusyEventAnalyzer::findMuons(susy::Event& ev, vector<susy::Muon*>& isoMuons, vector<susy::Muon*>& looseMuons) {

  map<TString, vector<susy::Muon> >::iterator muMap = ev.muons.find("muons");
  if(muMap != ev.muons.end()) {
    for(vector<susy::Muon>::iterator mu_it = muMap->second.begin(); mu_it != muMap->second.end(); mu_it++) {

      if((int)mu_it->bestTrackIndex() >= (int)(event.tracks).size() || (int)mu_it->bestTrackIndex() < 0) continue;

      if(isTightMuon(*mu_it, 
		     event.tracks, 
		     d0correction(event.vertices[0].position, event.tracks[mu_it->bestTrackIndex()]), 
		     dZcorrection(event.vertices[0].position, event.tracks[mu_it->bestTrackIndex()]))
	 ) {
	isoMuons.push_back(&*mu_it);
      }

      else if(isVetoMuon(*mu_it)) looseMuons.push_back(&*mu_it);

    }

  }

}

void SusyEventAnalyzer::findElectrons(susy::Event& ev, vector<susy::Electron*>& isoMuons, vector<susy::Electron*>& looseMuons, vector<susy::Electron*>& isoEles, vector<susy::Electron*>& looseEles) {

  map<TString, vector<susy::Electron> >::iterator eleMap = ev.electrons.find("gsfElectrons");
  if(eleMap != ev.electrons.end()) {
    for(vector<susy::Electron>::iterator ele_it = eleMap->second.begin(); ele_it != eleMap->second.end(); ele_it++) {

      if((int)ele_it->gsfTrackIndex >= (int)(event.tracks).size() || (int)ele_it->gsfTrackIndex < 0) continue;

      if(isTightElectron(*ele_it, 
			 event.superClusters, 
			 event.rho25, 
			 d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
			 dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]))) {
	
	bool overlapsMuon = false;
	
	for(unsigned int i = 0; i < isoMuons.size(); i++) {
	  if(deltaR(isoMuons[i]->momentum, ele_it->momentum) <= 0.5) {
	    overlapsMuon = true;
	    break;
	  }
	}

	for(unsigned int i = 0; i < looseMuons.size(); i++) {
	  if(deltaR(looseMuons[i]->momentum, ele_it->momentum) <= 0.5) {
	    overlapsMuon = true;
	    break;
	  }
	}

	if(!overlapsMuon) isoEles.push_back(&*ele_it);

      }

      else if(isLooseElectron(*ele_it,
			      event.superClusters, 
			      event.rho25, 
			      d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
			      dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]))) {	 

	bool overlapsMuon = false;
	
	for(unsigned int i = 0; i < isoMuons.size(); i++) {
	  if(deltaR(isoMuons[i]->momentum, ele_it->momentum) <= 0.5) {
	    overlapsMuon = true;
	    break;
	  }
	}

	for(unsigned int i = 0; i < looseMuons.size(); i++) {
	  if(deltaR(looseMuons[i]->momentum, ele_it->momentum) <= 0.5) {
	    overlapsMuon = true;
	    break;
	  }
	}

	if(!overlapsMuon) looseEles.push_back(&*ele_it);
	
      }
     
    }
  }

}

void SusyEventAnalyzer::FillMetFilter2D(susy::Event& ev, TH2F*& h) {

  for(int i = 0; i < susy::nMetFilters; i++) {
    for(int j = 0; j < susy::nMetFilters; j++) {

      if(!ev.passMetFilter(i) && !ev.passMetFilter(j)) h->Fill(i, j);

    }
  }

}

bool SusyEventAnalyzer::GetDiJetPt(susy::Event& ev, vector<susy::Photon*> candidates, float& diJetPt, float& leadpt, float& trailpt) {

  vector<TLorentzVector> jetP4s;

  map<TString, susy::PFJetCollection>::iterator iJets = ev.pfJets.find("ak5");
  if(iJets == ev.pfJets.end()) return false;
  susy::PFJetCollection& jetColl = iJets->second;

  bool worked = true;

  for(vector<susy::Photon*>::iterator it = candidates.begin(); it != candidates.end(); it++) {
      
    bool matched = false;
      
    for(vector<susy::PFJet>::iterator iJet = jetColl.begin(); iJet != jetColl.end(); iJet++) {
      float theJES = 1.0;
      map<TString, Float_t>::const_iterator iCorr1 = iJet->jecScaleFactors.find("L1FastL2L3");
      map<TString, Float_t>::const_iterator iCorr2 = iJet->jecScaleFactors.find("L2L3");
      TLorentzVector corrP4 = iJet->momentum;
      if(iCorr1 != iJet->jecScaleFactors.end()) {
	if(iCorr2 != iJet->jecScaleFactors.end()) {
	  theJES = iCorr1->second / iCorr2->second;
	  corrP4 = theJES * iJet->momentum;
	}
	else {
	  theJES = iCorr1->second;
	  corrP4 = theJES * iJet->momentum;
	}
      }

      float dEta = corrP4.Eta() - (*it)->caloPosition.Eta();
      float dPhi = TVector2::Phi_mpi_pi(corrP4.Phi() - (*it)->caloPosition.Phi());
      float dR = sqrt(dEta*dEta + dPhi*dPhi);

      if(corrP4.Et() > 20. &&
	 fabs(corrP4.Eta()) < 2.6 &&
	 dR < 0.3) {
	  
	jetP4s.push_back(corrP4);
	matched = true;
	break;
      }

    }
      
    if(!matched) {
      worked = false;
      jetP4s.clear();
      jetP4s.push_back(candidates[0]->momentum);
      jetP4s.push_back(candidates[1]->momentum);
      break;
    }
      
  }
  
  sort(jetP4s.begin(), jetP4s.end(), CorrPtGreater);

  diJetPt = (jetP4s[0] + jetP4s[1]).Pt();
  leadpt = jetP4s[0].Pt();
  trailpt = jetP4s[1].Pt();

  return worked;
}

bool SusyEventAnalyzer::PhotonMatchesElectron(susy::Event& ev, vector<susy::Photon*> candidates, int& bothMatchCounter) {

  if(ev.isRealData) return false;

  bool matchesLead = false;
  bool matchesTrail = false;

  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

    if(it->status == 3 && fabs(it->pdgId) == 11) {
      if(deltaR(it->momentum, candidates[0]->caloPosition) < 0.1) matchesLead = true;
      if(deltaR(it->momentum, candidates[1]->caloPosition) < 0.1) matchesTrail = true;
    }

  }

  if(matchesLead && matchesTrail) bothMatchCounter++;

  return (matchesLead || matchesTrail);
}

int SusyEventAnalyzer::FigureTTbarDecayMode(susy::Event& ev) {

  int decayMode = -1;

  int nElectronicWs = 0;
  int nMuonicWs = 0;
  int nTauonicWs = 0;
  
  int firstWindex = -1;
  
  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {
    
    //if(it->status != 3) continue;
    //if(it->momentum.Pt() < 20.) continue;
    
    bool isFromWfromStopOrTop = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
      (
       fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 6 || 
       fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000006
       );
    
    if(!isFromWfromStopOrTop) continue;
    
    if(firstWindex == it->motherIndex) continue;
    
    if(fabs(it->pdgId) == 11) nElectronicWs++;
    else if(fabs(it->pdgId) == 13) nMuonicWs++;
    else if(fabs(it->pdgId) == 15) nTauonicWs++;
    
    if(firstWindex < 0) firstWindex = it->motherIndex;
    
  }
  
  /* decayMode:
     0 hadronic
     1 semi-ele
     2 semi-mu
     3 semi-tau
     4 di-ele
     5 di-mu
     6 di-tau
     7 ele-mu
     8 ele-tau
     9 mu-tau
  */
  
  if((nElectronicWs + nMuonicWs + nTauonicWs) == 0) decayMode = 0;
  else if(nElectronicWs == 1 && (nMuonicWs + nTauonicWs) == 0) decayMode = 1;
  else if(nMuonicWs == 1 && (nElectronicWs + nTauonicWs) == 0) decayMode = 2;
  else if(nTauonicWs == 1 && (nElectronicWs + nMuonicWs) == 0) decayMode = 3;
  else if(nElectronicWs == 2 && (nMuonicWs + nTauonicWs) == 0) decayMode = 4;
  else if(nMuonicWs == 2 && (nElectronicWs + nTauonicWs) == 0) decayMode = 5;
  else if(nTauonicWs == 2 && (nElectronicWs + nMuonicWs) == 0) decayMode = 6;
  else if(nElectronicWs == 1 && nMuonicWs == 1 && nTauonicWs == 0) decayMode = 7;
  else if(nElectronicWs == 1 && nMuonicWs == 0 && nTauonicWs == 1) decayMode = 8;
  else if(nElectronicWs == 0 && nMuonicWs == 1 && nTauonicWs == 1) decayMode = 9;

  return decayMode;
}

void SusyEventAnalyzer::IncludeSyncFile(char* file) {
  FILE* syncList = fopen(file, "r");
  if(!syncList) return;

  char line[256];
  Int_t _run, _lumi;
  ULong_t _event;

  while(fgets(line, 255, syncList)) {
    sscanf(line, "%d %d %lu", &_run, &_lumi, &_event);
    syncRuns.push_back(_run);
    syncLumi.push_back(_lumi);
    syncEvents.push_back(_event);
  }
    
  fclose(syncList);
}

#endif