#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TObject.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>

#include "SusyEventAnalyzer.h"
#include "BtagWeight.h"
#include "EventQuality.h"

using namespace std;

bool sortTriggers(pair<TString, int> i, pair<TString, int> j) { return (i.second > j.second); }

void SusyEventAnalyzer::PileupWeights(TString puFile) {

  TFile * in = new TFile(puFile, "READ");
  TH1F * _data = (TH1F*)in->Get("pileup");
  
  TString output_code_t = FormatName(scan);

  TH1F * data = (TH1F*)_data->Clone("pu_data"+output_code_t); data->Sumw2();
  TH1F * mc = new TH1F("pu_mc"+output_code_t, "pu_mc"+output_code_t, 70, 0, 70); mc->Sumw2();
  TH1F * mc_nPVertex = new TH1F("mc_nPVertex"+output_code_t, "mc_nPVertex"+output_code_t, 70, 0, 70);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    int nPV = -1;
    susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
    bool foundInTimeBX = false;
    while((iBX != event.pu.end()) && !foundInTimeBX) {
      if(iBX->BX == 0) {
	nPV = iBX->trueNumInteractions;
	foundInTimeBX = true;
      }
      ++iBX;
    }
    
    if(foundInTimeBX) mc->Fill(nPV);

    // Now find the nPV from reconstruction
    int nPV_reco = GetNumberPV(event);
    mc_nPVertex->Fill(nPV_reco);

  } // end event loop

  TH1D * data_nonorm = (TH1D*)data->Clone("pu_data_nonorm"+output_code_t);
  TH1D * mc_nonorm = (TH1D*)mc->Clone("pu_mc_nonorm"+output_code_t);

  Double_t intData = data->Integral();
  Double_t intMC = mc->Integral();

  data->Scale(1./intData);
  mc->Scale(1./intMC);

  TH1F * weights = (TH1F*)data->Clone("puWeights"+output_code_t);
  weights->Divide(mc);

  TFile * out = new TFile("pileupReweighting"+output_code_t+".root", "RECREATE");
  out->cd();

  data->Write();
  data_nonorm->Write();
  mc->Write();
  mc_nonorm->Write();
  weights->Write();
  out->Write();
  out->Close();

  in->Close();

  return;
}

void SusyEventAnalyzer::CalculateBtagEfficiency() {

  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
    nCnt[i][j] = 0;
    }
  }
  
  TString output_code_t = FormatName(scan);

  // open histogram file and define histograms
  TFile * out = new TFile("btagEfficiency"+output_code_t+".root", "RECREATE");
  out->cd();

  TH1F * num_bjets = new TH1F("bjets"+output_code_t, "bjets"+output_code_t, 200, 0, 1000); num_bjets->Sumw2();
  TH1F * num_btags = new TH1F("btags"+output_code_t, "btags"+output_code_t, 200, 0, 1000); num_btags->Sumw2();
  TH1F * num_cjets = new TH1F("cjets"+output_code_t, "cjets"+output_code_t, 200, 0, 1000); num_cjets->Sumw2();
  TH1F * num_ctags = new TH1F("ctags"+output_code_t, "ctags"+output_code_t, 200, 0, 1000); num_ctags->Sumw2();
  TH1F * num_ljets = new TH1F("ljets"+output_code_t, "ljets"+output_code_t, 200, 0, 1000); num_ljets->Sumw2();
  TH1F * num_ltags = new TH1F("ltags"+output_code_t, "ltags"+output_code_t, 200, 0, 1000); num_ltags->Sumw2();

  TH2F * h_DR_jet_gg = new TH2F("DR_jet_gg", "#DeltaR between jets and lead/trailing #gamma#gamma candidates;#DeltaR_{lead #gamma, jet};#DeltaR_{trail #gamma, jet}", 50, 0, 5, 50, 0, 5);

  ScaleFactorInfo sf(btagger);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }

    nCnt[0][0]++; // events

    vector<susy::Photon*> photons;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Muon*> tightMuons, looseMuons;
    vector<susy::Electron*> tightEles, looseEles;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    float HT = 0.;

    findMuons(event, tightMuons, looseMuons, HT);
    findElectrons(event, tightMuons, looseMuons, tightEles, looseEles, HT);

    if(tightMuons.size() + tightEles.size() != 1) continue;

    bool passHLT = true;
    if(useTrigger) {
      if(tightEles.size() == 1) passHLT = PassTriggers(1);
      else if(tightMuons.size() == 1) passHLT = PassTriggers(2) || PassTriggers(3);
    }
    if(!passHLT) continue;

    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     tightMuons, looseMuons,
	     tightEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT, hadronicSystem);

    findPhotons(event, 
		photons,
		pfJets_corrP4,
		tightMuons, looseMuons,
		tightEles, looseEles,
		HT);

    ////////////////////

    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      map<TString, Float_t>::iterator s_it = pfJets[iJet]->jecScaleFactors.find("L1FastL2L3");
      if(s_it == pfJets[iJet]->jecScaleFactors.end()) {
	continue;
      }
      float scale = s_it->second;
      TLorentzVector corrP4 = scale * pfJets[iJet]->momentum;
      if(fabs(corrP4.Eta()) >= 2.4) continue;
	  
      if(fabs(pfJets[iJet]->algDefFlavour) == 5) {
	num_bjets->Fill(corrP4.Pt());
	if((btagger == "CSVL" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.898)) 
	  num_btags->Fill(corrP4.Pt());
      }
	  
      if(fabs(pfJets[iJet]->algDefFlavour) == 4) {
	num_cjets->Fill(corrP4.Pt());
	if((btagger == "CSVL" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.898)) 
	  num_ctags->Fill(corrP4.Pt());
      }
	  
      if(fabs(pfJets[iJet]->algDefFlavour) == 1 ||
	 fabs(pfJets[iJet]->algDefFlavour) == 2 ||
	 fabs(pfJets[iJet]->algDefFlavour) == 3 ||
	 fabs(pfJets[iJet]->algDefFlavour) == 21) {
	num_ljets->Fill(corrP4.Pt());
	if((btagger == "CSVL" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.898)) 
	  num_ltags->Fill(corrP4.Pt());
      }
	  
    } // for jets
	
  } // for entries

  TH1F * bEff = (TH1F*)num_btags->Clone("bEff"+output_code_t);
  bEff->Divide(num_bjets);

  TH1F * cEff = (TH1F*)num_ctags->Clone("cEff"+output_code_t);
  cEff->Divide(num_cjets);

  TH1F * lEff = (TH1F*)num_ltags->Clone("lEff"+output_code_t);
  lEff->Divide(num_ljets);

  out->Write();
  out->Close();

}

void SusyEventAnalyzer::Data() {

  TFile* out = new TFile("hist_"+outputName+"_"+btagger+".root", "RECREATE");
  out->cd();

  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
      nCnt[i][j] = 0;
    }
  }

  ///////////////////////////////////////////////////
  // Define histograms to be filled for all events
  ///////////////////////////////////////////////////

  TString metFilterNames[susy::nMetFilters] = {
    "CSCBeamHalo",
    "HcalNoise",
    "EcalDeadCellTP",
    "EcalDeadCellBE",
    "TrackingFailure",
    "EEBadSC",
    "HcalLaserOccupancy",
    "HcalLaserEventList",
    "HcalLaserRECOUserStep",
    "EcalLaserCorr",
    "ManyStripClus53X",
    "TooManyStripClus53X",
    "LogErrorTooManyClusters",
    "LogErrorTooManyTripletsPairs",
    "LogErrorTooManySeeds",
    "EERingOfFire",
    "InconsistentMuon",
    "GreedyMuon"};

  TH2F* h_metFilter = new TH2F("metFilter", "MET Filter Failures", susy::nMetFilters, 0, susy::nMetFilters, susy::nMetFilters, 0, susy::nMetFilters);
  for(int i = 0; i < susy::nMetFilters; i++) {
    h_metFilter->GetXaxis()->SetBinLabel(i+1, metFilterNames[i]);
    h_metFilter->GetYaxis()->SetBinLabel(i+1, metFilterNames[i]);
  }

  VTH2F h_diempt = BookTH2FVector("diempt", "di-EM pt vs nJets;diEMPt (GeV/c);nJets", 1000, 0., 2000., 200, 0., 200., nCategories, categories, nChannels, channels);
  VTH2F h_dijetpt = BookTH2FVector("dijetpt", "di-Jet pt vs nJets;diJetPt (GeV/c);nJets", 1000, 0., 2000., 200, 0., 200., nCategories, categories, nChannels, channels);
  
  TH2F * h_DR_jet_gg = new TH2F("DR_jet_gg", "#DeltaR between jets and lead/trailing #gamma#gamma candidates;#DeltaR_{lead #gamma, jet};#DeltaR_{trail #gamma, jet}", 50, 0, 5, 50, 0, 5);

  const int nDivisions_chi2 = 50;

  TH1F * h_met_varyCSVcut_ff_j[nDivisions_chi2];
  TH1F * h_met_varyCSVcut_gg_j[nDivisions_chi2];
  TH1F * h_met_squareCSVcut_ff_jj[nDivisions_chi2];
  TH1F * h_met_squareCSVcut_gg_jj[nDivisions_chi2];

  for(int i = 0; i < nDivisions_chi2; i++) {
    char tmp[10];
    double tmp_val = (double)i / (double)nDivisions_chi2;
    sprintf(tmp, "%f", tmp_val);
    TString tmp_t = tmp;

    h_met_varyCSVcut_ff_j[i] = new TH1F("met_varyCSVcut_ff_j_"+tmp_t, "MET for ff+j (maxCSV > "+tmp_t+")", 400, 0, 2000);
    h_met_varyCSVcut_gg_j[i] = new TH1F("met_varyCSVcut_gg_j_"+tmp_t, "MET for gg+j (maxCSV > "+tmp_t+")", 400, 0, 2000);
    h_met_squareCSVcut_ff_jj[i] = new TH1F("met_squareCSVcut_ff_jj_"+tmp_t, "MET for ff+jj (maxCSV & submaxCSV > "+tmp_t+")", 400, 0, 2000);
    h_met_squareCSVcut_gg_jj[i] = new TH1F("met_squareCSVcut_gg_jj_"+tmp_t, "MET for gg+jj (maxCSV & submaxCSV > "+tmp_t+")", 400, 0, 2000);
  }

  /////////////////////////////////
  // Reweighting trees
  /////////////////////////////////

  const int nTreeVariables = 67;

  TString varNames[nTreeVariables] = {
    "pfMET", "pfMET_x", "pfMET_y", "pfMET_phi",
    "pfMET_sysShift", "pfMET_sysShift_phi",
    "pfMET_t1", "pfMET_t1p2", "pfMET_t01", "pfMET_t01p2", "pfNoPUMET", "pfMVAMET",
    "Njets", "Nbtags", "Nphotons", "Nmuons", "Nelectrons",
    "HT", "HT_jets", "hadronic_pt", "w_mT", "m3",
    "ele_pt", "ele_phi", "ele_eta", "ele_mvaTrigV0", "ele_relIso",
    "muon_pt", "muon_phi", "muon_eta", "muon_relIso",
    "leadPhotonEt", "leadPhotonEta", "leadPhotonPhi", "leadChargedHadronIso", "leadSigmaIetaIeta", "lead_nPixelSeeds", "leadMVAregEnergy", "leadMVAregErr",
    "trailPhotonEt", "trailPhotonEta", "trailPhotonPhi", "trailChargedHadronIso", "trailSigmaIetaIeta", "trail_nPixelSeeds", "trailMVAregEnergy", "trailMVAregErr",
    "photon_invmass", "photon_dR", "photon_dPhi", "diEMpT", "diJetPt",
    "jet1_pt", "jet2_pt", "jet3_pt", "jet4_pt",
    "btag1_pt", "btag2_pt",
    "max_csv", "submax_csv", "min_csv",
    "nPV",
    "metFilterBit",
    "runNumber", "eventNumber", "lumiBlock", "jentry"};
    
  map<TString, float> treeMap;
  for(int i = 0; i < nTreeVariables; i++) treeMap[varNames[i]] = 0.;

  vector<TTree*> eventTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_EvtTree", "An event tree for final analysis");
    
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");

    eventTrees.push_back(tree);
  }
  
  ScaleFactorInfo sf(btagger);

  // to check duplicate events
  map<int, set<int> > allEvents;

  bool quitAfterProcessing = false;

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

  if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }
    
    if(useSyncFile) {
      bool sync = false;
      for(unsigned int i = 0; i < syncRuns.size(); i++) {
	//if(event.runNumber == syncRuns[i] && event.luminosityBlockNumber == syncLumi[i] && event.eventNumber == syncEvents[i]) {
	if(event.runNumber == syncRuns[i] && event.eventNumber == syncEvents[i]) {
	  sync = true;
	  //Print(*event);
	  break;
	}
      }
      if(!sync) continue;

      //if(nCnt[0][0] == (syncRuns.size() - 1)) quitAfterProcessing = true;
    }

    if(singleEvent) {
      if(event.runNumber != single_run || event.luminosityBlockNumber != single_lumi || event.eventNumber != single_event) continue;
      //Print(event);
      quitAfterProcessing = true;
    }

    FillMetFilter2D(event, h_metFilter);

    nCnt[0][0]++; // events

    if(useJson && event.isRealData && !IsGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;
    nCnt[1][0]++;

    if(event.isRealData) {
      if(event.passMetFilters() != 1 ||
	 event.passMetFilter(susy::kEcalLaserCorr) != 1 ||
	 event.passMetFilter(susy::kManyStripClus53X) != 1 ||
	 event.passMetFilter(susy::kTooManyStripClus53X) != 1) {
	nCnt[21][0]++;
	continue;
      }
    }

    vector<susy::Muon*> tightMuons, looseMuons;
    vector<susy::Electron*> tightEles, looseEles;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Photon*> photons;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) {
      nCnt[22][0]++;
      continue;
    }

    float HT = 0.;

    findMuons(event, tightMuons, looseMuons, HT);
    findElectrons(event, tightMuons, looseMuons, tightEles, looseEles, HT);
    
    if(tightMuons.size() + tightEles.size() == 0) {
      nCnt[23][0]++;
      continue;
    }

    if(tightMuons.size() + tightEles.size() > 1) {
      nCnt[29][0]++;
      continue;
    }

    if(looseMuons.size() + looseEles.size() != 0) {
      nCnt[24][0]++;
      continue;
    }

    bool passHLT = true;
    if(useTrigger) {
      if(tightEles.size() == 1) passHLT = PassTriggers(1);
      else if(tightMuons.size() == 1) passHLT = PassTriggers(2) || PassTriggers(3);
    }
    if(!passHLT) {
      nCnt[25][0]++;
      continue;
    }

    float HT_jets = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     tightMuons, looseMuons,
	     tightEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT_jets, hadronicSystem);

    findPhotons(event, 
		photons,
		pfJets_corrP4,
		tightMuons, looseMuons,
		tightEles, looseEles,
		HT);
		
    bool duplicateEvent = ! (allEvents[event.runNumber].insert(event.eventNumber)).second;
    if(event.isRealData && duplicateEvent) continue;

    SetTreeValues(treeMap,
		  event,
		  tightMuons, tightEles, 
		  pfJets, btags,
		  photons,
		  pfJets_corrP4, btags_corrP4,
		  csvValues,
		  hadronicSystem,
		  HT, HT_jets,
		  nPVertex,
		  0, 0,
		  jentry);

    ////////////////////

      for(unsigned int chan = 0; chan < nChannels; chan++) {

	if(pfJets.size() < nJetReq[chan]) continue;
	if(btags.size() < nBtagReq[chan]) continue;

	if(tightEles.size() != nEleReq[chan]) continue;
	if(tightMuons.size() != nMuonReq[chan]) continue;

	if(isQCDChannel[chan]) {
	  if(tightMuons.size() == nMuonReq[chan] && isIsolatedMuon(*tightMuons[0])) continue;
	  if((int)(tightEles.size()) == nEleReq[chan] && isIsolatedElectron(*tightEles[0], event.superClusters, event.rho25)) continue;
	}
	else {
	  if((int)(tightMuons.size()) == nMuonReq[chan] && useTrigger && !PassTriggers(2)) continue;
	  if((int)(tightMuons.size()) == nMuonReq[chan] && !isIsolatedMuon(*tightMuons[0])) continue;
	  if((int)(tightEles.size()) == nEleReq[chan] && !isIsolatedElectron(*tightEles[0], event.superClusters, event.rho25)) continue;
	}

	nCnt[2][chan]++;
	eventTrees[chan]->Fill();

      } // loop over jet/btag req channels

    ///////////////////////////////////
    
    if(quitAfterProcessing) break;
  } // for entries
  
  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "in_JSON              : " << nCnt[1][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "----------------" << channels[i] << " Requirement-------------" << endl;
    cout << "gg+" << channels[i] << " events              : " << nCnt[2][i] << endl;
    cout << "eg+" << channels[i] << " events              : " << nCnt[3][i] << endl;
    cout << "ff+" << channels[i] << " events              : " << nCnt[4][i] << endl;
    cout << "gf+" << channels[i] << " events              : " << nCnt[5][i] << endl;
    cout << "ee+" << channels[i] << " events              : " << nCnt[6][i] << endl;
    cout << "ef+" << channels[i] << " events              : " << nCnt[7][i] << endl;
  }
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
  cout << "fail MET filters         : " << nCnt[21][0] << endl;
  cout << "No primary vertex        : " << nCnt[22][0] << endl;
  cout << "zero tight leptons       : " << nCnt[23][0] << endl;
  cout << "2+ tight leptons         : " << nCnt[29][0] << endl;
  cout << "1+ loose leptons         : " << nCnt[24][0] << endl;
  cout << "fail HLT                 : " << nCnt[25][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;

  out->cd();
  out->Write();
  out->Close();

}

void SusyEventAnalyzer::Acceptance() {

  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
    nCnt[i][j] = 0;
    }
  }
  
  TString output_code_t = FormatName(scan);

  // open histogram file and define histograms
  TFile * out = new TFile("signal_contamination"+output_code_t+".root", "RECREATE");
  out->cd();

  TH1D * h_nEvents = new TH1D("nEvents"+output_code_t, "nEvents"+output_code_t, 1, 0, 1);

  const int nTreeVariables = 71;

  TString varNames[nTreeVariables] = {
    "pfMET", "pfMET_x", "pfMET_y", "pfMET_phi",
    "pfMET_sysShift", "pfMET_sysShift_phi",
    "pfMET_t1", "pfMET_t1p2", "pfMET_t01", "pfMET_t01p2", "pfNoPUMET", "pfMVAMET", "genMET",
    "Njets", "Nbtags", "Nphotons", "Nmuons", "Nelectrons",
    "HT", "HT_jets", "hadronic_pt", "w_mT", "m3",
    "ele_pt", "ele_phi", "ele_eta", "ele_mvaTrigV0", "ele_relIso",
    "muon_pt", "muon_phi", "muon_eta", "muon_relIso",
    "leadPhotonEt", "leadPhotonEta", "leadPhotonPhi", "leadChargedHadronIso", "leadSigmaIetaIeta", "lead_nPixelSeeds", "leadMVAregEnergy", "leadMVAregErr",
    "trailPhotonEt", "trailPhotonEta", "trailPhotonPhi", "trailChargedHadronIso", "trailSigmaIetaIeta", "trail_nPixelSeeds", "trailMVAregEnergy", "trailMVAregErr",
    "photon_invmass", "photon_dR", "photon_dPhi", "diEMpT", "diJetPt",
    "jet1_pt", "jet2_pt", "jet3_pt", "jet4_pt",
    "btag1_pt", "btag2_pt",
    "max_csv", "submax_csv", "min_csv",
    "nPV",
    "pileupWeight", "pileupWeightErr",
    "btagWeight", "btagWeightUp", "btagWeightDown", "btagWeightErr",
    "metFilterBit",
    "ttbarDecayMode"};
    
  map<TString, float> treeMap;
  for(int i = 0; i < nTreeVariables; i++) treeMap[varNames[i]] = 0.;

  vector<TTree*> eventTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("gg_"+channels[i]+"_EvtTree"+output_code_t, "An event tree for final analysis");
    
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");

    eventTrees.push_back(tree);
  }

  ScaleFactorInfo sf(btagger);
  TFile * btagEfficiency = new TFile("btagEfficiency"+output_code_t+".root", "READ");
  sf.SetTaggingEfficiencies((TH1F*)btagEfficiency->Get("lEff"+output_code_t), (TH1F*)btagEfficiency->Get("cEff"+output_code_t), (TH1F*)btagEfficiency->Get("bEff"+output_code_t));

  // get pileup weights
  TFile * puFile = new TFile("pileupReweighting"+output_code_t+".root", "READ");
  TH1F * puWeights = (TH1F*)puFile->Get("puWeights"+output_code_t);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;
  h_nEvents->Fill(0., (Double_t)nEntries);

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }

    nCnt[0][0]++; // events

    float numTrueInt = -1.;
    susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
    bool foundInTimeBX = false;
    while((iBX != event.pu.end()) && !foundInTimeBX) {
      if(iBX->BX == 0) {
	numTrueInt = iBX->trueNumInteractions;
	foundInTimeBX = true;
      }
      ++iBX;
    }

    float eventWeight = 0.;
    float eventWeightErr = 0.;
    if(numTrueInt >= 0.) {
      int binNum = puWeights->GetXaxis()->FindBin(numTrueInt);
      eventWeight = puWeights->GetBinContent(binNum);
      eventWeightErr = puWeights->GetBinError(binNum);
    }

    if(!doPileupReweighting) {
      eventWeight = 1.;
      eventWeightErr = 0.;
    }

    vector<susy::Muon*> tightMuons, looseMuons;
    vector<susy::Electron*> tightEles, looseEles;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Photon*> photons;
    
    vector<BtagInfo> tagInfos;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    float HT = 0.;

    findMuons(event, tightMuons, looseMuons, HT);
    findElectrons(event, tightMuons, looseMuons, tightEles, looseEles, HT);
    
    if(tightMuons.size() + tightEles.size() != 1) continue;
    if(looseMuons.size() + looseEles.size() != 0) continue;

    bool passHLT = true;
    if(useTrigger) {
      if(tightEles.size() == 1) passHLT = PassTriggers(1);
      else if(tightMuons.size() == 1) passHLT = PassTriggers(2) || PassTriggers(3);
    }
    if(!passHLT) continue;

    float HT_jets = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     tightMuons, looseMuons,
	     tightEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT_jets, hadronicSystem);

    findPhotons(event, 
		photons,
		pfJets_corrP4,
		tightMuons, looseMuons,
		tightEles, looseEles,
		HT);

    float btagWeight[nChannels];
    float btagWeightUp[nChannels];
    float btagWeightDown[nChannels];
    float btagWeightError[nChannels];
    for(int chan = 0; chan < nChannels; chan++) {
      BtagWeight * tagWeight = new BtagWeight(nBtagReq[chan]);
      pair<float, float> weightResult = tagWeight->weight(tagInfos, btags.size(), 0., false);
      btagWeight[chan] = weightResult.first;
      btagWeightError[chan] = weightResult.second;

      btagWeightUp[chan] = (tagWeight->weight(tagInfos, btags.size(), 1., true)).first;
      btagWeightDown[chan] = (tagWeight->weight(tagInfos, btags.size(), -1., true)).first;

      delete tagWeight;
    }
    tagInfos.clear();

    SetTreeValues(treeMap,
		  event,
		  tightMuons, tightEles, 
		  pfJets, btags,
		  photons,
		  pfJets_corrP4, btags_corrP4,
		  csvValues,
		  hadronicSystem,
		  HT, HT_jets,
		  nPVertex,
		  eventWeight, eventWeightErr,
		  0);

    ////////////////////

    for(int chan = 0; chan < nChannels; chan++) {
      
      if(pfJets.size() < nJetReq[chan]) continue;
      if(btags.size() < nBtagReq[chan]) continue;
      
      if(tightEles.size() != nEleReq[chan]) continue;
      if(tightMuons.size() != nMuonReq[chan]) continue;

      if(isQCDChannel[chan]) {
	if((int)(tightMuons.size()) == nMuonReq[chan] && isIsolatedMuon(*tightMuons[0])) continue;
	if((int)(tightEles.size()) == nEleReq[chan] && isIsolatedElectron(*tightEles[0], event.superClusters, event.rho25)) continue;
      }
      else {
	if((int)(tightMuons.size()) == nMuonReq[chan] && useTrigger && !PassTriggers(2)) continue;
	if((int)(tightMuons.size()) == nMuonReq[chan] && !isIsolatedMuon(*tightMuons[0])) continue;
	if((int)(tightEles.size()) == nEleReq[chan] && !isIsolatedElectron(*tightEles[0], event.superClusters, event.rho25)) continue;
      }

      treeMap["btagWeight"] = btagWeight[chan];
      treeMap["btagWeightErr"] = btagWeightError[chan];
      treeMap["btagWeightUp"] = btagWeightUp[chan];
      treeMap["btagWeightDown"] = btagWeightDown[chan];
      
      nCnt[2][chan]++;
      eventTrees[chan]->Fill();
      
    } // for channels
    
  } // for entries

  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "----------------" << channels[i] << " Requirement-------------" << endl;
    cout << "gg+" << channels[i] << " events              : " << nCnt[2][i] << endl;
    cout << "eg+" << channels[i] << " events              : " << nCnt[3][i] << endl;
    cout << "ff+" << channels[i] << " events              : " << nCnt[4][i] << endl;
    cout << "gf+" << channels[i] << " events              : " << nCnt[5][i] << endl;
  }
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
 
  puFile->Close();
  btagEfficiency->Close();

  out->Write();
  out->Close();

}

void SusyEventAnalyzer::ttggStudy() {
 
}

void SusyEventAnalyzer::SignalContent_gg() {
 
}

void SusyEventAnalyzer::PhotonInfo() {
  
}

void SusyEventAnalyzer::qcdStudy() {

  TFile* out = new TFile("hist_"+outputName+"_"+btagger+".root", "RECREATE");
  out->cd();

  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
      nCnt[i][j] = 0;
    }
  }

  ///////////////////////////////////////////////////
  // Define histograms to be filled for all events
  ///////////////////////////////////////////////////

  TH2D * h_tightEle_mva_iso = new TH2D("tightEle_mva_iso", "MVA vs relIso for tight electrons", 20, -1, 1, 50, 0, 10);

  TH1D * h_tightMuon_iso = new TH1D("tightMuon_iso", "relIso for tight muons", 200, 0, 20);
  TH1D * h_tightMuon_iso_isoHLT = new TH1D("tightMuon_iso_isoHLT", "relIso for tight muons from isoHLT", 200, 0, 20);
  TH1D * h_tightMuon_iso_nonisoHLT = new TH1D("tightMuon_iso_nonisoHLT", "relIso for tight muons from nonisoHLT", 200, 0, 20);

  ScaleFactorInfo sf(btagger);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

  if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }
    
    if(useJson && event.isRealData && !IsGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;

    if(event.isRealData) {
      if(event.passMetFilters() != 1 ||
	 event.passMetFilter(susy::kEcalLaserCorr) != 1 ||
	 event.passMetFilter(susy::kManyStripClus53X) != 1 ||
	 event.passMetFilter(susy::kTooManyStripClus53X) != 1) {
	nCnt[21][0]++;
	continue;
      }
    }

    vector<susy::Muon*> tightMuons, looseMuons;
    vector<susy::Electron*> tightEles, looseEles;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Photon*> photons;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) {
      nCnt[22][0]++;
      continue;
    }

    map<TString, vector<susy::Muon> >::iterator muMap = event.muons.find("muons");
    if(muMap != event.muons.end()) {
      for(vector<susy::Muon>::iterator mu_it = muMap->second.begin(); mu_it != muMap->second.end(); mu_it++) {

	if((int)mu_it->bestTrackIndex() >= (int)(event.tracks).size() || (int)mu_it->bestTrackIndex() < 0) continue;

	bool hasTracks = (int)mu_it->trackIndex < (int)event.tracks.size() && 
	  (int)mu_it->standAloneTrackIndex < (int)event.tracks.size() && 
	  (int)mu_it->combinedTrackIndex < (int)event.tracks.size() && 
	  (int)mu_it->bestTrackIndex() < (int)event.tracks.size() && 
	  (int)mu_it->bestTrackIndex() >= 0;

	if(!hasTracks) continue;
  
	bool isTight = mu_it->isGlobalMuon() && 
	  mu_it->isPFMuon() && 
	  event.tracks[mu_it->combinedTrackIndex].normChi2() < 10. && 
	  mu_it->nValidMuonHits > 0 && 
	  mu_it->nMatchedStations > 1 &&
	  fabs(d0correction(event.vertices[0].position, event.tracks[mu_it->bestTrackIndex()])) < 0.2 &&
	  fabs(dZcorrection(event.vertices[0].position, event.tracks[mu_it->bestTrackIndex()])) < 0.5 &&
	  event.tracks[mu_it->trackIndex].numberOfValidPixelHits > 0 && 
	  (mu_it->nPixelLayersWithMeasurement + mu_it->nStripLayersWithMeasurement) > 5 && 
	  mu_it->momentum.Pt() > 30. &&
	  fabs(mu_it->momentum.Eta()) < 2.1;
	// Skipped relIso < 0.12 cut

	bool isLoose = isVetoMuon(*mu_it);

	if(isTight) tightMuons.push_back(&*mu_it);
	else if(isLoose) looseMuons.push_back(&*mu_it);
	
      }
    }
    
    map<TString, vector<susy::Electron> >::iterator eleMap = event.electrons.find("gsfElectrons");
    if(eleMap != event.electrons.end()) {
      for(vector<susy::Electron>::iterator ele_it = eleMap->second.begin(); ele_it != eleMap->second.end(); ele_it++) {
	
	if((int)ele_it->gsfTrackIndex >= (int)(event.tracks).size() || (int)ele_it->gsfTrackIndex < 0) continue;
	
	bool overlapsMuon = false;
	for(unsigned int i = 0; i < tightMuons.size(); i++) {
	  if(deltaR(tightMuons[i]->momentum, ele_it->momentum) <= 0.5) {
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
	if(overlapsMuon) continue;

	if((int)ele_it->superClusterIndex >= (int)event.superClusters.size() || (int)ele_it->superClusterIndex < 0) continue;
	
	bool isTight = fabs(event.superClusters[ele_it->superClusterIndex].position.Eta()) < 2.5 &&
		       ele_it->momentum.Pt() > 30. &&
	               fabs(d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex])) < 0.02 &&
	               fabs(dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex])) < 1.0 &&
		       ele_it->passConversionVeto &&
		       ele_it->nMissingHits <= 0;
        // Skipped relIso < 0.1 and MVA > 0.5
	
	bool isLoose = isLooseElectron(*ele_it,
				       event.superClusters,
				       event.rho25,
				       d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]),
				       dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
	
	if(isTight) tightEles.push_back(&*ele_it);
	else if(isLoose) looseEles.push_back(&*ele_it);

      }
    }

    if(tightEles.size() + tightMuons.size() == 0) continue;
    if(tightEles.size() + tightMuons.size() > 1) continue;
    if(looseMuons.size() + looseEles.size() != 0) continue;

    bool passHLT = true;
    if(useTrigger) {
      if(tightEles.size() == 1) passHLT = PassTriggers(1);
      else if(tightMuons.size() == 1) passHLT = PassTriggers(2) || PassTriggers(3);
    }
    if(!passHLT) {
      nCnt[25][0]++;
      continue;
    }

    float HT = 0;
    float HT_jets = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     tightMuons, looseMuons,
	     tightEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT_jets, hadronicSystem);

    if(pfJets.size() < 3 || btags.size() < 1) continue;

    for(unsigned int i = 0; i < tightEles.size(); i++) {
      float ele_eta = fabs(event.superClusters[tightEles[i]->superClusterIndex].position.Eta());

      float ea;
      if(ele_eta < 1.0) ea = 0.13;
      else if(ele_eta < 1.479) ea = 0.14;
      else if(ele_eta < 2.0) ea = 0.07;
      else if(ele_eta < 2.2) ea = 0.09;
      else if(ele_eta < 2.3) ea = 0.11;
      else if(ele_eta < 2.4) ea = 0.11;
      else ea = 0.14;

      float ele_iso = max(0., (double)(tightEles[i]->photonIso + tightEles[i]->neutralHadronIso - event.rho25*ea));
      ele_iso += tightEles[i]->chargedHadronIso;

      h_tightEle_mva_iso->Fill(tightEles[i]->mvaTrig, ele_iso / tightEles[i]->momentum.Pt());
    }
        
    for(unsigned int i = 0; i < tightMuons.size(); i++) {
      float mu_iso = max(0., (double)(tightMuons[i]->sumNeutralHadronEt04 + tightMuons[i]->sumPhotonEt04 - 0.5*(tightMuons[i]->sumPUPt04)));
      mu_iso += tightMuons[i]->sumChargedHadronPt04;
      float mu_pt = tightMuons[i]->momentum.Pt();
      h_tightMuon_iso->Fill(mu_iso / mu_pt);
      if(PassTriggers(2) && !PassTriggers(3)) h_tightMuon_iso_isoHLT->Fill(mu_iso / mu_pt);
      if(!PassTriggers(2) && PassTriggers(3)) h_tightMuon_iso_nonisoHLT->Fill(mu_iso / mu_pt);
    }
    
    ///////////////////////////////////
    
  } // for entries
  
  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "in_JSON              : " << nCnt[1][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "----------------" << channels[i] << " Requirement-------------" << endl;
    cout << "gg+" << channels[i] << " events              : " << nCnt[2][i] << endl;
    cout << "eg+" << channels[i] << " events              : " << nCnt[3][i] << endl;
    cout << "ff+" << channels[i] << " events              : " << nCnt[4][i] << endl;
    cout << "gf+" << channels[i] << " events              : " << nCnt[5][i] << endl;
    cout << "ee+" << channels[i] << " events              : " << nCnt[6][i] << endl;
    cout << "ef+" << channels[i] << " events              : " << nCnt[7][i] << endl;
  }
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
  cout << "fail MET filters         : " << nCnt[21][0] << endl;
  cout << "No primary vertex        : " << nCnt[22][0] << endl;
  cout << "zero tight leptons       : " << nCnt[23][0] << endl;
  cout << "2+ tight leptons         : " << nCnt[29][0] << endl;
  cout << "1+ loose leptons         : " << nCnt[24][0] << endl;
  cout << "fail HLT                 : " << nCnt[25][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;

  out->cd();
  out->Write();
  out->Close();  
}
