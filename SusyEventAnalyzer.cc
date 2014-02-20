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
    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    float HT = 0.;

    findMuons(event, isoMuons, looseMuons, HT);
    findElectrons(event, isoMuons, looseMuons, isoEles, looseEles, HT);

    if(isoMuons.size() + isoEles.size() != 1) continue;

    bool passHLT = true;
    if(useTrigger) {
      if(isoEles.size() == 1) passHLT = PassTriggers(1);
      else if(isoMuons.size() == 1) passHLT = PassTriggers(2);
    }
    if(!passHLT) continue;

    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT, hadronicSystem);

    findPhotons(event, 
		photons,
		pfJets_corrP4,
		isoMuons, looseMuons,
		isoEles, looseEles,
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

  const int nTreeVariables = 63;

  TString varNames[nTreeVariables] = {
    "pfMET", "pfMET_x", "pfMET_y", "pfMET_phi",
    "pfMET_sysShift", "pfMET_sysShift_phi",
    "pfMET_t1", "pfMET_t1p2", "pfMET_t01", "pfMET_t01p2", "pfNoPUMET", "pfMVAMET",
    "Njets", "Nbtags", "Nphotons", "Nmuons", "Nelectrons",
    "HT", "HT_jets", "hadronic_pt", "w_mT",
    "ele_pt", "ele_phi", "ele_eta",
    "muon_pt", "muon_phi", "muon_eta",
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
    TTree * tree = new TTree("gg_"+channels[i]+"_EvtTree", "An event tree for final analysis");
    
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

    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
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

    findMuons(event, isoMuons, looseMuons, HT);
    findElectrons(event, isoMuons, looseMuons, isoEles, looseEles, HT);
    
    if(isoMuons.size() + isoEles.size() == 0) {
      nCnt[23][0]++;
      continue;
    }

    if(isoMuons.size() + isoEles.size() > 1) {
      nCnt[29][0]++;
      continue;
    }

    if(looseMuons.size() + looseEles.size() != 0) {
      nCnt[24][0]++;
      continue;
    }

    bool passHLT = true;
    if(useTrigger) {
      if(isoEles.size() == 1) passHLT = PassTriggers(1);
      else if(isoMuons.size() == 1) passHLT = PassTriggers(2);
    }
    if(!passHLT) {
      nCnt[25][0]++;
      continue;
    }

    float HT_jets = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT_jets, hadronicSystem);

    findPhotons(event, 
		photons,
		pfJets_corrP4,
		isoMuons, looseMuons,
		isoEles, looseEles,
		HT);
		
    bool duplicateEvent = ! (allEvents[event.runNumber].insert(event.eventNumber)).second;
    if(event.isRealData && duplicateEvent) continue;

    SetTreeValues(treeMap,
		  event,
		  isoMuons, isoEles, 
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
	if(nEleReq[chan] >= 0 && isoEles.size() != nEleReq[chan]) continue;
	if(nMuonReq[chan] >= 0 && isoMuons.size() != nMuonReq[chan]) continue;

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

  const int nTreeVariables = 68;

  TString varNames[nTreeVariables] = {
    "pfMET", "pfMET_x", "pfMET_y", "pfMET_phi",
    "pfMET_sysShift", "pfMET_sysShift_phi",
    "pfMET_t1", "pfMET_t1p2", "pfMET_t01", "pfMET_t01p2", "pfNoPUMET", "pfMVAMET", "genMET",
    "Njets", "Nbtags", "Nphotons", "Nmuons", "Nelectrons",
    "HT", "HT_jets", "hadronic_pt", "w_mT", "w_mT_genNeutrino",
    "ele_pt", "ele_phi", "ele_eta",
    "muon_pt", "muon_phi", "muon_eta",
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

    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Photon*> photons;
    
    vector<BtagInfo> tagInfos;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    float HT = 0.;

    findMuons(event, isoMuons, looseMuons, HT);
    findElectrons(event, isoMuons, looseMuons, isoEles, looseEles, HT);
    
    if(isoMuons.size() + isoEles.size() != 1) continue;
    if(looseMuons.size() + looseEles.size() != 0) continue;

    bool passHLT = true;
    if(useTrigger) {
      if(isoEles.size() == 1) passHLT = PassTriggers(1);
      else if(isoMuons.size() == 1) passHLT = PassTriggers(2);
    }
    if(!passHLT) continue;

    float HT_jets = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findJets(event, 
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT_jets, hadronicSystem);

    findPhotons(event, 
		photons,
		pfJets_corrP4,
		isoMuons, looseMuons,
		isoEles, looseEles,
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
		  isoMuons, isoEles, 
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
      if(nEleReq[chan] >= 0 && isoEles.size() != nEleReq[chan]) continue;
      if(nMuonReq[chan] >= 0 && isoMuons.size() != nMuonReq[chan]) continue;
      
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

// things to study after a preselection of gg+bj+X:
// jet eta and pt requirements
// njets requirements
// lepton pt and eta reqs
// lepton veto in hadronic mode

void SusyEventAnalyzer::ttggStudy() {
  /*
  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
    nCnt[i][j] = 0;
    }
  }
  
  TString output_code_t = FormatName(scan);

  // open histogram file and define histograms
  TFile * out = new TFile("ttggStudy"+output_code_t+".root", "RECREATE");
  out->cd();

  TH2F * h_DR_jet_gg = new TH2F("DR_jet_gg", "#DeltaR between jets and lead/trailing #gamma#gamma candidates;#DeltaR_{lead #gamma, jet};#DeltaR_{trail #gamma, jet}", 50, 0, 5, 50, 0, 5);

  // Jet trees -- filled once per gen jet in each event
  float jet_corrpt = 0.;
  float jet_eta = 0.;
  float jet_csv = 0.;
  float jet_dR_leadPhoton = 0.;
  float jet_dR_trailPhoton = 0.;
  bool jet_puTight = false;
  bool jet_puMedium = false;
  bool jet_puLoose = false;
  int jet_flavor = 0;
  int jet_algDef = 0;
  int jet_mother = 0;

  TTree * jetTree = new TTree("jetTree_gg"+output_code_t, "An event tree for final analysis");
  jetTree->Branch("corrPt", &jet_corrpt, "jet_corrpt/F");
  jetTree->Branch("eta", &jet_eta, "jet_eta/F");
  jetTree->Branch("csv", &jet_csv, "jet_csv/F");
  jetTree->Branch("deltaR_leadPhoton", &jet_dR_leadPhoton, "jet_dR_leadPhoton/F");
  jetTree->Branch("deltaR_trailPhoton", &jet_dR_trailPhoton, "jet_dR_trailPhoton/F");
  jetTree->Branch("puTight", &jet_puTight, "jet_puTight/O");
  jetTree->Branch("puMedium", &jet_puMedium, "jet_puMedium/O");
  jetTree->Branch("puLoose", &jet_puLoose, "jet_puLoose/O");
  jetTree->Branch("pdgId", &jet_flavor, "jet_flavor/I");
  jetTree->Branch("algDef", &jet_algDef, "jet_algDef/I");
  jetTree->Branch("motherId", &jet_mother, "jet_mother/I");

  TTree * chsJetTree = new TTree("chsJetTree_gg"+output_code_t, "An event tree for final analysis");
  chsJetTree->Branch("corrPt", &jet_corrpt, "jet_corrpt/F");
  chsJetTree->Branch("eta", &jet_eta, "jet_eta/F");
  chsJetTree->Branch("csv", &jet_csv, "jet_csv/F");
  chsJetTree->Branch("deltaR_leadPhoton", &jet_dR_leadPhoton, "jet_dR_leadPhoton/F");
  chsJetTree->Branch("deltaR_trailPhoton", &jet_dR_trailPhoton, "jet_dR_trailPhoton/F");
  chsJetTree->Branch("puTight", &jet_puTight, "jet_puTight/O");
  chsJetTree->Branch("puMedium", &jet_puMedium, "jet_puMedium/O");
  chsJetTree->Branch("puLoose", &jet_puLoose, "jet_puLoose/O");
  chsJetTree->Branch("pdgId", &jet_flavor, "jet_flavor/I");
  chsJetTree->Branch("algDef", &jet_algDef, "jet_algDef/I");
  chsJetTree->Branch("motherId", &jet_mother, "jet_mother/I");

  TTree * eventTree = new TTree("eventInfo"+output_code_t, "ttbar info");
  // nJets, nBtags, nTight/Loose leptons, event number, pt of tight/loose leptons, mvaNonTrigV0 for electrons...
  ULong_t eventNumber_ = 0;
  int nJets, nBtags, nVetoElectrons, nTightElectrons, nVetoMuons, nTightMuons, decayMode;
  eventTree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
  eventTree->Branch("nJets", &nJets, "nJets/I");
  eventTree->Branch("nBtags", &nBtags, "nBtags/I");
  eventTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
  eventTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
  eventTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
  eventTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
  eventTree->Branch("decayMode", &decayMode, "decayMode/I");

  Float_t ele_pt, ele_eta, ele_relIso, ele_dEtaIn, ele_dPhiIn, ele_sIetaIeta, ele_hOverE, ele_d0, ele_dz, ele_fabs, ele_mvaNonTrigV0, ele_dRLeadPhoton, ele_dRTrailPhoton;
  bool ele_conversionVeto, ele_isTight, ele_isVeto;
  int ele_nMissingHits, ele_genMatchID, ele_genMatchMotherID;
  TTree * electronTree = new TTree("eleTree"+output_code_t, "electron info");
  electronTree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
  electronTree->Branch("decayMode", &decayMode, "decayMode/I");
  electronTree->Branch("pt", &ele_pt, "ele_pt/F");
  electronTree->Branch("eta", &ele_eta, "ele_eta/F");
  electronTree->Branch("relIso", &ele_relIso, "ele_relIso/F");
  electronTree->Branch("dEtaIn", &ele_dEtaIn, "ele_dEtaIn/F");
  electronTree->Branch("dPhiIn", &ele_dPhiIn, "ele_dPhiIn/F");
  electronTree->Branch("sIetaIeta", &ele_sIetaIeta, "ele_sIetaIeta/F");
  electronTree->Branch("hOverE", &ele_hOverE, "ele_hOverE/F");
  electronTree->Branch("d0", &ele_d0, "ele_d0/F");
  electronTree->Branch("dz", &ele_dz, "ele_dz/F");
  electronTree->Branch("fabs", &ele_fabs, "ele_fabs/F");
  electronTree->Branch("mvaNonTrigV0", &ele_mvaNonTrigV0, "ele_mvaNonTrigV0/F");
  electronTree->Branch("dRLeadPhoton", &ele_dRLeadPhoton, "ele_dRLeadPhoton/F");
  electronTree->Branch("dRTrailPhoton", &ele_dRTrailPhoton, "ele_dRTrailPhoton/F");
  electronTree->Branch("conversionVeto", &ele_conversionVeto, "ele_conversionVeto/O");
  electronTree->Branch("isTight", &ele_isTight, "ele_isTight/O");
  electronTree->Branch("isVeto", &ele_isVeto, "ele_isVeto/O");
  electronTree->Branch("nMissingHits", &ele_nMissingHits, "ele_nMissingHits/I");
  electronTree->Branch("genMatchID", &ele_genMatchID, "ele_genMatchID/I");
  electronTree->Branch("genMatchMotherID", &ele_genMatchMotherID, "ele_genMatchMotherID/I");
  
  ScaleFactorInfo sf(btagger);
  TFile * btagEfficiency = new TFile("btagEfficiency"+output_code_t+".root", "READ");
  sf.SetTaggingEfficiencies((TH1F*)btagEfficiency->Get("lEff"+output_code_t), (TH1F*)btagEfficiency->Get("cEff"+output_code_t), (TH1F*)btagEfficiency->Get("bEff"+output_code_t));

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

    vector<susy::Photon*> candidate_pair;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    findPhotons_prioritizeCount(event, candidate_pair, event_type, useDPhiCut);
    
    if(event_type == 0) {
      nCnt[28][0]++;
      continue;
    }

    decayMode = FigureTTbarDecayMode(event);

    float HT = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    // already ahve a gamma-gamma candidate event at this point
    // study electrons
    map<TString, vector<susy::Electron> >::iterator eleMap = event.electrons.find("gsfElectrons");
    if(eleMap != event.electrons.end()) {
      for(vector<susy::Electron>::iterator ele_it = eleMap->second.begin(); ele_it != eleMap->second.end(); ele_it++) {

	if((int)ele_it->gsfTrackIndex >= (int)(event.tracks).size() || (int)ele_it->gsfTrackIndex < 0) continue;
	if((int)ele_it->superClusterIndex >= (int)event.superClusters.size() || (int)ele_it->superClusterIndex < 0) continue;
	if(ele_it->momentum.Pt() < 2.) continue;

	ele_pt = ele_it->momentum.Pt();
	ele_eta = fabs(event.superClusters[ele_it->superClusterIndex].position.Eta());

	float ea;
	if(ele_eta < 1.0) ea = 0.13;        // ± 0.001
	else if(ele_eta < 1.479) ea = 0.14; // ± 0.002
	else if(ele_eta < 2.0) ea = 0.07;   // ± 0.001
	else if(ele_eta < 2.2) ea = 0.09;   // ± 0.001
	else if(ele_eta < 2.3) ea = 0.11;   // ± 0.002
	else if(ele_eta < 2.4) ea = 0.11;   // ± 0.003
	else ea = 0.14;                     // ± 0.004

	ele_relIso = max(0., (double)(ele_it->photonIso + ele_it->neutralHadronIso - event.rho25*ea));
	ele_relIso += ele_it->chargedHadronIso;
	ele_relIso /= ele_pt;

	ele_dEtaIn = fabs(ele_it->deltaEtaSuperClusterTrackAtVtx);
	ele_dPhiIn = fabs(ele_it->deltaPhiSuperClusterTrackAtVtx);
	ele_sIetaIeta = fabs(ele_it->sigmaIetaIeta);
	ele_hOverE = ele_it->hcalOverEcalBc;
	ele_d0 = fabs(d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
	ele_dz = fabs(dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
	ele_fabs = fabs(1/(ele_it->ecalEnergy) - 1/(ele_it->ecalEnergy/ele_it->eSuperClusterOverP));
	ele_conversionVeto = ele_it->passConversionVeto;
	ele_nMissingHits = ele_it->nMissingHits;
	 
	ele_isTight = isTightElectron(*ele_it, 
				      event.superClusters, 
				      event.rho25, 
				      d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
				      dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));

	ele_isVeto = isVetoElectron(*ele_it,
				    event.superClusters, 
				    event.rho25, 
				    d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
				    dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));

	ele_mvaNonTrigV0 = ele_it->mvaNonTrig;
	
	ele_dRLeadPhoton = deltaR(candidate_pair[0]->momentum, ele_it->momentum);
	ele_dRTrailPhoton = deltaR(candidate_pair[1]->momentum, ele_it->momentum);
	
	bool foundGenMatch = false;
	for(vector<susy::Particle>::iterator genit = event.genParticles.begin(); genit != event.genParticles.end(); genit++) {

	  if(genit->status != 1) continue;
	  if(genit->pdgId == event.genParticles[genit->motherIndex].pdgId) continue;
	  if(genit->momentum.Pt() < 2.) continue;
	  if(deltaR(ele_it->momentum, genit->momentum) >= 0.1) continue;

	  ele_genMatchID = genit->pdgId;
	  ele_genMatchMotherID = event.genParticles[genit->motherIndex].pdgId;
	  foundGenMatch = true;
	  break;
	}
    
	if(!foundGenMatch) {
	  ele_genMatchID = 0;
	  ele_genMatchMotherID = 0;
	}
	     
	electronTree->Fill();
	
      }
    }

    findMuons(event, candidate_pair, isoMuons, looseMuons, HT);
    findElectrons(event, candidate_pair, isoEles, looseEles, HT);

    map<TString, susy::PFJetCollection>::iterator pfJets_it = event.pfJets.find("ak5");
    if(pfJets_it != event.pfJets.end()) {
      susy::PFJetCollection& jetColl = pfJets_it->second;
      
      for(vector<susy::PFJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {
	
	map<TString, Float_t>::iterator s_it = it->jecScaleFactors.find("L1FastL2L3");
	float scale = s_it->second;
	
	TLorentzVector corrP4 = scale * it->momentum;
	
	bool isGood = false;

	if((it->neutralHadronEnergy/it->momentum.Energy() < 0.99) &&
	   (it->neutralEmEnergy/it->momentum.Energy() < 0.99) &&
	   ((unsigned int)it->nConstituents > 1)) {
	  
	  if(fabs(corrP4.Eta()) < 2.4) {
	    if((it->chargedHadronEnergy > 0.0) &&
	       ((int)it->chargedMultiplicity > 0) &&
	       (it->chargedEmEnergy/it->momentum.Energy() < 0.99))
	      isGood = true;
	  }
	  else isGood = true;
	}

	if(!isGood) continue;

	bool foundParent = false;
	for(vector<susy::Particle>::iterator genit = event.genParticles.begin(); genit != event.genParticles.end(); genit++) {

	  if(genit->status != 3) continue;
	  if(genit->momentum.Pt() < 20.) continue;
	  if(deltaR(corrP4, genit->momentum) > 0.01) continue;

	  bool interestingMother = 
	    fabs(event.genParticles[genit->motherIndex].pdgId) == 6 ||
	    fabs(event.genParticles[genit->motherIndex].pdgId) == 24 ||
	    fabs(event.genParticles[genit->motherIndex].pdgId) == 23;

	  bool notARepeat = 
	    fabs(genit->pdgId) != 6 &&
	    fabs(genit->pdgId) != 24 &&
	    fabs(genit->pdgId) != 23;

	  if(interestingMother && notARepeat) {
	    jet_flavor = genit->pdgId;
	    jet_mother = event.genParticles[genit->motherIndex].pdgId;
	    foundParent = true;
	    break;
	  }
	  
	}
	if(!foundParent) {
	    jet_flavor = 0;
	    jet_mother = 0;
	}
	
	jet_corrpt = corrP4.Pt();
	jet_eta = corrP4.Eta();
	jet_csv = it->bTagDiscriminators[susy::kCSV];
	jet_dR_leadPhoton = deltaR(corrP4, candidate_pair[0]->caloPosition);
	jet_dR_trailPhoton = deltaR(corrP4, candidate_pair[1]->caloPosition);
	jet_puTight = it->passPuJetIdTight(susy::kPUJetIdFull);
	jet_puMedium = it->passPuJetIdMedium(susy::kPUJetIdFull);
	jet_puLoose = it->passPuJetIdLoose(susy::kPUJetIdFull);
	jet_algDef = it->algDefFlavour;

	jetTree->Fill();
	
      } // loop over jet coll
    } // if the jet coll exists
    
    pfJets_it = event.pfJets.find("ak5chs");
    if(pfJets_it != event.pfJets.end()) {
      susy::PFJetCollection& jetColl = pfJets_it->second;
      
      for(vector<susy::PFJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {
	
	map<TString, Float_t>::iterator s_it = it->jecScaleFactors.find("L1FastL2L3");
	float scale = s_it->second;
	
	TLorentzVector corrP4 = scale * it->momentum;
	
	bool isGood = false;

	if((it->neutralHadronEnergy/it->momentum.Energy() < 0.99) &&
	   (it->neutralEmEnergy/it->momentum.Energy() < 0.99) &&
	   ((unsigned int)it->nConstituents > 1)) {
	  
	  if(fabs(corrP4.Eta()) < 2.4) {
	    if((it->chargedHadronEnergy > 0.0) &&
	       ((int)it->chargedMultiplicity > 0) &&
	       (it->chargedEmEnergy/it->momentum.Energy() < 0.99))
	      isGood = true;
	  }
	  else isGood = true;
	}

	if(!isGood) continue;

	bool foundParent = false;
	for(vector<susy::Particle>::iterator genit = event.genParticles.begin(); genit != event.genParticles.end(); genit++) {

	  if(genit->status != 3) continue;
	  if(genit->momentum.Pt() < 20.) continue;
	  if(deltaR(corrP4, genit->momentum) > 0.01) continue;

	  bool interestingMother = 
	    fabs(event.genParticles[genit->motherIndex].pdgId) == 6 ||
	    fabs(event.genParticles[genit->motherIndex].pdgId) == 24 ||
	    fabs(event.genParticles[genit->motherIndex].pdgId) == 23;

	  bool notARepeat = 
	    fabs(genit->pdgId) != 6 &&
	    fabs(genit->pdgId) != 24 &&
	    fabs(genit->pdgId) != 23;

	  if(interestingMother && notARepeat) {
	    jet_flavor = genit->pdgId;
	    jet_mother = event.genParticles[genit->motherIndex].pdgId;
	    foundParent = true;
	    break;
	  }
	  
	}
	if(!foundParent) {
	    jet_flavor = 0;
	    jet_mother = 0;
	}
	
	jet_corrpt = corrP4.Pt();
	jet_eta = corrP4.Eta();
	jet_csv = it->bTagDiscriminators[susy::kCSV];
	jet_dR_leadPhoton = deltaR(corrP4, candidate_pair[0]->caloPosition);
	jet_dR_trailPhoton = deltaR(corrP4, candidate_pair[1]->caloPosition);
	jet_puTight = it->passPuJetIdTight(susy::kPUJetIdFull);
	jet_puMedium = it->passPuJetIdMedium(susy::kPUJetIdFull);
	jet_puLoose = it->passPuJetIdLoose(susy::kPUJetIdFull);
	jet_algDef = it->algDefFlavour;

	chsJetTree->Fill();
	
      } // loop over jet coll
    } // if the jet coll exists

    findJets(event, candidate_pair, 
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT, hadronicSystem,
	     h_DR_jet_gg);

    eventNumber_ = event.eventNumber;
    nJets = pfJets.size();
    nBtags = btags.size();
    nVetoElectrons = looseEles.size();
    nTightElectrons = isoEles.size();
    nVetoMuons = looseMuons.size();
    nTightMuons = isoMuons.size();
    
    eventTree->Fill();
  
  } // for entries

  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;

  btagEfficiency->Close();

  out->Write();
  out->Close();
  */

}

void SusyEventAnalyzer::SignalContent_gg() {
  /*
  char * tmp = getenv("CONDOR_SECTION");
  cout << "tmp = " << tmp << endl;
  int index = atoi(tmp);
  cout << "index = " << index << endl;
  
  // stop-bino scan
  int index1 = mst[int(index)/31];
  int index2 = mBino[int(index)%31];
  
  // bino/wino
  //int index1 = int(index)/17*100 + 400;
  //int index2 = int(index)%17*100 + 420;
  
  cout << "index1 = " << index1 << endl;
  cout << "index2 = " << index2 << endl;
  
  char output_file_name[100];
  
  ScaleFactorInfo sf(btagger);

  sprintf(output_file_name, "signalContent_mst_%d_m1_%d.root", index1, index2);

  TFile * out = new TFile(output_file_name, "RECREATE");
  out->cd();

  sprintf(output_file_name, "_mst_%d_m1_%d", index1, index2);
  TString code = output_file_name;

  TH2F * elePt_leadPt = new TH2F("elePt_leadPt"+code, "Electron pt vs leading photon pt", 400, 0, 2000, 400, 0, 2000);
  TH2F * elePt_subPt = new TH2F("elePt_subPt"+code, "Electron pt vs sub-leading photon pt", 400, 0, 2000, 400, 0, 2000);

  TH1F * h_w_mass = new TH1F("w_mass"+code, "jj mass", 1000, 0, 2000);
  TH1F * h_top_mass = new TH1F("top_mass"+code, "bjj mass", 1000, 0, 2000);

  // branching counts
  int nStops, tFromStop, bFromStop, wFromStop, nBinos;
  int nLeptonicW, nHadronicW, nLeptonicZ, nHadronicZ, nInvisibleZ;
  int gFromBino, zFromBino, wFromBino;

  // object counts
  int nPhotons_80, nPhotons_40, nPhotons_25, nPhotons;
  int nMuons, nElectrons;
  int nJets, nBjets;
  
  TTree * genTree = new TTree("genInfo"+code, "gen-level information tree");
  genTree->Branch("nStops", &nStops, "nStops/I");
  genTree->Branch("tFromStop", &tFromStop, "tFromStop/I");
  genTree->Branch("bFromStop", &bFromStop, "bFromStop/I");
  genTree->Branch("wFromStop", &wFromStop, "wFromStop/I");
  genTree->Branch("nBinos", &nBinos, "nBinos/I");
  genTree->Branch("nLeptonicW", &nLeptonicW, "nLeptonicW/I");
  genTree->Branch("nHadronicW", &nHadronicW, "nHadronicW/I");
  genTree->Branch("nLeptonicZ", &nLeptonicZ, "nLeptonicZ/I");
  genTree->Branch("nHadronicZ", &nHadronicZ, "nHadronicZ/I");
  genTree->Branch("nInvisibleZ", &nInvisibleZ, "nInvisibleZ/I");
  genTree->Branch("gFromBino", &gFromBino, "gFromBino/I");
  genTree->Branch("zFromBino", &zFromBino, "zFromBino/I");
  genTree->Branch("wFromBino", &wFromBino, "wFromBino/I");
  genTree->Branch("nPhotons_80", &nPhotons_80, "nPhotons_80/I");
  genTree->Branch("nPhotons_40", &nPhotons_40, "nPhotons_40/I");
  genTree->Branch("nPhotons_25", &nPhotons_25, "nPhotons_25/I");
  genTree->Branch("nPhotons", &nPhotons, "nPhotons/I");
  genTree->Branch("nMuons", &nMuons, "nMuons/I");
  genTree->Branch("nElectrons", &nElectrons, "nElectrons/I");
  genTree->Branch("nJets", &nJets, "nJets/I");
  genTree->Branch("nBjets", &nBjets, "nBjets/I");
  
  Float_t ele_pt, ele_eta, ele_relIso, ele_dEtaIn, ele_dPhiIn, ele_sIetaIeta, ele_hOverE, ele_d0, ele_dz, ele_fabs, ele_mvaNonTrigV0, ele_dRLeadPhoton, ele_dRTrailPhoton;
  bool ele_conversionVeto, ele_isTight, ele_isVeto;
  int ele_nMissingHits, ele_genMatchID, ele_genMatchMotherID, decayMode;
  ULong_t eventNumber_ = 0;
  TTree * electronTree = new TTree("eleTree", "electron info");
  electronTree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
  electronTree->Branch("decayMode", &decayMode, "decayMode/I");
  electronTree->Branch("pt", &ele_pt, "ele_pt/F");
  electronTree->Branch("eta", &ele_eta, "ele_eta/F");
  electronTree->Branch("relIso", &ele_relIso, "ele_relIso/F");
  electronTree->Branch("dEtaIn", &ele_dEtaIn, "ele_dEtaIn/F");
  electronTree->Branch("dPhiIn", &ele_dPhiIn, "ele_dPhiIn/F");
  electronTree->Branch("sIetaIeta", &ele_sIetaIeta, "ele_sIetaIeta/F");
  electronTree->Branch("hOverE", &ele_hOverE, "ele_hOverE/F");
  electronTree->Branch("d0", &ele_d0, "ele_d0/F");
  electronTree->Branch("dz", &ele_dz, "ele_dz/F");
  electronTree->Branch("fabs", &ele_fabs, "ele_fabs/F");
  electronTree->Branch("mvaNonTrigV0", &ele_mvaNonTrigV0, "ele_mvaNonTrigV0/F");
  electronTree->Branch("dRLeadPhoton", &ele_dRLeadPhoton, "ele_dRLeadPhoton/F");
  electronTree->Branch("dRTrailPhoton", &ele_dRTrailPhoton, "ele_dRTrailPhoton/F");
  electronTree->Branch("conversionVeto", &ele_conversionVeto, "ele_conversionVeto/O");
  electronTree->Branch("isTight", &ele_isTight, "ele_isTight/O");
  electronTree->Branch("isVeto", &ele_isVeto, "ele_isVeto/O");
  electronTree->Branch("nMissingHits", &ele_nMissingHits, "ele_nMissingHits/I");
  electronTree->Branch("genMatchID", &ele_genMatchID, "ele_genMatchID/I");
  electronTree->Branch("genMatchMotherID", &ele_genMatchMotherID, "ele_genMatchMotherID/I");

  float mu_pt, mu_eta, mu_iso, mu_relIso;
  TTree * muTree = new TTree("muTree"+code, "gen-matched muon reco info");
  muTree->Branch("pt", &mu_pt, "mu_pt/F");
  muTree->Branch("eta", &mu_eta, "mu_eta/F");
  muTree->Branch("iso", &mu_iso, "mu_iso/F");
  muTree->Branch("relIso", &mu_relIso, "mu_relIso/F");

  float g_pt, g_eta;
  TTree * photonTree = new TTree("photonTree"+code, "gen-matched photon reco info");
  photonTree->Branch("pt", &g_pt, "g_pt/F");
  photonTree->Branch("eta", &g_eta, "g_eta/F");

  float jet_pt, jet_eta, jet_csv, jet_jp, jet_tchp;
  int jet_flavor, jet_algDef, jet_phyDef, jet_mother;
  bool jet_pujid_lF, jet_pujid_lS, jet_pujid_lC, jet_pujid_mF, jet_pujid_mS, jet_pujid_mC, jet_pujid_tF, jet_pujid_tS, jet_pujid_tC;
  TTree * jetTree = new TTree("jetTree"+code, "gen-matched jet reco info");
  jetTree->Branch("pt", &jet_pt, "jet_pt/F");
  jetTree->Branch("eta", &jet_eta, "jet_eta/F");
  jetTree->Branch("csv", &jet_csv, "jet_csv/F");
  jetTree->Branch("jp", &jet_jp, "jet_jp/F");
  jetTree->Branch("tchp", &jet_tchp, "jet_tchp/F");
  jetTree->Branch("flavor", &jet_flavor, "jet_flavor/I");
  jetTree->Branch("algDef", &jet_algDef, "jet_algDef/I");
  jetTree->Branch("phyDef", &jet_phyDef, "jet_phyDef/I");
  jetTree->Branch("mother", &jet_mother, "jet_mother/I");
  jetTree->Branch("puJetId_lF", &jet_pujid_lF, "jet_pujid_lF/O");
  jetTree->Branch("puJetId_lS", &jet_pujid_lS, "jet_pujid_lS/O");
  jetTree->Branch("puJetId_lC", &jet_pujid_lC, "jet_pujid_lC/O");
  jetTree->Branch("puJetId_mF", &jet_pujid_mF, "jet_pujid_mF/O");
  jetTree->Branch("puJetId_mS", &jet_pujid_mS, "jet_pujid_mS/O");
  jetTree->Branch("puJetId_mC", &jet_pujid_mC, "jet_pujid_mC/O");
  jetTree->Branch("puJetId_tF", &jet_pujid_tF, "jet_pujid_tF/O");
  jetTree->Branch("puJetId_tS", &jet_pujid_tS, "jet_pujid_tS/O");
  jetTree->Branch("puJetId_tC", &jet_pujid_tC, "jet_pujid_tC/O");

  float btag_pt, btag_eta;
  int btag_flavor, btag_algDef, btag_phyDef, btag_mother;
  bool btag_pujid_lF, btag_pujid_lS, btag_pujid_lC, btag_pujid_mF, btag_pujid_mS, btag_pujid_mC, btag_pujid_tF, btag_pujid_tS, btag_pujid_tC;
  TTree * btagTree = new TTree("btagTree"+code, "gen-matched btag reco info");
  btagTree->Branch("pt", &btag_pt, "btag_pt/F");
  btagTree->Branch("eta", &btag_eta, "btag_eta/F");
  btagTree->Branch("flavor", &btag_flavor, "btag_flavor/I");
  btagTree->Branch("algDef", &btag_algDef, "btag_algDef/I");
  btagTree->Branch("phyDef", &btag_phyDef, "btag_phyDef/I");
  btagTree->Branch("mother", &btag_mother, "btag_mother/I");
  btagTree->Branch("puJetId_lF", &btag_pujid_lF, "btag_pujid_lF/O");
  btagTree->Branch("puJetId_lS", &btag_pujid_lS, "btag_pujid_lS/O");
  btagTree->Branch("puJetId_lC", &btag_pujid_lC, "btag_pujid_lC/O");
  btagTree->Branch("puJetId_mF", &btag_pujid_mF, "btag_pujid_mF/O");
  btagTree->Branch("puJetId_mS", &btag_pujid_mS, "btag_pujid_mS/O");
  btagTree->Branch("puJetId_mC", &btag_pujid_mC, "btag_pujid_mC/O");
  btagTree->Branch("puJetId_tF", &btag_pujid_tF, "btag_pujid_tF/O");
  btagTree->Branch("puJetId_tS", &btag_pujid_tS, "btag_pujid_tS/O");
  btagTree->Branch("puJetId_tC", &btag_pujid_tC, "btag_pujid_tC/O");
  
  float maxCSV, maxJP, maxTCHP;
  float submaxCSV, submaxJP, submaxTCHP;
  bool isDiPhoton, isSinglePhoton;
  bool isGGL, isGL;
  TTree * eventTree = new TTree("eventTree"+code, "event-wide info");
  eventTree->Branch("maxCSV", &maxCSV, "maxCSV/F");
  eventTree->Branch("submaxCSV", &submaxCSV, "submaxCSV/F");
  eventTree->Branch("maxJP", &maxJP, "maxJP/F");
  eventTree->Branch("submaxJP", &submaxJP, "submaxJP/F");
  eventTree->Branch("maxTCHP", &maxTCHP, "maxTCHP/F");
  eventTree->Branch("submaxTCHP", &submaxTCHP, "submaxTCHP/F");
  eventTree->Branch("isDiPhoton", &isDiPhoton, "isDiPhoton/O");
  eventTree->Branch("isSinglePhoton", &isSinglePhoton, "isSinglePhoton/O");
  eventTree->Branch("isGGL", &isGGL, "isGGL/O");
  eventTree->Branch("isGL", &isGL, "isGL/O");
  
  vector<TString> triggerNames_;
  vector<int> triggerCounts_;

  vector<pair<TString, int> > triggerFires;

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;
  
  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = "
	   << event.runNumber << ", event = " << event.eventNumber << ", integ lumi = " 
	   << event.intgRecLumi << endl;
    }

    for(susy::TriggerMap::iterator it = event.hltMap.begin(); it != event.hltMap.end(); it++) {
      if((int(it->second.second)) && it->first.Contains("HLT_") && it->second.first == 1) {

	bool alreadyFound = false;
	for(unsigned int k = 0; k < triggerFires.size(); k++) {
	  if((triggerFires[k].first).Contains(it->first)) {
	    triggerFires[k].second++;
	    alreadyFound = true;
	    break;
	  }
	}

	if(!alreadyFound) {
	  triggerFires.push_back(make_pair(it->first, 1));
	}
	
	
      }
    }



    nStops = tFromStop = bFromStop = nBinos = 0;
    nLeptonicW = nHadronicW = nLeptonicZ = nHadronicZ = nInvisibleZ = 0;
    gFromBino = zFromBino = wFromBino = 0;
    nPhotons_80 = nPhotons_40 = nPhotons_25 = nPhotons = 0;
    nMuons = nElectrons = 0;
    nJets = nBjets = 0;
    mu_pt = mu_eta = mu_iso = mu_relIso = 0;
    g_pt = g_eta = 0;
    jet_pt = jet_eta = jet_csv = jet_jp = jet_tchp = 0;
    jet_flavor = jet_algDef = jet_phyDef = jet_mother = 0;
    jet_pujid_lF = jet_pujid_lS = jet_pujid_lC = jet_pujid_mF = jet_pujid_mS = jet_pujid_mC = jet_pujid_tF = jet_pujid_tS = jet_pujid_tC = false;
    btag_pujid_lF = btag_pujid_lS = btag_pujid_lC = btag_pujid_mF = btag_pujid_mS = btag_pujid_mC = btag_pujid_tF = btag_pujid_tS = btag_pujid_tC = false;
    btag_pt = btag_eta = 0;
    btag_flavor = btag_algDef = btag_phyDef = btag_mother = 0;
    maxCSV = maxJP = maxTCHP = 0;
    submaxCSV = submaxJP = submaxTCHP = 0;
    isDiPhoton = false;
    isSinglePhoton = false;
    isGGL = false;
    isGL = false;

    vector<double> pfJets_CSV, pfJets_JP, pfJets_TCHP;

    vector<TLorentzVector> bMomenta, jMomenta;

    map<TString, susy::PFJetCollection>::iterator pfJets_it = event.pfJets.find("ak5");
    if(pfJets_it == event.pfJets.end()) cout << "Wat @ find ak5" << endl;

    vector<float> v_photon_pt;

    for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

      if(fabs(it->pdgId) == 1000006 && it->status == 3) nStops++;
      if(fabs(it->pdgId) == 1000022 && it->status == 3) nBinos++;
      if(fabs(it->pdgId) == 6 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000006 && it->status == 3) tFromStop++;
      if(fabs(it->pdgId) == 5 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000006 && it->status == 3) bFromStop++;
      if(fabs(it->pdgId) == 24 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000006 && it->status == 3) wFromStop++;
      if(fabs(it->pdgId) == 22 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && it->status == 1) gFromBino++;
      if(fabs(it->pdgId) == 23 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && it->status == 3) zFromBino++;
      if(fabs(it->pdgId) == 24 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && it->status == 3) wFromBino++;

      // Gen-level photons from binos in acceptance
      if(fabs(it->pdgId) == 22 && it->status == 1 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && fabs(it->momentum.Eta()) < 1.4442) {

	  if(it->momentum.Pt() > 25.0) nPhotons_25++;
	  if(it->momentum.Pt() > 40.0) nPhotons_40++;
	  if(it->momentum.Pt() > 80.0) nPhotons_80++;
	  nPhotons++;

	// find photons
	map<TString, vector<susy::Photon> >::iterator phoMap = event.photons.find("photons");
	if(phoMap != event.photons.end()) {
	  
	  for(vector<susy::Photon>::iterator pho_it = phoMap->second.begin();
	      pho_it != phoMap->second.end(); pho_it++) {
	    
	    if(deltaR(pho_it->momentum, it->momentum) >= 0.3) continue;

	    g_pt = pho_it->momentum.Pt();
	    g_eta = pho_it->caloPosition.Eta();

	    v_photon_pt.push_back(g_pt);

	    photonTree->Fill();
	    break;
	    
	  } // for photon
	} // if
    
      }

      float leading_photon_pt = -1.;
      float subleading_photon_pt = -1.;
      sort(v_photon_pt.begin(), v_photon_pt.end(), greater<float>());
      if(v_photon_pt.size() >= 1) leading_photon_pt = v_photon_pt[0];
      if(v_photon_pt.size() >= 2) subleading_photon_pt = v_photon_pt[1];

      bool isFromWfromStopOrTop = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
	 (
	  fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 6 || 
	  fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000006
	  );

      bool isFromWorZfromBino = (fabs(event.genParticles[it->motherIndex].pdgId) == 23 || fabs(event.genParticles[it->motherIndex].pdgId) == 24) &&
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000022;

      if((isFromWfromStopOrTop || isFromWorZfromBino) && 
	 fabs(it->momentum.Eta()) < 2.6 &&
	 !(fabs(it->momentum.Eta()) >= 1.4442 && fabs(it->momentum.Eta()) <= 1.566)) {

	if(fabs(it->pdgId) == 11 ||
	   fabs(it->pdgId) == 13) {

	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 24) nLeptonicW++;
	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 23 && nLeptonicZ == 0) nLeptonicZ++;

	  if(fabs(it->pdgId) == 13) {

	    nMuons++;

	    // lepton from stop->top->W; look at kinematics of reconstructed lepton
	    map<TString, vector<susy::Muon> >::iterator muMap = event.muons.find("muons");
	    if(muMap != event.muons.end()) {
	      for(vector<susy::Muon>::iterator mu_it = muMap->second.begin(); mu_it != muMap->second.end(); mu_it++) {
	
		if(deltaR(mu_it->momentum, it->momentum) >= 0.3) continue;

		mu_pt = mu_it->momentum.Pt();
		mu_eta = mu_it->momentum.Eta();
		mu_iso = max(0., (mu_it->sumNeutralHadronEt04 + mu_it->sumPhotonEt04 - 0.5*(mu_it->sumPUPt04)));
		mu_iso += mu_it->sumChargedHadronPt04;
		mu_relIso = mu_iso / mu_pt;
		
		muTree->Fill();
		break;

	      }
	    }
	  } // if muon

	  if(fabs(it->pdgId) == 11) {

	    nElectrons++;

	    map<TString, vector<susy::Electron> >::iterator eleMap = event.electrons.find("gsfElectrons");
	    if(eleMap != event.electrons.end()) {
	      for(vector<susy::Electron>::iterator ele_it = eleMap->second.begin(); ele_it != eleMap->second.end(); ele_it++) {

		if(deltaR(ele_it->momentum, it->momentum) >= 0.1) continue;
		if((int)ele_it->gsfTrackIndex >= (int)(event.tracks).size() || (int)ele_it->gsfTrackIndex < 0) continue;
		if((int)ele_it->superClusterIndex >= (int)event.superClusters.size() || (int)ele_it->superClusterIndex < 0) continue;

		ele_pt = ele_it->momentum.Pt();
		ele_eta = fabs(event.superClusters[ele_it->superClusterIndex].position.Eta());

		float ea;
		if(ele_eta < 1.0) ea = 0.13;        // ± 0.001
		else if(ele_eta < 1.479) ea = 0.14; // ± 0.002
		else if(ele_eta < 2.0) ea = 0.07;   // ± 0.001
		else if(ele_eta < 2.2) ea = 0.09;   // ± 0.001
		else if(ele_eta < 2.3) ea = 0.11;   // ± 0.002
		else if(ele_eta < 2.4) ea = 0.11;   // ± 0.003
		else ea = 0.14;                     // ± 0.004
		
		ele_relIso = max(0., (double)(ele_it->photonIso + ele_it->neutralHadronIso - event.rho25*ea));
		ele_relIso += ele_it->chargedHadronIso;
		ele_relIso /= ele_pt;
		
		ele_dEtaIn = fabs(ele_it->deltaEtaSuperClusterTrackAtVtx);
		ele_dPhiIn = fabs(ele_it->deltaPhiSuperClusterTrackAtVtx);
		ele_sIetaIeta = fabs(ele_it->sigmaIetaIeta);
		ele_hOverE = ele_it->hcalOverEcalBc;
		ele_d0 = fabs(d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
		ele_dz = fabs(dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
		ele_fabs = fabs(1/(ele_it->ecalEnergy) - 1/(ele_it->ecalEnergy/ele_it->eSuperClusterOverP));
		ele_conversionVeto = ele_it->passConversionVeto;
		ele_nMissingHits = ele_it->nMissingHits;
		
		ele_isTight = isTightElectron(*ele_it, 
					      event.superClusters, 
					      event.rho25, 
					      d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
					      dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
		
		ele_isVeto = isVetoElectron(*ele_it,
					    event.superClusters, 
					    event.rho25, 
					    d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
					    dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));
		
		ele_mvaNonTrigV0 = ele_it->mvaNonTrig;

		ele_genMatchID = it->pdgId;
		ele_genMatchMotherID = event.genParticles[it->motherIndex].pdgId;

		decayMode = FigureTTbarDecayMode(event);
		
		electronTree->Fill();
	
		break;
	      }
	      
	    }
	  } // if electron
	    
	} // if lepton from stop->top->W
	else if(fabs(it->pdgId) != 12 &&
		fabs(it->pdgId) != 14 &&
		fabs(it->pdgId) != 16) {
	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 24) nHadronicW++;
	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 23 && nHadronicZ == 0) nHadronicZ++;
	}
	else if(fabs(it->pdgId) == 12 &&
		fabs(it->pdgId) == 14 &&
		fabs(it->pdgId) == 16 &&
		nInvisibleZ == 0) nInvisibleZ++;

      }

      bool jetFromCascade = fabs(event.genParticles[it->motherIndex].pdgId) >= 1000001 && 
	fabs(event.genParticles[it->motherIndex].pdgId) <= 1000038;

      bool jetFromWfromTop = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 6;

      bool jetFromTop = fabs(event.genParticles[it->motherIndex].pdgId) == 6;

      bool jetFromWfromBino = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000022;

      bool jetFromZfromBino = fabs(event.genParticles[it->motherIndex].pdgId) == 23 && 
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000022;

      bool jetFromInteraction = jetFromCascade || jetFromWfromTop || jetFromTop || jetFromWfromBino || jetFromZfromBino;

      // Gen-level jets
      if((fabs(it->pdgId) == 1 ||
	  fabs(it->pdgId) == 2 ||
	  fabs(it->pdgId) == 3 ||
	  fabs(it->pdgId) == 4 ||
	  fabs(it->pdgId) == 5 ||
	  fabs(it->pdgId) == 21) &&
	 //it->status == 3 &&
	 jetFromInteraction &&
	 it->momentum.Pt() > 30.0 &&
	 fabs(it->momentum.Eta()) < 2.6) {

	nJets++;

	if(fabs(it->pdgId) == 5 && fabs(it->momentum.Eta()) < 2.4) {
	  nBjets++;
	}

	// Match this jet to a good reco pfjet
	if(pfJets_it != event.pfJets.end()) {
	  for(vector<susy::PFJet>::iterator jet = pfJets_it->second.begin(); jet != pfJets_it->second.end(); jet++) {
	    map<TString, Float_t>::iterator s_it = jet->jecScaleFactors.find("L1FastL2L3");
	    if(s_it == jet->jecScaleFactors.end()) { cout << "Wat @ jec" << endl; continue; }
	    TLorentzVector corrP4 = s_it->second * jet->momentum;
	    
	    if(deltaR(corrP4, it->momentum) >= 0.3) continue;

	    jet_pt = corrP4.Pt();
	    jet_eta = corrP4.Eta();
	    jet_csv = jet->bTagDiscriminators[susy::kCSV];
	    jet_jp = jet->bTagDiscriminators[susy::kJP];
	    jet_tchp = jet->bTagDiscriminators[susy::kTCHP];

	    jet_flavor = it->pdgId;
	    jet_algDef = jet->algDefFlavour;
	    jet_phyDef = jet->phyDefFlavour;
	    jet_mother = event.genParticles[it->motherIndex].pdgId;

	    jet_pujid_lF = jet->passPuJetIdLoose(susy::kPUJetIdFull);
	    jet_pujid_lS = jet->passPuJetIdLoose(susy::kPUJetIdSimple);
	    jet_pujid_lC = jet->passPuJetIdLoose(susy::kPUJetIdCutBased);

	    jet_pujid_mF = jet->passPuJetIdMedium(susy::kPUJetIdFull);
	    jet_pujid_mS = jet->passPuJetIdMedium(susy::kPUJetIdSimple);
	    jet_pujid_mC = jet->passPuJetIdMedium(susy::kPUJetIdCutBased);

	    jet_pujid_tF = jet->passPuJetIdTight(susy::kPUJetIdFull);
	    jet_pujid_tS = jet->passPuJetIdTight(susy::kPUJetIdSimple);
	    jet_pujid_tC = jet->passPuJetIdTight(susy::kPUJetIdCutBased);

	    jetTree->Fill();

	    bool is_btagged = false;
	    
	    if(fabs(corrP4.Eta()) < 2.4) {
	      pfJets_CSV.push_back(jet_csv);
	      pfJets_JP.push_back(jet_jp);
	      pfJets_TCHP.push_back(jet_tchp);
	   
	      if((btagger == "CSVL" && jet->bTagDiscriminators[susy::kCSV] > 0.244) ||
		 (btagger == "CSVM" && jet->bTagDiscriminators[susy::kCSV] > 0.679) ||
		 (btagger == "CSVT" && jet->bTagDiscriminators[susy::kCSV] > 0.898)) {

		btag_pt = jet_pt;
		btag_eta = jet_eta;

		btag_flavor = jet_flavor;
		btag_algDef = jet_algDef;
		btag_phyDef = jet_phyDef;
		btag_mother = jet_mother;

		btag_pujid_lF = jet->passPuJetIdLoose(susy::kPUJetIdFull);
		btag_pujid_lS = jet->passPuJetIdLoose(susy::kPUJetIdSimple);
		btag_pujid_lC = jet->passPuJetIdLoose(susy::kPUJetIdCutBased);

		btag_pujid_mF = jet->passPuJetIdMedium(susy::kPUJetIdFull);
		btag_pujid_mS = jet->passPuJetIdMedium(susy::kPUJetIdSimple);
		btag_pujid_mC = jet->passPuJetIdMedium(susy::kPUJetIdCutBased);

		btag_pujid_tF = jet->passPuJetIdTight(susy::kPUJetIdFull);
		btag_pujid_tS = jet->passPuJetIdTight(susy::kPUJetIdSimple);
		btag_pujid_tC = jet->passPuJetIdTight(susy::kPUJetIdCutBased);

		btagTree->Fill();

		bMomenta.push_back(corrP4);
		is_btagged = true;
	      }
	      
	    }

	    if(!is_btagged)  jMomenta.push_back(corrP4);

	    break;
	  }
	}

      }

    }

    sort(pfJets_CSV.begin(), pfJets_CSV.end(), greater<double>());
    if(pfJets_CSV.size() >= 1) maxCSV = pfJets_CSV[0];
    else maxCSV = -9999.;
    if(pfJets_CSV.size() >= 2) submaxCSV = pfJets_CSV[1];
    else submaxCSV = -9999.;

    sort(pfJets_JP.begin(), pfJets_JP.end(), greater<double>());
    if(pfJets_JP.size() >= 1) maxJP = pfJets_JP[0];
    else maxJP = -9999.;
    if(pfJets_JP.size() >= 2) submaxJP = pfJets_JP[1];
    else submaxJP = -9999.;

    sort(pfJets_TCHP.begin(), pfJets_TCHP.end(), greater<double>());
    if(pfJets_TCHP.size() >= 1) maxTCHP = pfJets_TCHP[0];
    else maxTCHP = -9999.;
    if(pfJets_TCHP.size() >= 2) submaxTCHP = pfJets_TCHP[1];
    else submaxTCHP = -9999.;

    isDiPhoton = nPhotons_25 > 1 && nPhotons_40 > 0;
    isSinglePhoton = nPhotons_80 > 0;

    isGL = isSinglePhoton && (nMuons > 0 || nElectrons > 0);
    isGGL = isDiPhoton && (nMuons > 0 || nElectrons > 0);

    genTree->Fill();
    eventTree->Fill();

    if(bMomenta.size() >= 1 && jMomenta.size() >= 2) {
      for(unsigned int ib = 0; ib < bMomenta.size(); ib++) {
	for(unsigned int j1 = 0; j1 < jMomenta.size() - 1; j1++) {
	  for(unsigned int j2 = j1; j2 < jMomenta.size(); j2++) {

	    double w_mass = (jMomenta[j1] + jMomenta[j2]).M();
	    h_w_mass->Fill(w_mass);

	    double top_mass = (bMomenta[ib] + jMomenta[j1] + jMomenta[j2]).M();
	    h_top_mass->Fill(top_mass);

	  }
	}
      }

    }
	 

  } // event loop

  cout << endl << endl << "Found trigger counts: " << endl << endl;

  sort(triggerFires.begin(), triggerFires.end(), sortTriggers);
  for(unsigned int i = 0; i < triggerFires.size(); i++) {
    cout << triggerFires[i].first << " -- " << triggerFires[i].second << " times" << endl;
  }

  out->cd();
  out->Write();
  out->Close();

  */
}

void SusyEventAnalyzer::PhotonInfo() {
  /*
  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
    nCnt[i][j] = 0;
    }
  }

  TFile* out = new TFile("photonInfo_"+outputName+"_"+btagger+".root", "RECREATE");
  out->cd();

  float eta, et, phi, hOverE, neutralHadIso, photonIso, chargedHadIso, r9, sIetaIeta, sIphiIphi, worstOtherVtxChargedHadronIso, MVAregEnergy;
  int nPixelSeeds;
  bool isGamma, isFake, isElectron, hasJetMatch;

  float jet_dR, jet_dPhi;

  float chargedHadronEnergy, neutralHadronEnergy, photonEnergy, electronEnergy, muonEnergy, HFHadronEnergy, HFEMEnergy, chargedEmEnergy, chargedMuEnergy, neutralEmEnergy;
  int chargedHadronMultiplicity, neutralHadronMultiplicity, photonMultiplicity, electronMultiplicity, muonMultiplicity, HFHadronMultiplicity, HFEMMultiplicity, chargedMultiplicity, neutralMultiplicity;
    
  TTree * photonTree = new TTree("photonTree", "photon info");
  photonTree->Branch("eta", &eta, "eta/F");
  photonTree->Branch("et", &et, "et/F");
  photonTree->Branch("phi", &phi, "phi/F");
  photonTree->Branch("hOverE", &hOverE, "hOverE/F");
  photonTree->Branch("neutralHadIso", &neutralHadIso, "neutralHadIso/F");
  photonTree->Branch("photonIso", &photonIso, "photonIso/F");
  photonTree->Branch("chargedHadIso", &chargedHadIso, "chargedHadIso/F");
  photonTree->Branch("r9", &r9, "r9/F");
  photonTree->Branch("sIetaIeta", &sIetaIeta, "sIetaIeta/F");
  photonTree->Branch("sIphiIphi", &sIphiIphi, "sIphiIphi/F");
  photonTree->Branch("worstOtherVtxChargedHadronIso", &worstOtherVtxChargedHadronIso, "worstOtherVtxChargedHadronIso/F");
  photonTree->Branch("MVAregEnergy", &MVAregEnergy, "MVAregEnergy/F");
  photonTree->Branch("nPixelSeeds", &nPixelSeeds, "nPixelSeeds/I");
  photonTree->Branch("isGamma", &isGamma, "isGamma/O");
  photonTree->Branch("isFake", &isFake, "isFake/O");
  photonTree->Branch("isElectron", &isElectron, "isElectron/O");
  photonTree->Branch("hasJetMatch", &hasJetMatch, "hasJetMatch/O");
  photonTree->Branch("jet_dR", &jet_dR, "jet_dR/F");
  photonTree->Branch("jet_dPhi", &jet_dPhi, "jet_dPhi/F");

  photonTree->Branch("jet_chargedHadronEnergy", &chargedHadronEnergy, "chargedHadronEnergy/F");
  photonTree->Branch("jet_neutralHadronEnergy", &neutralHadronEnergy, "neutralHadronEnergy/F");
  photonTree->Branch("jet_photonEnergy", &photonEnergy, "photonEnergy/F");
  photonTree->Branch("jet_electronEnergy", &electronEnergy, "electronEnergy/F");
  photonTree->Branch("jet_muonEnergy", &muonEnergy, "muonEnergy/F");
  photonTree->Branch("jet_HFHadronEnergy", &HFHadronEnergy, "HFHadronEnergy/F");
  photonTree->Branch("jet_HFEMEnergy", &HFEMEnergy, "HFEMEnergy/F");
  photonTree->Branch("jet_chargedEmEnergy", &chargedEmEnergy, "chargedEmEnergy/F");
  photonTree->Branch("jet_chargedMuEnergy", &chargedMuEnergy, "chargedMuEnergy/F");
  photonTree->Branch("jet_neutralEmEnergy", &neutralEmEnergy, "neutralEmEnergy/F");
  photonTree->Branch("jet_chargedHadronMultiplicity", &chargedHadronMultiplicity, "chargedHadronMultiplicity/I");
  photonTree->Branch("jet_neutralHadronMultiplicity", &neutralHadronMultiplicity, "neutralHadronMultiplicity/I");
  photonTree->Branch("jet_photonMultiplicity", &photonMultiplicity, "photonMultiplicity/I");
  photonTree->Branch("jet_electronMultiplicity", &electronMultiplicity, "electronMultiplicity/I");
  photonTree->Branch("jet_muonMultiplicity", &muonMultiplicity, "muonMultiplicity/I");
  photonTree->Branch("jet_HFHadronMultiplicity", &HFHadronMultiplicity, "HFHadronMultiplicity/I");
  photonTree->Branch("jet_HFEMMultiplicity", &HFEMMultiplicity, "HFEMMultiplicity/I");
  photonTree->Branch("jet_chargedMultiplicity", &chargedMultiplicity, "chargedMultiplicity/I");
  photonTree->Branch("jet_neutralMultiplicity", &neutralMultiplicity, "neutralMultiplicity/I");

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

    vector<susy::Photon*> candidate_pair;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    findPhotons_prioritizeCount(event, candidate_pair, event_type, useDPhiCut);

    if(event_type == 0) {
      nCnt[28][0]++;
      continue;
    }

    bool passHLT = useTrigger ? PassTriggers(abs(event_type)) : true;
    if(!passHLT) {
      const int nPos = 30 + abs(event_type);
      nCnt[nPos][0]++;
      continue;
    }

    map<TString, susy::PFJetCollection>::iterator iJets = event.pfJets.find("ak5");
    susy::PFJetCollection& jetColl = iJets->second;

    for(unsigned int i = 0; i < candidate_pair.size(); i++) {
      
      eta = candidate_pair[i]->caloPosition.Eta();
      phi = candidate_pair[i]->caloPosition.Phi();
      et = candidate_pair[i]->momentum.Et();
      hOverE = candidate_pair[i]->hadTowOverEm;
      neutralHadIso = neutralHadronIso_corrected(*candidate_pair[i], event.rho25);
      photonIso = photonIso_corrected(*candidate_pair[i], event.rho25);
      chargedHadIso = chargedHadronIso_corrected(*candidate_pair[i], event.rho25);
      r9 = candidate_pair[i]->r9;
      sIetaIeta = candidate_pair[i]->sigmaIetaIeta;
      sIphiIphi = candidate_pair[i]->sigmaIphiIphi;
      worstOtherVtxChargedHadronIso = candidate_pair[i]->worstOtherVtxChargedHadronIso;
      MVAregEnergy = candidate_pair[i]->MVAregEnergy;
      nPixelSeeds = candidate_pair[i]->nPixelSeeds;
      isGamma = is_eg(*candidate_pair[i], event.rho25) && nPixelSeeds == 0;
      isElectron = is_eg(*candidate_pair[i], event.rho25) && nPixelSeeds != 0;
      isFake = is_f(*candidate_pair[i], event.rho25) && nPixelSeeds == 0;

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
	
	float dEta = corrP4.Eta() - candidate_pair[i]->caloPosition.Eta();
	float dPhi = TVector2::Phi_mpi_pi(corrP4.Phi() - candidate_pair[i]->caloPosition.Phi());
	float dR = sqrt(dEta*dEta + dPhi*dPhi);
	
	if(corrP4.Et() > 20. &&
	   fabs(corrP4.Eta()) < 2.6 &&
	   dR < 0.3) {
	  
	  matched = true;

	  hasJetMatch = true;
	  jet_dR = dR;
	  jet_dPhi = dPhi;

	  chargedHadronEnergy = iJet->chargedHadronEnergy;
	  neutralHadronEnergy = iJet->neutralHadronEnergy;
	  photonEnergy = iJet->photonEnergy;
	  electronEnergy = iJet->electronEnergy;
	  muonEnergy = iJet->muonEnergy;
	  HFHadronEnergy = iJet->HFHadronEnergy;
	  HFEMEnergy = iJet->HFEMEnergy;
	  chargedEmEnergy = iJet->chargedEmEnergy;
	  chargedMuEnergy = iJet->chargedMuEnergy;
	  neutralEmEnergy = iJet->neutralEmEnergy;
	  chargedHadronMultiplicity = iJet->chargedHadronMultiplicity;
	  neutralHadronMultiplicity = iJet->neutralHadronMultiplicity;
	  photonMultiplicity = iJet->photonMultiplicity;
	  electronMultiplicity = iJet->electronMultiplicity;
	  muonMultiplicity = iJet->muonMultiplicity;
	  HFHadronMultiplicity = iJet->HFHadronMultiplicity;
	  HFEMMultiplicity = iJet->HFEMMultiplicity;
	  chargedMultiplicity = iJet->chargedMultiplicity;
	  neutralMultiplicity = iJet->neutralMultiplicity;

	  break;
	}
	
      }
      
      if(!matched) {
	
	hasJetMatch = false;
	jet_dR = -1;
	jet_dPhi = -1;
	
	chargedHadronEnergy = -1;
	neutralHadronEnergy = -1;
	photonEnergy = -1;
	electronEnergy = -1;
	muonEnergy = -1;
	HFHadronEnergy = -1;
	HFEMEnergy = -1;
	chargedEmEnergy = -1;
	chargedMuEnergy = -1;
	neutralEmEnergy = -1;
	chargedHadronMultiplicity = -1;
	neutralHadronMultiplicity = -1;
	photonMultiplicity = -1;
	electronMultiplicity = -1;
	muonMultiplicity = -1;
	HFHadronMultiplicity = -1;
	HFEMMultiplicity = -1;
	chargedMultiplicity = -1;
	neutralMultiplicity = -1;
	
      }
     
      photonTree->Fill();
 
    }
  
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
  }
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
  cout << "fail MET filters          : " << nCnt[21][0] << endl;
  cout << "no passing candidates     : " << nCnt[28][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "events with no dijetpt    : " << nCnt[46][0] << endl;

  out->cd();
  out->Write();
  out->Close();
  */
}
