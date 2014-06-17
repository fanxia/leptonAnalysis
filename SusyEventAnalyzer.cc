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

void SusyEventAnalyzer::PileupWeights(TString puFile, TString puFile_up, TString puFile_down) {

  TFile * in = new TFile(puFile, "READ");
  TFile * in_up = new TFile(puFile_up, "READ");
  TFile * in_down = new TFile(puFile_down, "READ");

  TH1F * _data = (TH1F*)in->Get("pileup");
  TH1F * _data_up = (TH1F*)in_up->Get("pileup");
  TH1F * _data_down = (TH1F*)in_down->Get("pileup");
  
  TString output_code_t = FormatName(scan);

  TH1F * data = (TH1F*)_data->Clone("pu_data"+output_code_t); data->Sumw2();
  TH1F * data_up = (TH1F*)_data_up->Clone("pu_data_up"+output_code_t); data_up->Sumw2();
  TH1F * data_down = (TH1F*)_data_down->Clone("pu_data_down"+output_code_t); data_down->Sumw2();

  TH1F * mc = new TH1F("pu_mc"+output_code_t, "pu_mc"+output_code_t, 70, 0, 70); mc->Sumw2();

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

  } // end event loop

  TH1D * data_nonorm = (TH1D*)data->Clone("pu_data_nonorm"+output_code_t);
  TH1D * data_up_nonorm = (TH1D*)data_up->Clone("pu_data_up_nonorm"+output_code_t);
  TH1D * data_down_nonorm = (TH1D*)data_down->Clone("pu_data_down_nonorm"+output_code_t);
  TH1D * mc_nonorm = (TH1D*)mc->Clone("pu_mc_nonorm"+output_code_t);

  Double_t intData = data->Integral();
  Double_t intData_up = data_up->Integral();
  Double_t intData_down = data_down->Integral();
  Double_t intMC = mc->Integral();

  data->Scale(1./intData);
  data_up->Scale(1./intData_up);
  data_down->Scale(1./intData_down);
  mc->Scale(1./intMC);

  TH1F * weights = (TH1F*)data->Clone("puWeights"+output_code_t);
  weights->Divide(mc);

  TH1F * weights_up = (TH1F*)data_up->Clone("puWeights_up"+output_code_t);
  weights_up->Divide(mc);

  TH1F * weights_down = (TH1F*)data_down->Clone("puWeights_down"+output_code_t);
  weights_down->Divide(mc);

  TFile * out = new TFile("pileupReweighting"+output_code_t+".root", "RECREATE");
  out->cd();

  data->Write();
  data_nonorm->Write();

  data_up->Write();
  data_up_nonorm->Write();

  data_down->Write();
  data_down_nonorm->Write();

  mc->Write();
  mc_nonorm->Write();

  weights->Write();
  weights_up->Write();
  weights_down->Write();

  out->Write();
  out->Close();

  in->Close();
  in_up->Close();
  in_down->Close();

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

  ScaleFactorInfo sf(btagger);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  vector<susy::PFJet*> pfJets, btags;
  vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
  vector<float> csvValues;
  vector<susy::Muon*> tightMuons, looseMuons;
  vector<susy::Electron*> tightEles, looseEles;
  vector<BtagInfo> tagInfos;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }

    nCnt[0][0]++; // events

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    float HT = 0.;

    tightMuons.clear();
    looseMuons.clear();
    tightEles.clear();
    looseEles.clear();
    pfJets.clear();
    btags.clear();
    pfJets_corrP4.clear();
    btags_corrP4.clear();
    csvValues.clear();
    tagInfos.clear();

    findMuons(event, tightMuons, looseMuons, HT, kSignal);
    if(tightMuons.size() > 1 || looseMuons.size() > 0) continue;

    findElectrons(event, tightMuons, looseMuons, tightEles, looseEles, HT, kSignal);
    if(tightEles.size() > 1 || looseEles.size() > 0) continue;

    if(tightMuons.size() + tightEles.size() != 1) continue;
    if(looseMuons.size() + looseEles.size() != 0) continue;

    bool passHLT = true;
    if(useTrigger) {
      if(tightEles.size() == 1) passHLT = PassTriggers(1);
      else if(tightMuons.size() == 1) passHLT = PassTriggers(2);
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

  TH1D * h_dR_gamma_ele = new TH1D("dR_gamma_ele", "dR between photons and electrons (N-1)", 100, 0, 10);
  TH1D * h_dR_gamma_muon = new TH1D("dR_gamma_muon", "dR between photons and muons (N-1)", 100, 0, 10);
  TH1D * h_dR_gamma_jet = new TH1D("dR_gamma_jet", "dR between photons and jets (N-1)", 100, 0, 10);
  TH1D * h_dR_gamma_photon = new TH1D("dR_gamma_photon", "dR between photons and other photons (N-1)", 100, 0, 10);

  /////////////////////////////////
  // Reweighting trees
  /////////////////////////////////

  const int nTreeVariables = 76;

  TString varNames[nTreeVariables] = {
    "pfMET", "pfMET_x", "pfMET_y", "pfMET_phi",
    "pfMET_sysShift", "pfMET_sysShift_phi",
    "pfMET_t1", "pfMET_t1p2", "pfMET_t01", "pfMET_t01p2", "pfNoPUMET", "pfMVAMET",
    "Njets", "Nbtags", "Nphotons", "Nmuons", "Nelectrons",
    "HT", "HT_jets", "hadronic_pt", 
    "w_mT", "w_mT_t1", "w_mT_t1p2", "w_mT_t01", "w_mT_t01p2", "w_mT_nopumet", "w_mT_mvamet",
    "m3",
    "ele_pt", "ele_phi", "ele_eta", "ele_mvaTrigV0", "ele_relIso",
    "muon_pt", "muon_phi", "muon_eta", "muon_relIso",
    "leadPhotonEt", "leadPhotonEta", "leadPhotonPhi", "leadChargedHadronIso", "leadSigmaIetaIeta", "lead_nPixelSeeds", "leadMVAregEnergy", "leadMVAregErr",
    "trailPhotonEt", "trailPhotonEta", "trailPhotonPhi", "trailChargedHadronIso", "trailSigmaIetaIeta", "trail_nPixelSeeds", "trailMVAregEnergy", "trailMVAregErr",
    "photon_invmass", "photon_dR", "photon_dPhi", "diEMpT", "diJetPt",
    "mLepGammaLead", "mLepGammaTrail", "mLepGammaGamma",
    "jet1_pt", "jet2_pt", "jet3_pt", "jet4_pt",
    "btag1_pt", "btag2_pt",
    "max_csv", "submax_csv", "min_csv",
    "nPV",
    "metFilterBit",
    "runNumber", "eventNumber", "lumiBlock", "jentry"};
    
  map<TString, float> treeMap;
  for(int i = 0; i < nTreeVariables; i++) treeMap[varNames[i]] = 0.;

  vector<TTree*> signalTrees, eQCDTrees, muQCDTrees, noSigmaIetaIetaTrees, noChHadIsoTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_signalTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    signalTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_eQCDTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    eQCDTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_muQCDTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    muQCDTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_noSigmaIetaIetaTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    noSigmaIetaIetaTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_noChHadIsoTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    noChHadIsoTrees.push_back(tree);
  }
    
  ScaleFactorInfo sf(btagger);

  bool quitAfterProcessing = false;

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  vector<susy::Muon*> tightMuons, looseMuons;
  vector<susy::Electron*> tightEles, looseEles;
  vector<susy::PFJet*> pfJets, btags;
  vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
  vector<float> csvValues;
  vector<susy::Photon*> photons;
  vector<BtagInfo> tagInfos;

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

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) {
      nCnt[22][0]++;
      continue;
    }

    for(int qcdMode = kSignal; qcdMode < kNumSearchModes; qcdMode++) {

      for(int photonMode = kSignalPhotons; photonMode < kNumPhotonModes; photonMode++) {

	if(qcdMode != kSignal && photonMode != kSignalPhotons) continue;

	float HT = 0.;
	
	tightMuons.clear();
	looseMuons.clear();
	tightEles.clear();
	looseEles.clear();
	pfJets.clear();
	btags.clear();
	pfJets_corrP4.clear();
	btags_corrP4.clear();
	csvValues.clear();
	photons.clear();
	tagInfos.clear();
	
	findMuons(event, tightMuons, looseMuons, HT, qcdMode);
	if(tightMuons.size() > 1 || looseMuons.size() > 0) {
	  nCnt[23][qcdMode]++;
	  continue;
	}
	
	findElectrons(event, tightMuons, looseMuons, tightEles, looseEles, HT, qcdMode);
	if(tightEles.size() > 1 || looseEles.size() > 0) {
	  nCnt[29][qcdMode]++;
	  continue;
	}
	
	if(tightMuons.size() + tightEles.size() != 1) {
	  nCnt[24][qcdMode]++;
	  continue;
	}
	if(looseMuons.size() + looseEles.size() != 0) {
	  nCnt[26][qcdMode]++;
	  continue;
	}
	
	bool passHLT = true;
	if(useTrigger) {
	  if(tightEles.size() == 1) passHLT = PassTriggers(1);
	  else if(tightMuons.size() == 1) passHLT = PassTriggers(2);
	}
	if(!passHLT) {
	  nCnt[25][qcdMode]++;
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
		    HT,
		    h_dR_gamma_ele,
		    h_dR_gamma_muon,
		    h_dR_gamma_jet,
		    h_dR_gamma_photon,
		    (photonMode != kNoSigmaIetaIeta), (photonMode != kNoChHadIso));
	
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
		      0, 0, 0, 0,
		      jentry);
	
	////////////////////
	
	for(unsigned int chan = 0; chan < nChannels; chan++) {
	  
	  if(pfJets.size() < nJetReq[chan]) continue;
	  if((nBtagInclusive[chan] && btags.size() < nBtagReq[chan]) || (!nBtagInclusive[chan] && btags.size() != nBtagReq[chan])) continue;
	  
	  if(tightEles.size() != nEleReq[chan]) continue;
	  if(tightMuons.size() != nMuonReq[chan]) continue;
	  
	  if(photonMode == kSignalPhotons) {
	    if(qcdMode == kSignal) {
	      nCnt[2][chan]++;
	      signalTrees[chan]->Fill();
	    }
	    else if(qcdMode == kElectronQCD) {
	      nCnt[3][chan]++;
	      eQCDTrees[chan]->Fill();
	    }
	    else if(qcdMode == kMuonQCD) {
	      nCnt[4][chan]++;
	      muQCDTrees[chan]->Fill();
	    }
	  }
	  
	  if(photonMode == kNoSigmaIetaIeta && qcdMode == kSignal) {
	    nCnt[5][chan]++;
	    noSigmaIetaIetaTrees[chan]->Fill();
	  }

	  if(photonMode == kNoChHadIso && qcdMode == kSignal) {
	    nCnt[6][chan]++;
	    noChHadIsoTrees[chan]->Fill();
	  }

	} // loop over jet/btag req channels
    
	///////////////////////////////////
    
      } // for photon modes

    } // for qcd modes

    if(quitAfterProcessing) break;
  } // for entries
  
  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "in_JSON              : " << nCnt[1][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "--------------- " << channels[i] << " Requirement ----------------" << endl;
    cout << "Signal          " << channels[i] << " events : " << nCnt[2][i] << endl;
    cout << "eQCD            " << channels[i] << " events : " << nCnt[3][i] << endl;
    cout << "muQCD           " << channels[i] << " events : " << nCnt[4][i] << endl;
    cout << "noSigmaIetaIeta " << channels[i] << " events : " << nCnt[5][i] << endl;
    cout << "noChHadIso      " << channels[i] << " events : " << nCnt[6][i] << endl;
  }
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
  cout << "fail MET filters         : " << nCnt[21][0] << endl;
  cout << "No primary vertex        : " << nCnt[22][0] << endl;
  cout << "Fail signal HLT          : " << nCnt[25][0] << endl;
  cout << "Fail eQCD HLT            : " << nCnt[25][1] << endl;
  cout << "Fail muQCD HLT           : " << nCnt[25][2] << endl;
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

  TH2D * h_ttA_phaseSpace = new TH2D("ttA_phaseSpace"+output_code_t, "ttA_phaseSpace"+output_code_t, 500, 0, 1000, 500, 0, 5);
  TH2D * h_ttbar_phaseSpace = new TH2D("ttbar_phaseSpace"+output_code_t, "ttbar_phaseSpace"+output_code_t, 500, 0, 1000, 500, 0, 5);

  TH1D * h_dR_gamma_ele = new TH1D("dR_gamma_ele", "dR between photons and electrons (N-1)", 100, 0, 10);
  TH1D * h_dR_gamma_muon = new TH1D("dR_gamma_muon", "dR between photons and muons (N-1)", 100, 0, 10);
  TH1D * h_dR_gamma_jet = new TH1D("dR_gamma_jet", "dR between photons and jets (N-1)", 100, 0, 10);
  TH1D * h_dR_gamma_photon = new TH1D("dR_gamma_photon", "dR between photons and other photons (N-1)", 100, 0, 10);

  const int nTreeVariables = 85;

  TString varNames[nTreeVariables] = {
    "pfMET", "pfMET_x", "pfMET_y", "pfMET_phi",
    "pfMET_sysShift", "pfMET_sysShift_phi",
    "pfMET_t1", "pfMET_t1p2", "pfMET_t01", "pfMET_t01p2", "pfNoPUMET", "pfMVAMET", "genMET",
    "Njets", "Nbtags", "Nphotons", "Nmuons", "Nelectrons",
    "HT", "HT_jets", "hadronic_pt", 
    "w_mT", "w_mT_t1", "w_mT_t1p2", "w_mT_t01", "w_mT_t01p2", "w_mT_nopumet", "w_mT_mvamet", "w_mT_genmet",
    "m3",
    "ele_pt", "ele_phi", "ele_eta", "ele_mvaTrigV0", "ele_relIso",
    "muon_pt", "muon_phi", "muon_eta", "muon_relIso",
    "leadPhotonEt", "leadPhotonEta", "leadPhotonPhi", "leadChargedHadronIso", "leadSigmaIetaIeta", "lead_nPixelSeeds", "leadMVAregEnergy", "leadMVAregErr",
    "trailPhotonEt", "trailPhotonEta", "trailPhotonPhi", "trailChargedHadronIso", "trailSigmaIetaIeta", "trail_nPixelSeeds", "trailMVAregEnergy", "trailMVAregErr",
    "photon_invmass", "photon_dR", "photon_dPhi", "diEMpT", "diJetPt",
    "mLepGammaLead", "mLepGammaTrail", "mLepGammaGamma",
    "jet1_pt", "jet2_pt", "jet3_pt", "jet4_pt",
    "btag1_pt", "btag2_pt",
    "max_csv", "submax_csv", "min_csv",
    "nPV",
    "pileupWeight", "pileupWeightErr", "pileupWeightUp", "pileupWeightDown",
    "btagWeight", "btagWeightUp", "btagWeightDown", "btagWeightErr",
    "metFilterBit",
    "ttbarDecayMode",
    "overlaps_ttA",
    "TopPtReweighting"};
    
  map<TString, float> treeMap;
  for(int i = 0; i < nTreeVariables; i++) treeMap[varNames[i]] = 0.;

  vector<TTree*> signalTrees, signalTrees_JECup, signalTrees_JECdown;
  vector<TTree*> eQCDTrees, muQCDTrees;
  vector<TTree*> noSigmaIetaIetaTrees, noChHadIsoTrees;

  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_signalTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    signalTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_signalTree_JECup", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    signalTrees_JECup.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_signalTree_JECdown", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    signalTrees_JECdown.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_eQCDTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    eQCDTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_muQCDTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    muQCDTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_noSigmaIetaIetaTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    noSigmaIetaIetaTrees.push_back(tree);
  }
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree(channels[i]+"_noChHadIsoTree", "An event tree for final analysis");
    for(int j = 0; j < nTreeVariables; j++) tree->Branch(varNames[j], &treeMap[varNames[j]], varNames[j]+"/F");
    noChHadIsoTrees.push_back(tree);
  }
  
  ScaleFactorInfo sf(btagger);
  TFile * btagEfficiency = new TFile("btagEfficiency"+output_code_t+".root", "READ");
  sf.SetTaggingEfficiencies((TH1F*)btagEfficiency->Get("lEff"+output_code_t), (TH1F*)btagEfficiency->Get("cEff"+output_code_t), (TH1F*)btagEfficiency->Get("bEff"+output_code_t));

  // get pileup weights
  TFile * puFile = new TFile("pileupReweighting"+output_code_t+".root", "READ");
  TH1F * puWeights = (TH1F*)puFile->Get("puWeights"+output_code_t);
  TH1F * puWeights_up = (TH1F*)puFile->Get("puWeights_up"+output_code_t);
  TH1F * puWeights_down = (TH1F*)puFile->Get("puWeights_down"+output_code_t);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;
  h_nEvents->Fill(0., (Double_t)nEntries);

  vector<susy::Muon*> tightMuons, looseMuons;
  vector<susy::Electron*> tightEles, looseEles;
  vector<susy::PFJet*> pfJets, btags;
  vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
  vector<float> csvValues;
  vector<susy::Photon*> photons;
  vector<BtagInfo> tagInfos;

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
    float eventWeightUp = 0.;
    float eventWeightDown = 0.;
    if(numTrueInt >= 0.) {
      int binNum = puWeights->GetXaxis()->FindBin(numTrueInt);
      eventWeight = puWeights->GetBinContent(binNum);
      eventWeightErr = puWeights->GetBinError(binNum);
      eventWeightUp = puWeights_up->GetBinContent(binNum);
      eventWeightDown = puWeights_down->GetBinContent(binNum);
    }

    if(!doPileupReweighting) {
      eventWeight = 1.;
      eventWeightErr = 0.;
      eventWeightUp = 1.;
      eventWeightDown = 1.;
    }

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;
    
    ttA_phaseSpace(event, h_ttA_phaseSpace);
    ttbar_phaseSpace(event, h_ttbar_phaseSpace);
    
    for(int qcdMode = kSignal; qcdMode < kNumSearchModes; qcdMode++) {
      for(int jetSyst = kCentral; jetSyst < kNumJetSytematics; jetSyst++) {
	for(int photonMode = kSignalPhotons; photonMode < kNumPhotonModes; photonMode++) {

	  if(qcdMode != kSignal && photonMode != kSignalPhotons) continue;
	  if(jetSyst != kCentral && photonMode != kSignalPhotons) continue;
	  if(qcdMode != kSignal && jetSyst != kCentral) continue;

	  if(jetSyst == kJERup || jetSyst == kJERdown) continue;
	  
	  float HT = 0.;
	  
	  tightMuons.clear();
	  looseMuons.clear();
	  tightEles.clear();
	  looseEles.clear();
	  pfJets.clear();
	  btags.clear();
	  pfJets_corrP4.clear();
	  btags_corrP4.clear();
	  csvValues.clear();
	  photons.clear();
	  tagInfos.clear();
	  
	  findMuons(event, tightMuons, looseMuons, HT, qcdMode);
	  if(tightMuons.size() > 1 || looseMuons.size() > 0) continue;
	  
	  findElectrons(event, tightMuons, looseMuons, tightEles, looseEles, HT, qcdMode);
	  if(tightEles.size() > 1 || looseEles.size() > 0) continue;
	  
	  if(tightMuons.size() + tightEles.size() != 1) continue;
	  if(looseMuons.size() + looseEles.size() != 0) continue;
	  
	  bool passHLT = true;
	  if(useTrigger) {
	    if(tightEles.size() == 1) passHLT = PassTriggers(1);
	    else if(tightMuons.size() == 1) passHLT = PassTriggers(2);
	  }
	  if(!passHLT) continue;
	  
	  float HT_jets = 0.;
	  TLorentzVector hadronicSystem(0., 0., 0., 0.);
	  
	  findJets_inMC(event, 
			tightMuons, looseMuons,
			tightEles, looseEles,
			pfJets, btags,
			sf,
			tagInfos, csvValues, 
			pfJets_corrP4, btags_corrP4, 
			HT_jets, hadronicSystem,
			jetSyst);
	  
	  findPhotons(event, 
		      photons,
		      pfJets_corrP4,
		      tightMuons, looseMuons,
		      tightEles, looseEles,
		      HT,
		      h_dR_gamma_ele,
		      h_dR_gamma_muon,
		      h_dR_gamma_jet,
		      h_dR_gamma_photon,
		      (photonMode != kNoSigmaIetaIeta), (photonMode != kNoChHadIso));
	  
	  float btagWeight[nChannels];
	  float btagWeightUp[nChannels];
	  float btagWeightDown[nChannels];
	  float btagWeightError[nChannels];
	  for(int chan = 0; chan < nChannels; chan++) {
	    BtagWeight * tagWeight = new BtagWeight(nBtagReq[chan]);
	    pair<float, float> weightResult = tagWeight->weight(tagInfos, btags.size(), 0., false, nBtagInclusive[chan]);
	    btagWeight[chan] = weightResult.first;
	    btagWeightError[chan] = weightResult.second;
	    
	    btagWeightUp[chan] = (tagWeight->weight(tagInfos, btags.size(), 1., true, nBtagInclusive[chan])).first;
	    btagWeightDown[chan] = (tagWeight->weight(tagInfos, btags.size(), -1., true, nBtagInclusive[chan])).first;
	    
	    delete tagWeight;
	  }
	  
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
			eventWeight, eventWeightErr, eventWeightUp, eventWeightDown,
			0);
	  
	  ////////////////////
	  
	  for(unsigned int chan = 0; chan < nChannels; chan++) {
	    
	    if(pfJets.size() < nJetReq[chan]) continue;
	    if((nBtagInclusive[chan] && btags.size() < nBtagReq[chan]) || (!nBtagInclusive[chan] && btags.size() != nBtagReq[chan])) continue;
	    
	    if(tightEles.size() != nEleReq[chan]) continue;
	    if(tightMuons.size() != nMuonReq[chan]) continue;
	    
	    treeMap["btagWeight"] = btagWeight[chan];
	    treeMap["btagWeightErr"] = btagWeightError[chan];
	    treeMap["btagWeightUp"] = btagWeightUp[chan];
	    treeMap["btagWeightDown"] = btagWeightDown[chan];

	    if(photonMode == kSignalPhotons) {
	      if(qcdMode == kSignal) {
		if(jetSyst == kCentral) {
		  nCnt[2][chan]++;
		  signalTrees[chan]->Fill();
		}
		else if(jetSyst == kJECup) signalTrees_JECup[chan]->Fill();
		else if(jetSyst == kJECdown) signalTrees_JECdown[chan]->Fill();
	      }
	      
	      else if(qcdMode == kElectronQCD) {
		nCnt[3][chan]++;
		eQCDTrees[chan]->Fill();
	      }
	      
	      else if(qcdMode == kMuonQCD) {
		nCnt[4][chan]++;
		muQCDTrees[chan]->Fill();
	      }
	      
	    }

	    if(photonMode == kNoSigmaIetaIeta && qcdMode == kSignal && jetSyst == kCentral) {
	      nCnt[5][chan]++;
	      noSigmaIetaIetaTrees[chan]->Fill();
	    }

	    if(photonMode == kNoChHadIso && qcdMode == kSignal && jetSyst == kCentral) {
	      nCnt[6][chan]++;
	      noChHadIsoTrees[chan]->Fill();
	    }

	  } // for channels

	} // for photon modes
	  
      } // for jet systematic modes
	
    } // for qcd modes

  } // for entries

  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "---------------- " << channels[i] << " Requirement ----------------" << endl;
    cout << "Signal " << channels[i] << " events : " << nCnt[2][i] << endl;
    if(nEleReq[i] == 1) cout << "eQCD   " << channels[i] << " events : " << nCnt[3][i] << endl;
    if(nMuonReq[i] == 1) cout << "muQCD   " << channels[i] << " events : " << nCnt[4][i] << endl;
  }
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
 
  puFile->Close();
  btagEfficiency->Close();

  out->Write();
  out->Close();

}

