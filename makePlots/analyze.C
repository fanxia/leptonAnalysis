#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>
#include <vector>
#include <stdio.h>
#include <stdarg.h>
#include <exception>

#include "analyze.h"

using namespace std;

void analyze(TString input, bool addMC, int channel, int intLumi_int, double metCut, int nPhotons_req, int nBtagReq, bool displayKStest, bool blinded, int photonMode) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * in = new TFile(input, "READ");

  TTree * ggTree = (TTree*)in->Get(channels[channel]+"_signalTree");
  if(photonMode == 1) ggTree = (TTree*)in->Get(channels[channel]+"_noSigmaIetaIetaTree");
  if(photonMode == 2) ggTree = (TTree*)in->Get(channels[channel]+"_noChHadIsoTree");

  TTree * qcdTree = (TTree*)in->Get(qcdChannels[channel]);
  if(photonMode == 1) qcdTree = (TTree*)in->Get(qcdChannels_noSigmaIetaIeta[channel]);
  if(photonMode == 2) qcdTree = (TTree*)in->Get(qcdChannels_noChHadIso[channel]);

  TFile * fSigA = new TFile("../acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get(channels[channel]+"_signalTree");
  TTree * sigaTree_JECup = (TTree*)fSigA->Get(channels[channel]+"_signalTree_JECup");
  TTree * sigaTree_JECdown = (TTree*)fSigA->Get(channels[channel]+"_signalTree_JECdown");

  TFile * fSigB = new TFile("../acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSigB->Get(channels[channel]+"_signalTree");
  TTree * sigbTree_JECup = (TTree*)fSigB->Get(channels[channel]+"_signalTree_JECup");
  TTree * sigbTree_JECdown = (TTree*)fSigB->Get(channels[channel]+"_signalTree_JECdown");

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  PlotMaker * pMaker = new PlotMaker(intLumi_int, channels[channel], blinded);
  pMaker->LoadLeptonSFs("../data/lepton_SF_8TeV_53x_baseline.root");
  pMaker->LoadPhotonSFs("../data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");

  bool loadSuccess = true;
  
  pMaker->SetPhotonMode(photonMode);

  Double_t ttbar_hadronic_xsec = 245.8 * 0.457;
  Double_t ttbar_semiLep_xsec  = 245.8 * 0.438;
  Double_t ttbar_fullLep_xsec  = 245.8 * 0.105;

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic", 
					  ttbar_hadronic_xsec, ttbar_hadronic_xsec * 0.025, ttbar_hadronic_xsec * 0.034, ttbar_hadronic_xsec * 0.026, ttbar_hadronic_xsec * 0.026,
					  true, true,
					  channel, 0, kGray, "t#bar{t} inclusive", "ttInclusive");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep", 
					  ttbar_semiLep_xsec, ttbar_semiLep_xsec * 0.025, ttbar_semiLep_xsec * 0.034, ttbar_semiLep_xsec * 0.026, ttbar_semiLep_xsec * 0.026,
					  true, true,
					  channel, 0, kGray, "t#bar{t} inclusive", "ttInclusive");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep", 
					  ttbar_fullLep_xsec, ttbar_fullLep_xsec * 0.025, ttbar_fullLep_xsec * 0.034, ttbar_fullLep_xsec * 0.026, ttbar_fullLep_xsec * 0.026,
					  true, true,
					  channel, 0, kGray, "t#bar{t} inclusive", "ttInclusive");

  /*
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WJetsToLNu.root", "WJetsToLNu", 
					  12234.4 * 3, 79.0, 39.7, 414.7, 414.7,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  */

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W1JetsToLNu.root", "W1JetsToLNu", 
					  12234.4 * 3 * 6662. / 37509., 79.0 * 6662. / 37509., 39.7 * 6662. / 37509., 414.7 * 6662. / 37509., 414.7 * 6662. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W2JetsToLNu.root", "W2JetsToLNu", 
					  12234.4 * 3 * 2159. / 37509., 79.0 * 2159. / 37509., 39.7 * 2159. / 37509., 414.7 * 2159. / 37509., 414.7 * 2159. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W3JetsToLNu.root", "W3JetsToLNu", 
					  12234.4 * 3 * 640. / 37509., 79.0 * 640. / 37509., 39.7 * 640. / 37509., 414.7 * 640. / 37509., 414.7 * 640. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W4JetsToLNu.root", "W4JetsToLNu", 
					  12234.4 * 3 * 264. / 37509., 79.0 * 264. / 37509., 39.7 * 264. / 37509., 414.7 * 264. / 37509., 414.7 * 264. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");


  /*
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dyJetsToLL.root", "dyJetsToLL", 
                                          1177.3 * 3, 5.9, 3.6, 38.8, 38.8,
                                          false, false,
                                          channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  */
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy1JetsToLL.root", "dy1JetsToLL", 
					  666.7 * 1177.3 * 3 / 3503.71, 5.9, 3.6, 38.8, 38.8,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy2JetsToLL.root", "dy2JetsToLL", 
					  215.1 * 1177.3 * 3 / 3503.71, 5.9 * 215.1 / 666.7, 3.6 * 215.1 / 666.7, 38.8 * 215.1 / 666.7, 38.8 * 215.1 / 666.7,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy3JetsToLL.root", "dy3JetsToLL", 
					  66.07 * 1177.3 * 3 / 3503.71, 5.9 * 66.07 / 666.7, 3.6 * 66.07 / 666.7, 38.8 * 66.07 / 666.7, 38.8 * 66.07 / 666.7,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy4JetsToLL.root", "dy4JetsToLL", 
					  27.38 * 1177.3 * 3 / 3503.71, 5.9 * 27.38 / 666.7, 3.6 * 27.38 / 666.7, 38.8 * 27.38 / 666.7, 38.8 * 27.38 / 666.7,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_s.root", "TBar_s", 
					  1.76, 0.01, 0.01, 0.08, 0.08,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_t.root", "TBar_t", 
					  30.7, 0.7, 0.7, 0.9, 1.1,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_tW.root", "TBar_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_s.root", "T_s", 
					  3.79, 0.07, 0.07, 0.13, 0.13,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_t.root", "T_t", 
					  56.4, 2.1, 0.3, 1.1, 1.1,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_tW.root", "T_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WW.root", "WW",
					  57.1097, 2.3, 2.3, 2.0, 2.0,
					  false, false,
					  channel, 4, kCyan, "Diboson", "diboson");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WZ.root", "WZ",
					  32.3161, 1.3, 1.3, 1.3, 1.3,
					  false, false,
					  channel, 4, kCyan, "Diboson", "diboson");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ZZ.root", "ZZ",
					  8.25561, 0.3, 0.3, 0.3, 0.3,
					  false, false,
					  channel, 4, kCyan, "Diboson", "diboson");
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TTWJets.root", "TTWJets", 
					  0.232, 0.067, 0.067, 0.03, 0.03,
					  false, false,
					  channel, 5, kAzure-2, "t#bar{t} + W/Z", "ttV");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TTZJets.root", "TTZJets", 
					  0.2057, 0., 0., 0.019, 0.024,
					  false, false,
					  channel, 5, kAzure-2, "t#bar{t} + W/Z", "ttV");

  // http://arxiv.org/abs/1102.1967
  //loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttgjets.root", "ttgjets", 
  //2.166, 2.166 * .25, 2.166 * .25, 2.166 * 0.076, 2.166 * 0.099,
  //false, true,
  //channel, 6, 8, "t#bar{t} + #gamma", "ttgamma");

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/WhizardMCTeeTeeGamma#2_to_5_All_ttbar_decay_channels
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttA_2to5.root", "ttA_2to5", 
					  .9081 * 2, .9081 * .5, .9081 * .5, .9081 * 2 * 0.076, .9081 * 2 * 0.099, 
					  false, true,
					  channel, 6, 8, "t#bar{t} + #gamma", "ttgamma");
  pMaker->SetUseWHIZARD(true);

  //loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttGG.root", "ttGG", 0.146, channel, 6, kCyan+3, "t#bar{t} + #gamma#gamma");

  if(!loadSuccess) return;

  pMaker->SetTrees(ggTree, qcdTree,
		   sigaTree, sigaTree_JECup, sigaTree_JECdown,
		   sigbTree, sigbTree_JECup, sigbTree_JECdown);

  pMaker->SetDisplayKStest(displayKStest);

  const int nMetBins_0g = 17;
  Double_t xbins_met_0g[nMetBins_0g+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 300, 650};
  const int nKinematicBins_0g = 28;
  Double_t xbins_kinematic_0g[nKinematicBins_0g+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};

  const int nMetBins_1g = 10;
  Double_t xbins_met_1g[nMetBins_1g+1] = {0, 10, 20, 30, 40, 50, 75, 100, 150, 300, 600};
  const int nKinematicBins_1g = 20;
  Double_t xbins_kinematic_1g[nKinematicBins_1g+1] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000};

  const int nMetBins_2g = 6;
  Double_t xbins_met_2g[nMetBins_2g+1] = {0, 25, 50, 75, 100, 150, 300};
  const int nKinematicBins_2g = 15;
  Double_t xbins_kinematic_2g[nKinematicBins_2g+1] = {0, 25, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000, 1500, 2000};

  // has to start with nphotons then met, then HT

  if(nPhotons_req == 0) {
    pMaker->BookHistogram("Nphotons", 4, 0., 4.);
    pMaker->BookHistogram("pfMET", nMetBins_0g, xbins_met_0g);
    pMaker->BookHistogram("HT", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("Njets", 20, 0., 20.);
    pMaker->BookHistogram("Nbtags", 20, 0., 20.);
    pMaker->BookHistogram("max_csv", 20, 0., 1.);
    pMaker->BookHistogram("submax_csv", 20, 0., 1.);
    pMaker->BookHistogram("HT_jets", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("hadronic_pt", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("jet1_pt", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("jet2_pt", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("jet3_pt", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("btag1_pt", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("w_mT", nKinematicBins_0g, xbins_kinematic_0g);    // 13
    pMaker->BookHistogram("m3", nKinematicBins_0g, xbins_kinematic_0g);
    pMaker->BookHistogram("ele_pt", nKinematicBins_0g, xbins_kinematic_0g);  // 15
    pMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    pMaker->BookHistogram("muon_pt", nKinematicBins_0g, xbins_kinematic_0g); // 17
    pMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);                  // 18

    /*
    pMaker->BookHistogram2D("Njets", "Nbtags", 15, 0., 15., 7, 0., 7.);
    pMaker->BookHistogram2D("HT", "pfMET", 20, 0., 1200., 20, 0., 350.);
    pMaker->BookHistogram2D("w_mT", "Njets", 60, 0., 600., 15, 0., 15.);
    pMaker->BookHistogram2D("w_mT", "Nbtags", 60, 0., 600., 7, 0., 7.);
    pMaker->BookHistogram2D("w_mT", "pfMET", 60, 0., 600., 20, 0., 350.);
    pMaker->BookHistogram2D("w_mT", "HT", 60, 0., 600., 20, 0., 1200.);
    */
  }
    
  if(nPhotons_req == 1) {
    pMaker->BookHistogram("Nphotons", 4, 0., 4.);
    pMaker->BookHistogram("pfMET", nMetBins_1g, xbins_met_1g);
    pMaker->BookHistogram("HT", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("Njets", 20, 0., 20.);
    pMaker->BookHistogram("Nbtags", 20, 0., 20.);
    pMaker->BookHistogram("max_csv", 20, 0., 1.);
    pMaker->BookHistogram("submax_csv", 20, 0., 1.);
    pMaker->BookHistogram("HT_jets", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("hadronic_pt", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("jet1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("jet2_pt", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("jet3_pt", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("btag1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("w_mT", nKinematicBins_1g, xbins_kinematic_1g);    // 13
    pMaker->BookHistogram("m3", nKinematicBins_1g, xbins_kinematic_1g);
    pMaker->BookHistogram("ele_pt", nKinematicBins_1g, xbins_kinematic_1g);  // 15
    pMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    pMaker->BookHistogram("muon_pt", nKinematicBins_1g, xbins_kinematic_1g); // 17
    pMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    pMaker->BookHistogram("leadPhotonEt", nKinematicBins_1g, xbins_kinematic_1g); // 19
    pMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
    pMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
    pMaker->BookHistogram("leadSigmaIetaIeta", 40, 0., 0.02);
    pMaker->BookHistogram("leadChargedHadronIso", 35, 0, 15.0);
    pMaker->BookHistogram("mLepGammaLead", nKinematicBins_1g, xbins_kinematic_1g);
    
    pMaker->BookHistogram2D("leadSigmaIetaIeta", "pfMET", 80, 0., 0.04, 20, 0., 350.);
    pMaker->BookHistogram2D("leadChargedHadronIso", "pfMET", 70, 0., 15., 20, 0., 350.);
  }
  
  if(nPhotons_req == 2) {
    pMaker->BookHistogram("Nphotons", 4, 0., 4.);
    pMaker->BookHistogram("pfMET", nMetBins_2g, xbins_met_2g);
    pMaker->BookHistogram("HT", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("Njets", 20, 0., 20.);
    pMaker->BookHistogram("Nbtags", 20, 0., 20.);
    pMaker->BookHistogram("max_csv", 20, 0., 1.);
    pMaker->BookHistogram("submax_csv", 20, 0., 1.);
    pMaker->BookHistogram("HT_jets", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("hadronic_pt", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("jet1_pt", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("jet2_pt", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("jet3_pt", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("btag1_pt", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("w_mT", nKinematicBins_2g, xbins_kinematic_2g);    // 13
    pMaker->BookHistogram("m3", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("ele_pt", nKinematicBins_2g, xbins_kinematic_2g);  // 15
    pMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    pMaker->BookHistogram("muon_pt", nKinematicBins_2g, xbins_kinematic_2g); // 17
    pMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    pMaker->BookHistogram("leadPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 19
    pMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
    pMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
    pMaker->BookHistogram("leadSigmaIetaIeta", 40, 0., 0.02);
    pMaker->BookHistogram("leadChargedHadronIso", 35, 0, 15.0);
    pMaker->BookHistogram("mLepGammaLead", nKinematicBins_2g, xbins_kinematic_2g);
    
    pMaker->BookHistogram2D("leadSigmaIetaIeta", "pfMET", 80, 0., 0.04, 20, 0., 350.);
    pMaker->BookHistogram2D("leadChargedHadronIso", "pfMET", 70, 0., 15., 20, 0., 350.);

    pMaker->BookHistogram("trailPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 24
    pMaker->BookHistogram("trailPhotonPhi", 63, -3.14159, 3.14159);          // 25
    pMaker->BookHistogram("trailPhotonEta", 40, -1.5, 1.5);                  // 26
    pMaker->BookHistogram("trailSigmaIetaIeta", 40, 0, 0.02);
    pMaker->BookHistogram("trailChargedHadronIso", 35, 0, 15.0);
    pMaker->BookHistogram("diEMpT", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("diJetPt", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("photon_invmass", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("photon_dR", 50, 0., 5.);
    pMaker->BookHistogram("photon_dPhi", 35, 0., 3.14159);
    pMaker->BookHistogram("mLepGammaTrail", nKinematicBins_2g, xbins_kinematic_2g);
    pMaker->BookHistogram("mLepGammaGamma", nKinematicBins_2g, xbins_kinematic_2g);
  }

  pMaker->FillHistograms(metCut, nPhotons_req, nBtagReq, channel);
  pMaker->SubtractMCFromQCD();

  pMaker->NormalizeQCD();

  pMaker->CreateFSRPlot(fSigA, fSigB);

  pMaker->CreateTable();
  pMaker->CreateAllDatacards(channel, nPhotons_req, nBtagReq);
  pMaker->SaveBackgroundOutput();

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("mcPlots_"+channels[channel]+".root", "RECREATE");

  bool needsQCD = (channel < 2);

  pMaker->Create2DPlots(needsQCD, true, out);

  pMaker->CreatePlot("Nphotons", false, needsQCD, "Number of #gamma's", "Number of Events",
		     0, 4, 2.e-2, 3.e6,
		     0.5, 1.5,
		     true, true, false,
		     out);

  // pfMET
  if(nPhotons_req == 0) pMaker->CreatePlot("pfMET", true, needsQCD, "#slash{E}_{T} (GeV)", "Number of Events",
					   0., 300., 7.e-3, 2.5e4,
					   0.7, 1.3,
					   true, true, true,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("pfMET", true, needsQCD, "#slash{E}_{T} (GeV)", "Number of Events",
						0., 300., 7.e-4, 2.5e2,
						0.5, 1.5,
						true, true, true,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("pfMET", true, needsQCD, "#slash{E}_{T} (GeV)", "Number of Events",
						0., 300., 7.e-6, 9.,
						0., 2.1,
						true, true, true,
						out);

  // Njets
  if(nPhotons_req == 0) pMaker->CreatePlot("Njets", false, needsQCD, "nJets", "Number of Events",
					   2, 14, 2.e-3, 9.e5,
					   0.6, 2.3,
					   true, true, false,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("Njets", false, needsQCD, "nJets", "Number of Events",
						2, 14, 2.e-3, 9.e4,
						0., 2.3,
						true, true, false,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("Njets", false, needsQCD, "nJets", "Number of Events",
						2, 14, 2.e-4, 9.e2,
						0., 2.3,
						true, true, false,
						out);
  
  // Nbtags
  if(nPhotons_req == 0) pMaker->CreatePlot("Nbtags", false, needsQCD, "nBtags", "Number of Events",
					   0, 8, 2.e-3, 3.e6, 
					   0.6, 1.6,
					   true, true, false,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("Nbtags", false, needsQCD, "nBtags", "Number of Events",
						0, 8, 2.e-3, 9.e4,
						0., 1.8,
						true, true, false,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("Nbtags", false, needsQCD, "nBtags", "Number of Events",
						0, 8, 2.e-3, 9.e2,
						0., 1.8,
						true, true, false,
						out);
  
  // max_csv
  if(nPhotons_req == 0) pMaker->CreatePlot("max_csv", false, needsQCD, "max csv", "Number of Events",
					   0.65, 1., 2.e-2, 3.e6,
					   0.7, 1.3,
					   true, false, false,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("max_csv", false, needsQCD, "max csv", "Number of Events",
						0.65, 1., 2.e-2, 5.e3,
						0., 1.8,
						true, false, false,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("max_csv", false, needsQCD, "max csv", "Number of Events",
						0.65, 1., 2.e-3, 8.e2,
						0., 2.1,
						true, false, false,
						out);
  
  // submax_csv
  if(nPhotons_req == 0) pMaker->CreatePlot("submax_csv", false, needsQCD, "sub-max csv", "Number of Events",
					   0, 1, 2.e-2, 3.e5,
					   0.6, 1.6,
					   true, false, false,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("submax_csv", false, needsQCD, "sub-max csv", "Number of Events",
						0, 1, 2.e-2, 3.e3,
						0.6, 1.8,
						true, false, false,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("submax_csv", false, needsQCD, "sub-max csv", "Number of Events",
						0, 1, 2.e-4, 3.e1,
						0., 1.8,
						true, false, false,
						out);
  // HT_jets
  if(nPhotons_req == 0) pMaker->CreatePlot("HT_jets", true, needsQCD, "HT (jets only) (GeV/c^{2})", "Number of Events",
					   0, 2000, 2.e-4, 4.e3,
					   0.5, 1.6,
					   true, true, false,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("HT_jets", true, needsQCD, "HT (jets only) (GeV/c^{2})", "Number of Events",
						0, 2000, 2.e-5, 4.e1,
						0., 2.5,
						true, true, false,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("HT_jets", true, needsQCD, "HT (jets only) (GeV/c^{2})", "Number of Events",
						0, 2000, 2.e-6, 4.e-1,
						0., 2.5,
						true, true, false,
						out);
  
  // MHT
  if(nPhotons_req == 0) pMaker->CreatePlot("hadronic_pt", true, needsQCD, "MHT (GeV/c)", "Number of Events",
					   0, 1500, 2.e-5, 8.e3,
					   0., 8.5,
					   true, true, true,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("hadronic_pt", true, needsQCD, "MHT (GeV/c)", "Number of Events",
						0, 1500, 2.e-5, 8.e1,
						0., 4.5,
						true, true, true,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("hadronic_pt", true, needsQCD, "MHT (GeV/c)", "Number of Events",
						0, 1500, 2.e-5, 8.e-1,
						0., 6.,
						true, true, true,
						out);
  
  // HT
  if(nPhotons_req == 0) pMaker->CreatePlot("HT", true, needsQCD, "HT (GeV)", "Number of Events",
					   0, 2000, 2.e-4, 2.8e3,
					   0.5, 1.5,
					   true, true, false,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("HT", true, needsQCD, "HT (GeV)", "Number of Events",
						0, 2000, 2.e-4, 4.8e1,
						0., 1.8,
						true, true, false,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("HT", true, needsQCD, "HT (GeV)", "Number of Events",
						0, 2000, 2.e-6, 2.8e-1,
						0., 2.5,
						true, true, false,
						out);
  
  // jet1_pt
  if(nPhotons_req == 0) pMaker->CreatePlot("jet1_pt", true, needsQCD, "Pt of leading jet", "Number of Events",
					   0, 1500, 2.e-4, 8.e3,
					   0., 2.3,
					   true, true, true,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("jet1_pt", true, needsQCD, "Pt of leading jet", "Number of Events",
						0, 1500, 2.e-5, 8.e1,
						0., 2.9,
						true, true, true,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("jet1_pt", true, needsQCD, "Pt of leading jet", "Number of Events",
						0, 1500, 2.e-6, 8.e-1,
						0., 2.3,
						true, true, true,
						out);

  if(nPhotons_req == 0) pMaker->CreatePlot("jet2_pt", true, needsQCD, "Pt of sub-leading jet", "Number of Events",
					   0, 1200, 2.e-4, 1.3e4,
					   0., 4.5,
					   true, true, true,
					   out);
  else if(nPhotons_req == 1) pMaker->CreatePlot("jet2_pt", true, needsQCD, "Pt of sub-leading jet", "Number of Events",
						0, 1200, 2.e-5, 1.3e2,
						0., 5.5,
						true, true, true,
						out);
  else if(nPhotons_req == 2) pMaker->CreatePlot("jet2_pt", true, needsQCD, "Pt of sub-leading jet", "Number of Events",
						0, 1200, 2.e-6, 1.3e0,
						0., 4.5,
						true, true, true,
						out);

  pMaker->CreatePlot("jet3_pt",
		     true, needsQCD,
		     "Pt of third-leading jet", "Number of Events",
		     0, 800, 
		     2.e-5, 8.e4,
		     0., 4.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("btag1_pt",
		     true, needsQCD,
		     "Pt of leading btag", "Number of Events",
		     0, 1400, 
		     2.e-4, 8.e3,
		     0., 4.8,
		     true, true, false,
		     out);
  
  pMaker->CreatePlot("w_mT",
		     true, needsQCD,
		     "Transverse Mass", "Number of Events",
		     0, 1000, 
		     2.e-4, 8.e3,
		     0., 3.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("m3",
		     true, needsQCD,
		     "M3 (GeV/c^{2})", "Number of Events",
		     0, 2000, 
		     2.e-4, 2.8e3,
		     0., 3.8,
		     true, true, false,
		     out);

  pMaker->CreatePlot("ele_pt",
		     true, needsQCD,
		     "Pt of electron", "Number of Events",
		     0, 1500, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     true, true, true,
		     out);

  pMaker->CreatePlot("ele_eta",
		     false, needsQCD,
		     "#Eta of electron", "Number of Events",
		     -2.4, 2.4, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     false, false, false,
		     out);

  pMaker->CreatePlot("muon_pt",
		     true, needsQCD,
		     "Pt of #mu", "Number of Events",
		     0, 1500, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     true, true, true,
		     out);

  pMaker->CreatePlot("muon_eta",
		     false, needsQCD,
		     "#Eta of #mu", "Number of Events",
		     -2.4, 2.4, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     false, false, false,
		     out);

  if(nPhotons_req >= 1) {
    pMaker->CreatePlot("leadPhotonEta",
		       false, needsQCD,
		       "#eta of leading #gamma", "Number of Events",
		       -1.5, 1.5, 
		       2.e-3, 3.e4,
		       0., 2.1,
		       false, false, false,
		       out);

    pMaker->CreatePlot("leadPhotonPhi",
		       false, needsQCD,
		       "#phi of leading #gamma", "Number of Events",
		       -3.2, 3.2, 
		       2.e-3, 3.e4,
		       0., 2.1,
		       false, false, false,
		       out);

    pMaker->CreatePlot("leadPhotonEt",
		       true, needsQCD,
		       "Et of leading #gamma", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);

    pMaker->CreatePlot("leadSigmaIetaIeta",
		       false, needsQCD,
		       "Lead #gamma #sigma_{i#etai#eta}", "Number of Events",
		       0, 0.02,
		       2.3e-3, 5.e4,
		       0., 5.1,
		       true, false, false,
		       out);

    pMaker->CreatePlot("leadChargedHadronIso",
		       false, needsQCD,
		       "Lead #gamma ChHadIso", "Number of Events",
		       0, 15,
		       2.3e-4, 2.5e4,
		       0., 1.8,
		       true, true, false,
		       out);

    pMaker->CreatePlot("mLepGammaLead",
		       true, needsQCD,
		       "m_{#ell#gamma_{lead}}", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);

  }

  if(nPhotons_req >= 2) {
    pMaker->CreatePlot("photon_dR",
		       false, needsQCD,
		       "#DeltaR_{#gamma#gamma}", "Number of Events",
		       0.5, 5., 
		       2.e-2, 3.e5,
		       0., 2.1,
		       true, false, false,
		       out);

    pMaker->CreatePlot("photon_dPhi",
		       false, needsQCD,
		       "#Delta#phi_{#gamma#gamma}", "Number of Events",
		       0., 3.14159, 
		       2.e-2, 3.e5,
		       0., 2.1,
		       true, false, false,
		       out);

    pMaker->CreatePlot("trailPhotonEta",
		       false, needsQCD,
		       "#eta of trailing #gamma", "Number of Events",
		       -1.5, 1.5, 
		       2.e-3, 3.e4,
		       0., 2.1,
		       false, false, false,
		       out);
    
    pMaker->CreatePlot("trailPhotonPhi",
		       false, needsQCD,
		       "#phi of trailing #gamma", "Number of Events",
		       -3.2, 3.2, 
		       2.e-3, 3.e4,
		       0., 2.1,
		       false, false, false,
		       out);
    
    pMaker->CreatePlot("trailSigmaIetaIeta",
		       false, needsQCD,
		       "#sigma_{i#etai#eta} of Trail #gamma", "Number of Events",
		       0, 0.02,
		       2.3e-3, 5.e4,
		       0., 5.1,
		       true, false, false,
		       out);
    
    pMaker->CreatePlot("trailChargedHadronIso",
		       false, needsQCD,
		       "Ch. Hadron Iso. of Trail #gamma", "Number of Events",
		       0, 15,
		       2.3e-4, 5.e3,
		       0., 1.8,
		       true, true, false,
		       out);
    
    pMaker->CreatePlot("photon_invmass",
		       true, needsQCD,
		       "m_{#gamma#gamma} (GeV/c^{2})", "Number of Events",
		       0, 2000, 
		       2.e-3, 3.e4,
		       0., 11.5,
		       true, true, true,
		       out);
    
    pMaker->CreatePlot("trailPhotonEt",
		       true, needsQCD,
		       "Et of trailing #gamma", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);
    
    pMaker->CreatePlot("diEMpT",
		       true, needsQCD,
		       "di-EM Pt", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);
    
    pMaker->CreatePlot("diJetPt",
		       true, needsQCD,
		       "di-Jet Pt", "Number of Events",
		       0, 1400, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);
   
    pMaker->CreatePlot("mLepGammaTrail",
		       true, needsQCD,
		       "m_{#ell#gamma_{trail}}", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);

    pMaker->CreatePlot("mLepGammaGamma",
		       true, needsQCD,
		       "m_{#ell#gamma#gamma}", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);

 
  }

  pMaker->PlotKolmogorovValues();
  
  delete pMaker;
  
  out->Close();

  in->Close();
  fSigA->Close();
  fSigB->Close();

}

void fitPhotons(TString input, bool addMC, int channel, int intLumi_int, double metCut, int nPhotons_req, int nBtagReq, bool displayKStest, bool blinded, int photonMode) {

  if(nPhotons_req != 1) return;

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * in = new TFile(input, "READ");
  TTree * ggTree = (TTree*)in->Get(channels[channel]+"_signalTree");
  if(photonMode == 1) ggTree = (TTree*)in->Get(channels[channel]+"_noSigmaIetaIetaTree");
  if(photonMode == 2) ggTree = (TTree*)in->Get(channels[channel]+"_noChHadIsoTree");

  TTree * qcdTree = (TTree*)in->Get(qcdChannels[channel]);
  if(photonMode == 1) qcdTree = (TTree*)in->Get(qcdChannels_noSigmaIetaIeta[channel]);
  if(photonMode == 2) qcdTree = (TTree*)in->Get(qcdChannels_noChHadIso[channel]);

  TFile * fSigA = new TFile("../acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get(channels[channel]+"_signalTree");
  TTree * sigaTree_JECup = (TTree*)fSigA->Get(channels[channel]+"_signalTree_JECup");
  TTree * sigaTree_JECdown = (TTree*)fSigA->Get(channels[channel]+"_signalTree_JECdown");

  TFile * fSigB = new TFile("../acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSigB->Get(channels[channel]+"_signalTree");
  TTree * sigbTree_JECup = (TTree*)fSigB->Get(channels[channel]+"_signalTree_JECup");
  TTree * sigbTree_JECdown = (TTree*)fSigB->Get(channels[channel]+"_signalTree_JECdown");

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  PlotMaker * pMaker = new PlotMaker(intLumi_int, channels[channel], blinded);
  pMaker->LoadLeptonSFs("../data/lepton_SF_8TeV_53x_baseline.root");
  pMaker->LoadPhotonSFs("../data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");

  bool loadSuccess = true;
  
  pMaker->SetPhotonMode(photonMode);

  Double_t ttbar_hadronic_xsec = 245.8 * 0.457;
  Double_t ttbar_semiLep_xsec  = 245.8 * 0.438;
  Double_t ttbar_fullLep_xsec  = 245.8 * 0.105;

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic", 
					  ttbar_hadronic_xsec, ttbar_hadronic_xsec * 0.025, ttbar_hadronic_xsec * 0.034, ttbar_hadronic_xsec * 0.026, ttbar_hadronic_xsec * 0.026,
					  true, true,
					  channel, 0, kGray, "t#bar{t} inclusive", "ttInclusive");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep", 
					  ttbar_semiLep_xsec, ttbar_semiLep_xsec * 0.025, ttbar_semiLep_xsec * 0.034, ttbar_semiLep_xsec * 0.026, ttbar_semiLep_xsec * 0.026,
					  true, true,
					  channel, 0, kGray, "t#bar{t} inclusive", "ttInclusive");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep", 
					  ttbar_fullLep_xsec, ttbar_fullLep_xsec * 0.025, ttbar_fullLep_xsec * 0.034, ttbar_fullLep_xsec * 0.026, ttbar_fullLep_xsec * 0.026,
					  true, true,
					  channel, 0, kGray, "t#bar{t} inclusive", "ttInclusive");

  /*
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WJetsToLNu.root", "WJetsToLNu", 
					  12234.4 * 3, 79.0, 39.7, 414.7, 414.7,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  */

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W1JetsToLNu.root", "W1JetsToLNu", 
					  12234.4 * 3 * 6662. / 37509., 79.0 * 6662. / 37509., 39.7 * 6662. / 37509., 414.7 * 6662. / 37509., 414.7 * 6662. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W2JetsToLNu.root", "W2JetsToLNu", 
					  12234.4 * 3 * 2159. / 37509., 79.0 * 2159. / 37509., 39.7 * 2159. / 37509., 414.7 * 2159. / 37509., 414.7 * 2159. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W3JetsToLNu.root", "W3JetsToLNu", 
					  12234.4 * 3 * 640. / 37509., 79.0 * 640. / 37509., 39.7 * 640. / 37509., 414.7 * 640. / 37509., 414.7 * 640. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_W4JetsToLNu.root", "W4JetsToLNu", 
					  12234.4 * 3 * 264. / 37509., 79.0 * 264. / 37509., 39.7 * 264. / 37509., 414.7 * 264. / 37509., 414.7 * 264. / 37509.,
					  false, false,
					  channel, 1, kOrange-3, "W + Jets", "vJets");


  /*
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dyJetsToLL.root", "dyJetsToLL", 
                                          1177.3 * 3, 5.9, 3.6, 38.8, 38.8,
                                          false, false,
                                          channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  */
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy1JetsToLL.root", "dy1JetsToLL", 
					  666.7 * 1177.3 * 3 / 3503.71, 5.9, 3.6, 38.8, 38.8,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy2JetsToLL.root", "dy2JetsToLL", 
					  215.1 * 1177.3 * 3 / 3503.71, 5.9 * 215.1 / 666.7, 3.6 * 215.1 / 666.7, 38.8 * 215.1 / 666.7, 38.8 * 215.1 / 666.7,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy3JetsToLL.root", "dy3JetsToLL", 
					  66.07 * 1177.3 * 3 / 3503.71, 5.9 * 66.07 / 666.7, 3.6 * 66.07 / 666.7, 38.8 * 66.07 / 666.7, 38.8 * 66.07 / 666.7,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy4JetsToLL.root", "dy4JetsToLL", 
					  27.38 * 1177.3 * 3 / 3503.71, 5.9 * 27.38 / 666.7, 3.6 * 27.38 / 666.7, 38.8 * 27.38 / 666.7, 38.8 * 27.38 / 666.7,
					  false, false,
					  channel, 2, kYellow, "Z/#gamma* + Jets", "vJets");
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_s.root", "TBar_s", 
					  1.76, 0.01, 0.01, 0.08, 0.08,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_t.root", "TBar_t", 
					  30.7, 0.7, 0.7, 0.9, 1.1,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_tW.root", "TBar_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_s.root", "T_s", 
					  3.79, 0.07, 0.07, 0.13, 0.13,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_t.root", "T_t", 
					  56.4, 2.1, 0.3, 1.1, 1.1,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_tW.root", "T_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7,
					  false, false,
					  channel, 3, kRed, "Single Top", "singleTop");

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WW.root", "WW",
					  57.1097, 2.3, 2.3, 2.0, 2.0,
					  false, false,
					  channel, 4, kCyan, "Diboson", "diboson");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WZ.root", "WZ",
					  32.3161, 1.3, 1.3, 1.3, 1.3,
					  false, false,
					  channel, 4, kCyan, "Diboson", "diboson");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ZZ.root", "ZZ",
					  8.25561, 0.3, 0.3, 0.3, 0.3,
					  false, false,
					  channel, 4, kCyan, "Diboson", "diboson");
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TTWJets.root", "TTWJets", 
					  0.232, 0.067, 0.067, 0.03, 0.03,
					  false, false,
					  channel, 5, kAzure-2, "t#bar{t} + W/Z", "ttV");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TTZJets.root", "TTZJets", 
					  0.2057, 0., 0., 0.019, 0.024,
					  false, false,
					  channel, 5, kAzure-2, "t#bar{t} + W/Z", "ttV");

  // http://arxiv.org/abs/1102.1967
  //loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttgjets.root", "ttgjets", 
  //2.166, 2.166 * .25, 2.166 * .25, 2.166 * 0.076, 2.166 * 0.099,
  //false, true,
  //channel, 6, 8, "t#bar{t} + #gamma", "ttgamma");

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/WhizardMCTeeTeeGamma#2_to_5_All_ttbar_decay_channels
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttA_2to5.root", "ttA_2to5", 
					  .9081 * 2, .9081 * .5, .9081 * .5, .9081 * 2 * 0.076, .9081 * 2 * 0.099, 
					  false, true,
					  channel, 6, 8, "t#bar{t} + #gamma", "ttgamma");
  pMaker->SetUseWHIZARD(true);

  //loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttGG.root", "ttGG", 0.146, channel, 6, kCyan+3, "t#bar{t} + #gamma#gamma");

  if(!loadSuccess) return;

  pMaker->SetTrees(ggTree, qcdTree,
		   sigaTree, sigaTree_JECup, sigaTree_JECdown,
		   sigbTree, sigbTree_JECup, sigbTree_JECdown);

  pMaker->SetDisplayKStest(displayKStest);

  const int nMetBins = 17;
  Double_t xbins_met[nMetBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 300, 650};

  const int nKinematicBins = 28;
  Double_t xbins_kinematic[nKinematicBins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};
  
  // has to start with nphotons then met, then HT
  pMaker->BookHistogram("Nphotons", 4, 0., 4.);
  pMaker->BookHistogram("pfMET", nMetBins, xbins_met);
  pMaker->BookHistogram("HT", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("Njets", 20, 0., 20.);
  pMaker->BookHistogram("Nbtags", 20, 0., 20.);
  pMaker->BookHistogram("max_csv", 20, 0., 1.);
  pMaker->BookHistogram("submax_csv", 20, 0., 1.);
  pMaker->BookHistogram("HT_jets", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("hadronic_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("jet1_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("jet2_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("jet3_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("btag1_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("w_mT", nKinematicBins, xbins_kinematic);    // 13
  pMaker->BookHistogram("m3", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("ele_pt", nKinematicBins, xbins_kinematic);  // 15
  pMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
  pMaker->BookHistogram("muon_pt", nKinematicBins, xbins_kinematic); // 17
  pMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);                  // 18

  pMaker->BookHistogram2D("Njets", "Nbtags", 15, 0., 15., 7, 0., 7.);
  pMaker->BookHistogram2D("HT", "pfMET", 20, 0., 1200., 20, 0., 350.);

  pMaker->BookHistogram2D("w_mT", "Njets", 60, 0., 600., 15, 0., 15.);
  pMaker->BookHistogram2D("w_mT", "Nbtags", 60, 0., 600., 7, 0., 7.);
  pMaker->BookHistogram2D("w_mT", "pfMET", 60, 0., 600., 20, 0., 350.);
  pMaker->BookHistogram2D("w_mT", "HT", 60, 0., 600., 20, 0., 1200.);  

  pMaker->BookHistogram("leadPhotonEt", nKinematicBins, xbins_kinematic); // 19
  pMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
  pMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
  if(channel == 3 && metCut < 0.) pMaker->BookHistogram("leadSigmaIetaIeta", 70, 0., 0.04);
  else pMaker->BookHistogram("leadSigmaIetaIeta", 80, 0., 0.04);
  pMaker->BookHistogram("leadChargedHadronIso", 70, 0, 15.0);
  pMaker->BookHistogram("mLepGammaLead", nKinematicBins, xbins_kinematic);

  pMaker->FillHistograms(metCut, nPhotons_req, nBtagReq, channel);
  pMaker->SubtractMCFromQCD();

  double ttjetsSF, ttjetsSFerror;
  double ttgammaSF, ttgammaSFerror;

  pMaker->NormalizeQCD();

  pMaker->FitSigmaIetaIeta(0.005, 0.025, nPhotons_req,
			   1., 0., 1., 0.,
			   1., 0., 1., 0.,
			   ttjetsSF, ttjetsSFerror, ttgammaSF, ttgammaSFerror);
  /*
  pMaker->FitChHadIso(0., 15.0, nPhotons_req,
		      1., 0., 1., 0.,
		      1., 0., 1., 0.,
		      ttjetsSF, ttjetsSFerror, ttgammaSF, ttgammaSFerror);
  */
  delete pMaker;
  
  in->Close();
  fSigA->Close();
  fSigB->Close();

}
