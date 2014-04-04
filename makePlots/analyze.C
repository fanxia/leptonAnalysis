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

void analyze(TString input, bool addMC, int channel, int intLumi_int, double metCut, int nPhotons_req, int nBtagReq, bool displayKStest, bool blinded) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  prep_signal(channels[channel], nPhotons_req);

  TFile * in = new TFile(input, "READ");
  TTree * ggTree = (TTree*)in->Get(channels[channel]+"_signalTree");
  TTree * qcdTree = (TTree*)in->Get(channels[channel]+"_eQCDTree");

  TFile * fSigA = new TFile("../acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get(channels[channel]+"_signalTree");

  TFile * fSigB = new TFile("../acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSigB->Get(channels[channel]+"_signalTree");

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  // Make the correlation plot for MET filters
  TH2D * metFilter = (TH2D*)in->Get("metFilter");
  if(channel == 0) {
    metFilter->GetXaxis()->SetLabelSize(0.035);
    metFilter->GetYaxis()->SetLabelSize(0.015);
    metFilter->GetZaxis()->SetLabelSize(0.02);

    metFilter->Draw("colz");
    metFilter->SetMarkerColor(kWhite);
    metFilter->Draw("text same");
    can->SetLogz(true);
    can->SaveAs("metFilter"+gifOrPdf);

    can->SetLogz(false);
  }

  PlotMaker * pMaker = new PlotMaker(intLumi_int, channels[channel], blinded);

  bool loadSuccess = true;
  
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic", 
					  245.8 * 0.457, 2.5, 3.4, 2.6,
					  channel, 0, kGray, "t#bar{t} inclusive");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep", 
					  245.8 * 0.438, 2.5, 3.4, 2.6,
					  channel, 0, kGray, "t#bar{t} inclusive");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep", 
					  245.8 * 0.105, 2.5, 3.4, 2.6,
					  channel, 0, kGray, "t#bar{t} inclusive");

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_WJetsToLNu.root", "WJetsToLNu", 36257.2, channel, 1, kOrange-3, "W + Jets");

  //loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dyJetsToLL.root", "dyJetsToLL", 3503.71, channel, 2, kYellow, "Z/#gamma* + Jets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy1JetsToLL.root", "dy1JetsToLL", 666.7, channel, 2, kYellow, "Z/#gamma* + Jets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy2JetsToLL.root", "dy2JetsToLL", 215.1, channel, 2, kYellow, "Z/#gamma* + Jets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy3JetsToLL.root", "dy3JetsToLL", 66.07, channel, 2, kYellow, "Z/#gamma* + Jets");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_dy4JetsToLL.root", "dy4JetsToLL", 27.38, channel, 2, kYellow, "Z/#gamma* + Jets");

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_s.root", "TBar_s", 1.76, channel, 3, kRed, "Single Top");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_t.root", "TBar_t", 30.7, channel, 3, kRed, "Single Top");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TBar_tW.root", "TBar_tW", 11.1, channel, 3, kRed, "Single Top");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_s.root", "T_s", 3.79, channel, 3, kRed, "Single Top");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_t.root", "T_t", 56.4, channel, 3, kRed, "Single Top");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_T_tW.root", "T_tW", 11.1, channel, 3, kRed, "Single Top");

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TTWJets.root", "TTWJets", 0.232, channel, 4, kAzure-2, "t#bar{t} + W/Z");
  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_TTZJets.root", "TTZJets", 0.2057, channel, 4, kAzure-2, "t#bar{t} + W/Z");

  loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttgjets.root", "ttgjets", 14.0, channel, 5, 8, "t#bar{t} + #gamma");

  //loadSuccess |= pMaker->LoadMCBackground("inputs/signal_contamination_ttGG.root", "ttGG", 0.146, channel, 6, kCyan+3, "t#bar{t} + #gamma#gamma");

  if(!loadSuccess) return;

  pMaker->SetTrees(ggTree, qcdTree,
		   sigaTree, sigbTree);

  pMaker->SetDisplayKStest(displayKStest);

  const int nMetBins = 17;
  Double_t xbins_met[nMetBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 300, 650};

  const int nKinematicBins = 28;
  Double_t xbins_kinematic[nKinematicBins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350, 
						400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};
  

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
  pMaker->BookHistogram("w_mT", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("m3", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("ele_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("muon_pt", nKinematicBins, xbins_kinematic);

  if(nPhotons_req >= 1) {
    pMaker->BookHistogram("leadPhotonEt", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);
    pMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
    pMaker->BookHistogram("leadSigmaIetaIeta", 40, 0., 0.02);
    pMaker->BookHistogram("leadChargedHadronIso", 75, 0, 15.0);
  }

  if(nPhotons_req >= 2) {
    pMaker->BookHistogram("trailPhotonEt", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("trailPhotonPhi", 63, -3.14159, 3.14159);
    pMaker->BookHistogram("trailPhotonEta", 40, -1.5, 1.5);
    pMaker->BookHistogram("trailSigmaIetaIeta", 40, 0, 0.02);
    pMaker->BookHistogram("trailChargedHadronIso", 75, 0, 15.0);
    pMaker->BookHistogram("diEMpT", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("diJetPt", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("photon_invmass", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("photon_dR", 50, 0., 5.);
    pMaker->BookHistogram("photon_dPhi", 35, 0., 3.14159);
  }

  pMaker->FillHistograms(metCut, nPhotons_req, nBtagReq);
  pMaker->SubtractMCFromQCD();

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("mcPlots_"+channels[channel]+".root", "RECREATE");

  pMaker->CreatePlot("Nphotons",
		     false,
		     "Number of #gamma's", "Number of Events",
		     0, 4,
		     2.e-3, 3.e6,
		     0., 1.3,
		     true, true, false,
		     out);

  pMaker->CreatePlot("pfMET",
		     true,
		     "#slash{E}_{T} (GeV)", "Number of Events",
		     xbins_met[0], xbins_met[nMetBins],
		     7.e-3, 2.5e4,
		     0., 2.8,
		     true, true, true,
		     out);

  pMaker->CreatePlot("Njets",
		     false,
		     "nJets", "Number of Events",
		     0, 14, 
		     2.e-2, 3.e6,
		     0., 1.8,
		     true, true, false,
		     out);
  
  pMaker->CreatePlot("Nbtags",
		     false,
		     "nBtags", "Number of Events",
		     0, 6, 
		     2.e-2, 3.e6,
		     0., 1.8,
		     true, true, false,
		     out);

   pMaker->CreatePlot("max_csv",
		      false,
		     "max csv", "Number of Events",
		     0.65, 1., 
		     2.e-2, 3.e6,
		     0., 1.8,
		     true, false, false,
		     out);

  pMaker->CreatePlot("submax_csv",
		     false,
		     "sub-max csv", "Number of Events",
		     0, 4, 
		     2.e-2, 3.e5,
		     0., 1.8,
		     true, false, false,
		     out);

  pMaker->CreatePlot("HT_jets",
		     true,
		     "HT (jets only) (GeV/c^{2})", "Number of Events",
		     0, 2000, 
		     2.e-4, 4.e3,
		     0., 3.5,
		     true, true, false,
		     out);

  pMaker->CreatePlot("hadronic_pt",
		     true,
		     "MHT (GeV/c)", "Number of Events",
		     0, 1500, 
		     2.e-4, 5.8e3,
		     0., 8.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("HT",
		     true,
		     "HT (GeV)", "Number of Events",
		     0, 2000, 
		     2.e-4, 2.8e3,
		     0., 3.8,
		     true, true, false,
		     out);

  pMaker->CreatePlot("jet1_pt",
		     true,
		     "Pt of leading jet", "Number of Events",
		     0, 1500, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     true, true, true,
		     out);

  pMaker->CreatePlot("jet2_pt",
		     true,
		     "Pt of sub-leading jet", "Number of Events",
		     0, 1200, 
		     2.e-4, 1.3e4,
		     0., 4.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("jet3_pt",
		     true,
		     "Pt of third-leading jet", "Number of Events",
		     0, 800, 
		     2.e-5, 8.e4,
		     0., 4.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("btag1_pt",
		     true,
		     "Pt of leading btag", "Number of Events",
		     0, 1400, 
		     2.e-4, 8.e3,
		     0., 4.8,
		     true, true, false,
		     out);
  
  pMaker->CreatePlot("w_mT",
		     true,
		     "Transverse Mass", "Number of Events",
		     0, 1000, 
		     2.e-4, 8.e3,
		     0., 3.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("m3",
		     true,
		     "M3 (GeV/c^{2})", "Number of Events",
		     0, 2000, 
		     2.e-4, 2.8e3,
		     0., 3.8,
		     true, true, false,
		     out);

  pMaker->CreatePlot("ele_pt",
		     true,
		     "Pt of electron", "Number of Events",
		     0, 1500, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     true, true, true,
		     out);

  pMaker->CreatePlot("muon_pt",
		     true,
		     "Pt of #mu", "Number of Events",
		     0, 1500, 
		     2.e-4, 8.e3,
		     0., 2.3,
		     true, true, true,
		     out);

  if(nPhotons_req >= 1) {
    pMaker->CreatePlot("leadPhotonEta",
		       false,
		       "#eta of leading #gamma", "Number of Events",
		       -1.5, 1.5, 
		       2.e-3, 3.e4,
		       0., 2.1,
		       false, false, false,
		       out);

    pMaker->CreatePlot("leadPhotonPhi",
		       false,
		       "#phi of leading #gamma", "Number of Events",
		       -3.2, 3.2, 
		       2.e-3, 3.e4,
		       0., 2.1,
		       false, false, false,
		       out);

    pMaker->CreatePlot("leadPhotonEt",
		       true,
		       "Et of leading #gamma", "Number of Events",
		       0, 1200, 
		       2.e-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);

    pMaker->CreatePlot("leadSigmaIetaIeta",
		       false,
		       "Lead #gamma #sigma_{i#etai#eta}", "Number of Events",
		       0, 0.015,
		       2.3-3, 5.e4,
		       0., 5.1,
		       true, true, true,
		       out);

    pMaker->CreatePlot("leadChargedHadronIso",
		       false,
		       "Lead #gamma ChHadIso", "Number of Events",
		       0, 15,
		       2.3e-4, 5.e3,
		       0., 1.8,
		       true, true, false,
		       out);

  }

  if(nPhotons_req >= 2) {
    pMaker->CreatePlot("photon_dR",
		       false,
		       "#DeltaR_{#gamma#gamma}", "Number of Events",
		       0.5, 5., 
		       2.e-2, 3.e5,
		       0., 2.1,
		       true, false, false,
		       out);

    pMaker->CreatePlot("photon_dPhi",
		       false,
		       "#Delta#phi_{#gamma#gamma}", "Number of Events",
		       0., 3.14159, 
		       2.e-2, 3.e5,
		       0., 2.1,
		       true, false, false,
		       out);

  pMaker->CreatePlot("trailPhotonEta",
		     false,
		     "#eta of trailing #gamma", "Number of Events",
		     -1.5, 1.5, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out);

  pMaker->CreatePlot("trailPhotonPhi",
		     false,
		     "#phi of trailing #gamma", "Number of Events",
		     -3.2, 3.2, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out);
  
  pMaker->CreatePlot("trailSigmaIetaIeta",
		     false,
		     "#sigma_{i#etai#eta} of Trail #gamma", "Number of Events",
		     0, 0.015,
		     2.3-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out);
  
  pMaker->CreatePlot("trailChargedHadronIso",
		     false,
		     "Ch. Hadron Iso. of Trail #gamma", "Number of Events",
		     0, 15,
		     2.3e-4, 5.e3,
		     0., 1.8,
		     true, true, false,
		     out);

  pMaker->CreatePlot("photon_invmass",
		     true,
		     "m_{#gamma#gamma} (GeV/c^{2})", "Number of Events",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("trailPhotonEt",
		     true,
		     "Et of trailing #gamma", "Number of Events",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out);

  pMaker->CreatePlot("diEMpT",
		     true,
		     "di-EM Pt", "Number of Events",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out);
  
  pMaker->CreatePlot("diJetPt",
		     true,
		     "di-Jet Pt", "Number of Events",
		     0, 1400, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out);

  }

  //pMaker->CreateTable();

  pMaker->PlotKolmogorovValues();
  
  delete pMaker;
  
  out->Close();

  in->Close();
  fSigA->Close();
  fSigB->Close();

}
