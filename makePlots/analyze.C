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

const TString ffColor = "kOrange+10";
const TString eeColor = "kBlue";
const TString egColor = "kGreen";

void analyze(TString input, bool addMC, int channel, int intLumi_int, double metCut, int nPhotons_req, bool displayKStest, bool blinded) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  const int nChannels = 8;
  TString channels[nChannels] = {"ele", "muon",
				 "ele_b", "muon_b",
				 "ele_jjj", "muon_jjj",
				 "ele_bjj", "muon_bjj"};
  
  prep_signal(channels[channel]);

  TFile * in = new TFile(input, "READ");

  TTree * ggTree = (TTree*)in->Get("gg_"+channels[channel]+"_EvtTree");

  TFile * fTTHadronic = new TFile("inputs/signal_contamination_ttJetsHadronic.root", "READ");
  TTree * ttHadronicTree = (TTree*)fTTHadronic->Get("gg_"+channels[channel]+"_EvtTree_ttJetsHadronic");
  TH1D * nGen_ttHadronic = (TH1D*)fTTHadronic->Get("nEvents_ttJetsHadronic");

  TFile * fTTSemiLep = new TFile("inputs/signal_contamination_ttJetsSemiLep.root", "READ");
  TTree * ttSemiLepTree = (TTree*)fTTSemiLep->Get("gg_"+channels[channel]+"_EvtTree_ttJetsSemiLep");
  TH1D * nGen_ttSemiLep = (TH1D*)fTTSemiLep->Get("nEvents_ttJetsSemiLep");

  TFile * fTTFullLep = new TFile("inputs/signal_contamination_ttJetsFullLep.root", "READ");
  TTree * ttFullLepTree = (TTree*)fTTFullLep->Get("gg_"+channels[channel]+"_EvtTree_ttJetsFullLep");
  TH1D * nGen_ttFullLep = (TH1D*)fTTFullLep->Get("nEvents_ttJetsFullLep");

  TFile * fWJets = new TFile("inputs/signal_contamination_WJetsToLNu.root", "READ");
  TTree * wjetsTree = (TTree*)fWJets->Get("gg_"+channels[channel]+"_EvtTree_WJetsToLNu");
  TH1D * nGen_wjets = (TH1D*)fWJets->Get("nEvents_WJetsToLNu");

  TFile * fDYJets = new TFile("inputs/signal_contamination_dyJetsToLL.root", "READ");
  TTree * dyjetsTree = (TTree*)fDYJets->Get("gg_"+channels[channel]+"_EvtTree_dyJetsToLL");
  TH1D * nGen_dyjets = (TH1D*)fDYJets->Get("nEvents_dyJetsToLL");

  TFile * fTTGJets = new TFile("inputs/signal_contamination_ttgjets.root", "READ");
  TTree * ttgjetsTree = (TTree*)fTTGJets->Get("gg_"+channels[channel]+"_EvtTree_ttgjets");
  TH1D * nGen_ttgjets = (TH1D*)fTTGJets->Get("nEvents_ttgjets");

  TFile * fTTGG = new TFile("inputs/signal_contamination_ttGG.root", "READ");
  TTree * ttggTree = (TTree*)fTTGG->Get("gg_"+channels[channel]+"_EvtTree_ttGG");
  TH1D * nGen_ttgg = (TH1D*)fTTGG->Get("nEvents_ttGG");

  TFile * fSigA = new TFile("../acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get("gg_"+channels[channel]+"_EvtTree_mst_460_m1_175");

  TFile * fSigB = new TFile("../acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSigB->Get("gg_"+channels[channel]+"_EvtTree_mst_560_m1_325");

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

  pMaker->SetTrees(ggTree,
		   ttHadronicTree, ttSemiLepTree, ttFullLepTree,
		   wjetsTree, dyjetsTree,
		   ttgjetsTree, ttggTree,
		   sigaTree, sigbTree);

  pMaker->GetNGen(nGen_ttHadronic, nGen_ttSemiLep, nGen_ttFullLep,
		  nGen_wjets, nGen_dyjets,
		  nGen_ttgjets, nGen_ttgg);

  pMaker->SetDisplayKStest(displayKStest);

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("mcPlots_"+channels[channel]+".root", "RECREATE");

  pMaker->CreatePlot("photon_dR",
		     50, 0., 5.,
		     "#DeltaR_{#gamma#gamma}", "Number of Events",
		     0.5, 5., 
		     2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("photon_dPhi",
		     35, 0., 3.14159,
		     "#Delta#phi_{#gamma#gamma}", "Number of Events",
		     0., 3.14159, 
		     2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("leadPhotonEta",
		     40, -1.5, 1.5,
		     "#eta of leading #gamma", "Number of Events",
		     -1.5, 1.5, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("trailPhotonEta",
		     40, -1.5, 1.5,
		     "#eta of trailing #gamma", "Number of Events",
		     -1.5, 1.5, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("leadPhotonPhi",
		     63, -3.14159, 3.14159,
		     "#phi of leading #gamma", "Number of Events",
		     -3.2, 3.2, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("trailPhotonPhi",
		     63, -3.14159, 3.14159,
		     "#phi of trailing #gamma", "Number of Events",
		     -3.2, 3.2, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("Njets",
		     20, 0., 20.,
		     "nJets", "Number of Events",
		     0, 9, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("Nbtags",
		     20, 0., 20.,
		     "nBtags", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("Nelectrons",
		     20, 0., 20.,
		     "nElectrons", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("Nmuons",
		     20, 0., 20.,
		     "nMuons", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("max_csv",
		     20, 0., 1.,
		     "max csv", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("submax_csv",
		     20, 0., 1.,
		     "sub-max csv", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut, nPhotons_req);

  const int nKinematicBins = 41;
  Double_t xbins_kinematic[nKinematicBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
						110, 120, 130, 140, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};

  pMaker->CreatePlot("HT_jets",
		     nKinematicBins, xbins_kinematic,
		     "HT (jets only) (GeV/c^{2})",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("hadronic_pt",
		     nKinematicBins, xbins_kinematic,
		     "Jet System Pt (GeV/c)",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("invmass",
		     nKinematicBins, xbins_kinematic,
		     "m_{#gamma#gamma} (GeV/c^{2})",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("HT",
		     nKinematicBins, xbins_kinematic,
		     "HT (GeV)",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("jet1_pt",
		     nKinematicBins, xbins_kinematic,
		     "Pt of leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("jet2_pt",
		     nKinematicBins, xbins_kinematic,
		     "Pt of sub-leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("jet3_pt",
		     nKinematicBins, xbins_kinematic,
		     "Pt of third-leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("jet4_pt",
		     nKinematicBins, xbins_kinematic,
		     "Pt of fourth-leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("btag1_pt",
		     nKinematicBins, xbins_kinematic,
		     "Pt of leading btag",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut, nPhotons_req);
  
  pMaker->CreatePlot("btag2_pt",
		     nKinematicBins, xbins_kinematic,
		     "Pt of sub-leading btag",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("leadPhotonEt",
		     nKinematicBins, xbins_kinematic,
		     "Et of leading #gamma",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("trailPhotonEt",
		     nKinematicBins, xbins_kinematic,
		     "Et of trailing #gamma",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreatePlot("diEMpT",
		     nKinematicBins, xbins_kinematic,
		     "di-EM Pt",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut, nPhotons_req);

  const int ndijetptbins = 31;
  Double_t dijetptbins[ndijetptbins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 200, 300, 400, 600, 1000, 1400};
  pMaker->CreatePlot("diJetPt",
		     ndijetptbins, dijetptbins,
		     "di-Jet Pt",
		     0, 1400, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut, nPhotons_req);

  const int nMetBins = 16;
  Double_t xbins_met[nMetBins+1] = {
    0,
    5,
    10,
    15,
    20,
    25,
    30,
    35,
    40,
    45,
    50,
    60,
    70,
    80,
    100,
    150,
    300};
  //650};

  pMaker->CreatePlot("pfMET",
		     nMetBins, xbins_met,
		     "#slash{E}_{T} (GeV)",
		     xbins_met[0], xbins_met[nMetBins],
		     7.e-4, 25000.,
		     0., 9.1,
		     true, true, true,
		     out, metCut, nPhotons_req);

  pMaker->CreateTable();

  pMaker->PlotKolmogorovValues();

  delete pMaker;
    
  out->Close();

  in->Close();
  fTTHadronic->Close();
  fTTSemiLep->Close();
  fTTFullLep->Close();
  fWJets->Close();
  fDYJets->Close();
  fTTGJets->Close();
  fTTGG->Close();
  fSigA->Close();
  fSigB->Close();

}
