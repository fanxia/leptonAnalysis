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

void analyze(TString input, bool addMC, int channel, int intLumi_int, double metCut, int nPhotons_req, int nBtagReq, bool displayKStest, bool blinded) {

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
  
  prep_signal(channels[channel], nPhotons_req);

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

  TFile * fTBar_s = new TFile("inputs/signal_contamination_TBar_s.root", "READ");
  TTree * tbar_sTree = (TTree*)fTBar_s->Get("gg_"+channels[channel]+"_EvtTree_TBar_s");
  TH1D * nGen_tbar_s = (TH1D*)fTBar_s->Get("nEvents_TBar_s");

  TFile * fTBar_t = new TFile("inputs/signal_contamination_TBar_t.root", "READ");
  TTree * tbar_tTree = (TTree*)fTBar_t->Get("gg_"+channels[channel]+"_EvtTree_TBar_t");
  TH1D * nGen_tbar_t = (TH1D*)fTBar_t->Get("nEvents_TBar_t");

  TFile * fTBar_tW = new TFile("inputs/signal_contamination_TBar_tW.root", "READ");
  TTree * tbar_tWTree = (TTree*)fTBar_tW->Get("gg_"+channels[channel]+"_EvtTree_TBar_tW");
  TH1D * nGen_tbar_tW = (TH1D*)fTBar_tW->Get("nEvents_TBar_tW");

  TFile * fT_s = new TFile("inputs/signal_contamination_T_s.root", "READ");
  TTree * t_sTree = (TTree*)fT_s->Get("gg_"+channels[channel]+"_EvtTree_T_s");
  TH1D * nGen_t_s = (TH1D*)fT_s->Get("nEvents_T_s");

  TFile * fT_t = new TFile("inputs/signal_contamination_T_t.root", "READ");
  TTree * t_tTree = (TTree*)fT_t->Get("gg_"+channels[channel]+"_EvtTree_T_t");
  TH1D * nGen_t_t = (TH1D*)fT_t->Get("nEvents_T_t");

  TFile * fT_tW = new TFile("inputs/signal_contamination_T_tW.root", "READ");
  TTree * t_tWTree = (TTree*)fT_tW->Get("gg_"+channels[channel]+"_EvtTree_T_tW");
  TH1D * nGen_t_tW = (TH1D*)fT_tW->Get("nEvents_T_tW");

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
		   tbar_sTree, tbar_tTree, tbar_tWTree,
		   t_sTree, t_tTree, t_tWTree,
		   wjetsTree, dyjetsTree,
		   ttgjetsTree, ttggTree,
		   sigaTree, sigbTree);

  pMaker->SetNGen(nGen_ttHadronic, nGen_ttSemiLep, nGen_ttFullLep,
		  nGen_tbar_s, nGen_tbar_t, nGen_tbar_tW,
		  nGen_t_s, nGen_t_t, nGen_t_tW,
		  nGen_wjets, nGen_dyjets,
		  nGen_ttgjets, nGen_ttgg);

  pMaker->SetDisplayKStest(displayKStest);

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

  const int nKinematicBins = 41;
  Double_t xbins_kinematic[nKinematicBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
						110, 120, 130, 140, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};

  const int ndijetptbins = 31;
  Double_t dijetptbins[ndijetptbins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 200, 300, 400, 600, 1000, 1400};

  // has to start with nphotons then met
  pMaker->BookHistogram("Nphotons", 4, 0., 4.);
  pMaker->BookHistogram("pfMET", nMetBins, xbins_met);
  pMaker->BookHistogram("Njets", 20, 0., 20.);
  pMaker->BookHistogram("Nbtags", 20, 0., 20.);
  pMaker->BookHistogram("max_csv", 20, 0., 1.);
  pMaker->BookHistogram("submax_csv", 20, 0., 1.);
  pMaker->BookHistogram("HT_jets", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("hadronic_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("HT", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("jet1_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("jet2_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("jet3_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("btag1_pt", nKinematicBins, xbins_kinematic);
  pMaker->BookHistogram("w_mT", nKinematicBins, xbins_kinematic);

  if(nPhotons_req >= 1) {
    pMaker->BookHistogram("leadPhotonEt", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);
    pMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
  }

  if(nPhotons_req >= 2) {
    pMaker->BookHistogram("trailPhotonEt", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("diEMpT", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("diJetPt", ndijetptbins, dijetptbins);
    pMaker->BookHistogram("photon_invmass", nKinematicBins, xbins_kinematic);
    pMaker->BookHistogram("trailPhotonPhi", 63, -3.14159, 3.14159);
    pMaker->BookHistogram("trailPhotonEta", 40, -1.5, 1.5);
    pMaker->BookHistogram("photon_dR", 50, 0., 5.);
    pMaker->BookHistogram("photon_dPhi", 35, 0., 3.14159);
  }

  pMaker->FillHistograms(metCut, nPhotons_req, nBtagReq);

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("mcPlots_"+channels[channel]+".root", "RECREATE");

  pMaker->CreatePlot("pfMET",
		     true,
		     "#slash{E}_{T} (GeV)", "Number of Events",
		     xbins_met[0], xbins_met[nMetBins],
		     7.e-4, 25000.,
		     0., 9.1,
		     true, true, true,
		     out);

  pMaker->CreatePlot("Njets",
		     false,
		     "nJets", "Number of Events",
		     0, 9, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out);
  
  pMaker->CreatePlot("Nbtags",
		     false,
		     "nBtags", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out);

   pMaker->CreatePlot("max_csv",
		      false,
		     "max csv", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out);

  pMaker->CreatePlot("submax_csv",
		     false,
		     "sub-max csv", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out);

  pMaker->CreatePlot("HT_jets",
		     true,
		     "HT (jets only) (GeV/c^{2})", "Number of Events",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("hadronic_pt",
		     true,
		     "Jet System Pt (GeV/c)", "Number of Events",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("HT",
		     true,
		     "HT (GeV)", "Number of Events",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 5.1,
		     true, true, true,
		     out);

  pMaker->CreatePlot("jet1_pt",
		     true,
		     "Pt of leading jet", "Number of Events",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("jet2_pt",
		     true,
		     "Pt of sub-leading jet", "Number of Events",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("jet3_pt",
		     true,
		     "Pt of third-leading jet", "Number of Events",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out);

  pMaker->CreatePlot("btag1_pt",
		     true,
		     "Pt of leading btag", "Number of Events",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out);
  
  pMaker->CreatePlot("w_mT",
		     true,
		     "Transverse Mass", "Number of Events",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
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
  fTTHadronic->Close();
  fTTSemiLep->Close();
  fTTFullLep->Close();
  fTBar_s->Close();
  fTBar_t->Close();
  fTBar_tW->Close();
  fT_s->Close();
  fT_t->Close();
  fT_tW->Close();
  fWJets->Close();
  fDYJets->Close();
  fTTGJets->Close();
  fTTGG->Close();
  fSigA->Close();
  fSigB->Close();

}
