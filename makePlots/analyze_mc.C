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

#include "analyze_mc.h"

using namespace std;

const TString ffColor = "kOrange+10";
const TString eeColor = "kBlue";
const TString egColor = "kGreen";

void analyze(TString input, bool addMC, int channel, int intLumi_int, double metCut, bool useTTbar, bool useTTMBD, bool displayKStest, bool blinded) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  const int nChannels = 8;
  TString channels[nChannels] = {
    "nojet",
    "j", "b",
    "jj", "bj",
    "muJets", // gg+mu+bj + X (dilep veto)
    "eleJets",
    "hadronic" // gg+5j1b + X (lep veto)
  };

  prep_signal(channels[channel]);

  TFile * in = new TFile(input, "READ");

  // pixel veto
  Float_t fakeRate = 0.019;
  Float_t fakeRate_err = 0.0001;

  // electron conversion veto
  //Float_t fakeRate = 0.08;
  //Float_t fakeRate_err = 0.02;

  Float_t fakeRate_sys = 0.0006;

  // Scale the backgrounds
  Float_t egScale = fakeRate/(1. - fakeRate);
  Float_t egScaleErr = fakeRate_err/(1. - fakeRate)/(1. - fakeRate);

  TTree * ggTree = (TTree*)in->Get("gg_"+channels[channel]+"_EvtTree");
  TTree * egTree = (TTree*)in->Get("eg_"+channels[channel]+"_EvtTree");

  TFile * fTTGJets = new TFile("inputs/signal_contamination_ttgjets.root", "READ");
  TTree * ttgjetsTree = (TTree*)fTTGJets->Get("gg_"+channels[channel]+"_EvtTree_ttgjets");

  TFile * fTTMBD = new TFile("inputs/signal_contamination_TTJets.root");
  TTree * ttMBDTree = (TTree*)fTTMBD->Get("gg_"+channels[channel]+"_EvtTree_TTJets");

  TFile * fQCD30to40 = new TFile("inputs/signal_contamination_qcd30to40.root", "READ");
  TTree * qcd30to40Tree = (TTree*)fQCD30to40->Get("gg_"+channels[channel]+"_EvtTree_qcd30to40");

  TFile * fQCD40 = new TFile("inputs/signal_contamination_qcd40.root", "READ");
  TTree * qcd40Tree = (TTree*)fQCD40->Get("gg_"+channels[channel]+"_EvtTree_qcd40");

  TFile * fGJet20to40 = new TFile("inputs/signal_contamination_GJet20to40.root", "READ");
  TTree * gjet20to40Tree = (TTree*)fGJet20to40->Get("gg_"+channels[channel]+"_EvtTree_GJet20to40");

  TFile * fGJet40 = new TFile("inputs/signal_contamination_GJet40.root", "READ");
  TTree * gjet40Tree = (TTree*)fGJet40->Get("gg_"+channels[channel]+"_EvtTree_GJet40");

  TFile * fDiPhotonJets = new TFile("inputs/signal_contamination_DiPhotonJets.root", "READ");
  TTree * diphotonjetsTree = (TTree*)fDiPhotonJets->Get("gg_"+channels[channel]+"_EvtTree_DiPhotonJets");

  TFile * fDiphoBox10to25 = new TFile("inputs/signal_contamination_DiPhotonBox10To25.root");
  TTree * diphoBox10to25Tree = (TTree*)fDiphoBox10to25->Get("gg_"+channels[channel]+"_EvtTree_DiPhotonBox10To25");

  TFile * fDiphoBox25to250 = new TFile("inputs/signal_contamination_DiPhotonBox25To250.root");
  TTree * diphoBox25to250Tree = (TTree*)fDiphoBox25to250->Get("gg_"+channels[channel]+"_EvtTree_DiPhotonBox25To250");

  TFile * fDiphoBox250toInf = new TFile("inputs/signal_contamination_DiPhotonBox250ToInf.root");
  TTree * diphoBox250toInfTree = (TTree*)fDiphoBox250toInf->Get("gg_"+channels[channel]+"_EvtTree_DiPhotonBox250ToInf");

  TFile * fTTHadronic = new TFile("inputs/signal_contamination_ttJetsHadronic.root", "READ");
  TTree * ttHadronicTree = (TTree*)fTTHadronic->Get("gg_"+channels[channel]+"_EvtTree_ttJetsHadronic");
  
  TFile * fTTSemiLep = new TFile("inputs/signal_contamination_ttJetsSemiLep.root", "READ");
  TTree * ttSemiLepTree = (TTree*)fTTSemiLep->Get("gg_"+channels[channel]+"_EvtTree_ttJetsSemiLep");

  TFile * fTTFullLep = new TFile("inputs/signal_contamination_ttJetsFullLep.root", "READ");
  TTree * ttFullLepTree = (TTree*)fTTFullLep->Get("gg_"+channels[channel]+"_EvtTree_ttJetsFullLep");

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

  TH2D * DR_jet_gg = (TH2D*)in->Get("DR_jet_gg");
  if(channel == 0) {
    DR_jet_gg->Draw("colz");
    TLine * vertline = new TLine(0.5, 0, 0.5, 5);
    vertline->SetLineColor(kRed);
    vertline->SetLineWidth(3);
    vertline->Draw("same");
    TLine * horiline = new TLine(0, 0.5, 5, 0.5);
    horiline->SetLineColor(kRed);
    horiline->SetLineWidth(3);
    horiline->Draw("same");
    can->SetLogz(true);
    can->SaveAs("DR_jet_gg_"+channels[channel]+gifOrPdf);

    can->SetLogz(false);
  }

  PlotMaker * pMaker = new PlotMaker(intLumi_int, egScale, egScaleErr, channels[channel], blinded);
  pMaker->SetTrees(ggTree, egTree,
		   qcd30to40Tree, qcd40Tree,
		   gjet20to40Tree, gjet40Tree,
		   diphotonjetsTree,
		   diphoBox10to25Tree, diphoBox25to250Tree, diphoBox250toInfTree,
		   ttHadronicTree, ttSemiLepTree, ttFullLepTree,
		   ttgjetsTree,
		   ttMBDTree,
		   sigaTree, sigbTree);

  pMaker->SetUseTTbar(useTTbar);
  pMaker->SetUseTTMBD(useTTMBD);
  pMaker->SetDisplayKStest(displayKStest);

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("mcPlots_"+channels[channel]+".root", "RECREATE");

  pMaker->CreatePlot("photon_dR", true,
		     50, 0., 5.,
		     "#DeltaR_{#gamma#gamma}", "Number of Events",
		     0.5, 5., 
		     2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false,
		     out, metCut);

  pMaker->CreatePlot("minDR_leadPhoton_jets", true,
		     50, 0., 5.,
		     "min(#DeltaR_{#gamma_{lead}, jets})", "Number of Events",
		     0.5, 5., 
		     2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false,
		     out, metCut);

  pMaker->CreatePlot("minDR_trailPhoton_jets", true,
		     50, 0., 5.,
		     "min(#DeltaR_{#gamma_{trail}, jets})", "Number of Events",
		     0.5, 5., 
		     2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false,
		     out, metCut);

  pMaker->CreatePlot("photon_dPhi", true,
		     35, 0., 3.14159,
		     "#Delta#phi_{#gamma#gamma}", "Number of Events",
		     0., 3.14159, 
		     2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false,
		     out, metCut);

  pMaker->CreatePlot("minDPhi_gMET", true,
		     35, 0., 3.14159,
		     "min(#Delta#phi_{#gamma, MET})", "Number of Events",
		     0., 3.14159, 
		     2.e-2, 3.e5,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("minDPhi_jMET", true,
		     35, 0., 3.14159,
		     "min(#Delta#phi_{jets, MET})", "Number of Events",
		     0., 3.14159, 
		     2.e-2, 3.e5,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("leadPhotonEta", true,
		     40, -1.5, 1.5,
		     "#eta of leading #gamma", "Number of Events",
		     -1.5, 1.5, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("trailPhotonEta", true,
		     40, -1.5, 1.5,
		     "#eta of trailing #gamma", "Number of Events",
		     -1.5, 1.5, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("leadPhotonPhi", true,
		     63, -3.14159, 3.14159,
		     "#phi of leading #gamma", "Number of Events",
		     -3.2, 3.2, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("trailPhotonPhi", true,
		     63, -3.14159, 3.14159,
		     "#phi of trailing #gamma", "Number of Events",
		     -3.2, 3.2, 
		     2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("Njets", false,
		     20, 0., 20.,
		     "nJets", "Number of Events",
		     0, 9, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut);

  pMaker->CreatePlot("Nbtags", false,
		     20, 0., 20.,
		     "nBtags", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut);

  pMaker->CreatePlot("Nelectrons", false,
		     20, 0., 20.,
		     "nElectrons", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut);

  pMaker->CreatePlot("Nmuons", false,
		     20, 0., 20.,
		     "nMuons", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut);

  pMaker->CreatePlot("max_csv", true,
		     20, 0., 1.,
		     "max csv", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut);

  pMaker->CreatePlot("submax_csv", true,
		     20, 0., 1.,
		     "sub-max csv", "Number of Events",
		     0, 4, 
		     2.e-3, 3.e6,
		     0., 2.1,
		     true, true, false,
		     out, metCut);

  const int nKinematicBins = 41;
  Double_t xbins_kinematic[nKinematicBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
						110, 120, 130, 140, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};

  pMaker->CreatePlot("HT_jets", true,
		     nKinematicBins, xbins_kinematic,
		     "HT (jets only) (GeV/c^{2})",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("hadronic_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Jet System Pt (GeV/c)",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("invmass", true,
		     nKinematicBins, xbins_kinematic,
		     "m_{#gamma#gamma} (GeV/c^{2})",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("HT", true,
		     nKinematicBins, xbins_kinematic,
		     "HT (GeV)",
		     0, 2000, 
		     2.e-3, 3.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("jet1_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("jet2_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of sub-leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("jet3_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of third-leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("jet4_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of fourth-leading jet",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("btag1_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of leading btag",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut);
  
  pMaker->CreatePlot("btag2_pt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of sub-leading btag",
		     0, 1400, 
		     2.e-3, 8.e3,
		     0., 4.5,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("leadPhotonEt", true,
		     nKinematicBins, xbins_kinematic,
		     "Et of leading #gamma",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("trailPhotonEt", true,
		     nKinematicBins, xbins_kinematic,
		     "Et of trailing #gamma",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("diEMpT", true,
		     nKinematicBins, xbins_kinematic,
		     "di-EM Pt",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  const int ndijetptbins = 31;
  Double_t dijetptbins[ndijetptbins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 200, 300, 400, 600, 1000, 1400};
  pMaker->CreatePlot("diJetPt", true,
		     ndijetptbins, dijetptbins,
		     "di-Jet Pt",
		     0, 1400, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("leadMatchedJetPt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of jet matched to leading #gamma",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("trailMatchedJetPt", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt of jet matched to trailing #gamma",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("leadptOverInvmass", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt(lead #gamma) / m_{#gamma#gamma}",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreatePlot("trailptOverInvmass", true,
		     nKinematicBins, xbins_kinematic,
		     "Pt(trail #gamma) / m_{#gamma#gamma}",
		     0, 1200, 
		     2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true,
		     out, metCut);

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

  pMaker->CreatePlot("pfMET", true,
		     nMetBins, xbins_met,
		     "#slash{E}_{T} (GeV)",
		     xbins_met[0], xbins_met[nMetBins],
		     7.e-4, 25000.,
		     0., 9.1,
		     true, true, true,
		     out, metCut);

  pMaker->CreateTable();

  pMaker->PlotKolmogorovValues();

  delete pMaker;
    
  out->Close();

  in->Close();
  fTTGJets->Close();
  fQCD30to40->Close();
  fQCD40->Close();
  fGJet20to40->Close();
  fGJet40->Close();
  fTTHadronic->Close();
  fTTSemiLep->Close();
  fSigA->Close();
  fSigB->Close();

}
