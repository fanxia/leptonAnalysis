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

void mvaTreeMaker(TString input, int channelNumber) {

  const int nChannels = 8;
  TString channels[nChannels] = {
    "nojet",
    "j", "b",
    "jj",
    "bj",
    "muJets", // gg+mu+bj + X (dilep veto)
    "eleJets",
    "hadronic" // gg+5j1b + X (lep veto)
  };

  TString channel = channels[channelNumber];

  TFile * in = new TFile(input, "READ");

  TTree * ffTree = (TTree*)in->Get("ff_"+channel+"_EvtTree");
  TTree * gfTree = (TTree*)in->Get("gf_"+channel+"_EvtTree");
  TTree * egTree = (TTree*)in->Get("eg_"+channel+"_EvtTree");
  TTree * ggTree = (TTree*)in->Get("gg_"+channel+"_EvtTree");

  TString emOrJet = "jet";
  TH2D * diempt_gg = (TH2D*)in->Get("di"+emOrJet+"pt_gg_"+channel); diempt_gg->Sumw2();
  TH2D * diempt_ff = (TH2D*)in->Get("di"+emOrJet+"pt_ff_"+channel); diempt_ff->Sumw2();
  TH2D * diempt_gf = (TH2D*)in->Get("di"+emOrJet+"pt_gf_"+channel); diempt_gf->Sumw2();
  
  // Make ratios
  TH1D * diempt_gg_0 = (TH1D*)diempt_gg->ProjectionX("gg_px0", 1, 1, "e");
  TH1D * diempt_gg_1 = (TH1D*)diempt_gg->ProjectionX("gg_px1", 2, 2, "e");
  TH1D * diempt_gg_2 = (TH1D*)diempt_gg->ProjectionX("gg_px2", 3, -1, "e");

  Float_t gg_total = diempt_gg->Integral();
  Float_t ff_total = diempt_ff->Integral();
  Float_t gf_total = diempt_gf->Integral();

  TH1D * ratio_ff_0 = GetWeights(diempt_gg_0, diempt_ff->ProjectionX("ff_px0", 1, 1, "e"), gg_total, ff_total);
  TH1D * ratio_ff_1 = GetWeights(diempt_gg_1, diempt_ff->ProjectionX("ff_px1", 2, 2, "e"), gg_total, ff_total);
  TH1D * ratio_ff_2 = GetWeights(diempt_gg_2, diempt_ff->ProjectionX("ff_px2", 3, -1, "e"), gg_total, ff_total);

  TH1D * ratio_gf_0 = GetWeights(diempt_gg_0, diempt_gf->ProjectionX("gf_px0", 1, 1, "e"), gg_total, gf_total);
  TH1D * ratio_gf_1 = GetWeights(diempt_gg_1, diempt_gf->ProjectionX("gf_px1", 2, 2, "e"), gg_total, gf_total);
  TH1D * ratio_gf_2 = GetWeights(diempt_gg_2, diempt_gf->ProjectionX("gf_px2", 3, -1, "e"), gg_total, gf_total);

  if(channel == "muJets" || channel == "eleJets") {

    // Doesn't really matter what histo goes in this function, just returns 1 +- 0

    ratio_ff_0 = GetFlatWeights(diempt_gg_0);
    ratio_ff_1 = GetFlatWeights(diempt_gg_0);
    ratio_ff_2 = GetFlatWeights(diempt_gg_0);

    ratio_gf_0 = GetFlatWeights(diempt_gg_0);
    ratio_gf_1 = GetFlatWeights(diempt_gg_0);
    ratio_gf_2 = GetFlatWeights(diempt_gg_0);
  }
    

  TH1D * trialWeights_ff;
  GetTrialWeights(ggTree, ffTree, channel, trialWeights_ff);

  TH1D * trialWeights_gf;
  GetTrialWeights(ggTree, gfTree, channel, trialWeights_gf);

  TFile * out = new TFile("mvaTree_"+channel+".root", "RECREATE");
  out->cd();

  TTree * ffNewTree = new TTree("ff_"+channel+"_EvtTree", "An event tree for MVA analysis");
  TTree * gfNewTree = new TTree("gf_"+channel+"_EvtTree", "An event tree for MVA analysis");
  TTree * egNewTree = new TTree("eg_"+channel+"_EvtTree", "An event tree for MVA analysis");
  TTree * ggNewTree = new TTree("gg_"+channel+"_EvtTree", "An event tree for MVA analysis");

  // Define here which variables you want to keep.
  vector<TString> floatNames;
  floatNames.push_back("diJetPt");
  floatNames.push_back("pfMET");
  floatNames.push_back("HT");
  floatNames.push_back("invmass");
  floatNames.push_back("leadPhotonEt");
  floatNames.push_back("leadPhotonEta");
  floatNames.push_back("leadPhotonPhi");
  floatNames.push_back("trailPhotonEt");
  floatNames.push_back("trailPhotonEta");
  floatNames.push_back("trailPhotonPhi");
  floatNames.push_back("photon_dR");
  floatNames.push_back("photon_dPhi");
  floatNames.push_back("jet1_pt");
  floatNames.push_back("jet2_pt");
  floatNames.push_back("diEMpT");
  floatNames.push_back("max_csv");
  floatNames.push_back("leadptOverInvmass");
  floatNames.push_back("trailptOverInvmass");
  floatNames.push_back("leadMatchedJetPt");
  floatNames.push_back("trailMatchedJetPt");
  floatNames.push_back("btag1_pt");
  floatNames.push_back("btag2_pt");
  floatNames.push_back("jet3_pt");
  floatNames.push_back("jet4_pt");
  floatNames.push_back("HT_jets");
  floatNames.push_back("hadronic_pt");
  floatNames.push_back("minDR_leadPhoton_jets");
  floatNames.push_back("minDR_trailPhoton_jets");
  floatNames.push_back("minDPhi_gMET");
  floatNames.push_back("minDPhi_jMET");
  floatNames.push_back("submax_csv");
  floatNames.push_back("w_mT");
  

  vector<TString> intNames;
  intNames.push_back("Njets");
  intNames.push_back("Nbtags");
  intNames.push_back("Nelectrons");
  intNames.push_back("Nmuons");

  vector<Float_t> floatVariablesFF, floatVariablesGF, floatVariablesEG, floatVariablesGG;
  floatVariablesFF.resize(floatNames.size());
  floatVariablesGF.resize(floatNames.size());
  floatVariablesEG.resize(floatNames.size());
  floatVariablesGG.resize(floatNames.size());

  vector<Int_t> intVariablesFF, intVariablesGF, intVariablesEG, intVariablesGG;
  intVariablesFF.resize(intNames.size());
  intVariablesGF.resize(intNames.size());
  intVariablesEG.resize(intNames.size());
  intVariablesGG.resize(intNames.size());

  for(unsigned int i = 0; i < floatNames.size(); i++) {
    ffTree->SetBranchAddress(floatNames[i], &(floatVariablesFF[i]));
    ffNewTree->Branch(floatNames[i], &(floatVariablesFF[i]), floatNames[i]+"/F");
    
    gfTree->SetBranchAddress(floatNames[i], &(floatVariablesGF[i]));
    gfNewTree->Branch(floatNames[i], &(floatVariablesGF[i]), floatNames[i]+"/F");

    egTree->SetBranchAddress(floatNames[i], &(floatVariablesEG[i]));
    egNewTree->Branch(floatNames[i], &(floatVariablesEG[i]), floatNames[i]+"/F");

    ggTree->SetBranchAddress(floatNames[i], &(floatVariablesGG[i]));
    ggNewTree->Branch(floatNames[i], &(floatVariablesGG[i]), floatNames[i]+"/F");
  }

  for(unsigned int i = 0; i < intNames.size(); i++) {
    ffTree->SetBranchAddress(intNames[i], &(intVariablesFF[i]));
    ffNewTree->Branch(intNames[i], &(intVariablesFF[i]), intNames[i]+"/I");

    gfTree->SetBranchAddress(intNames[i], &(intVariablesGF[i]));
    gfNewTree->Branch(intNames[i], &(intVariablesGF[i]), intNames[i]+"/I");

    egTree->SetBranchAddress(intNames[i], &(intVariablesEG[i]));
    egNewTree->Branch(intNames[i], &(intVariablesEG[i]), intNames[i]+"/I");
    
    ggTree->SetBranchAddress(intNames[i], &(intVariablesGG[i]));
    ggNewTree->Branch(intNames[i], &(intVariablesGG[i]), intNames[i]+"/I");
  }

  Float_t diemptWeight, diemptWeightErr;

  ffNewTree->Branch("weight", &diemptWeight);
  ffNewTree->Branch("weightError", &diemptWeightErr);

  for(int j = 0; j < ffTree->GetEntries(); j++) {
    ffTree->GetEntry(j);
    evaluateWeight(intVariablesFF[0], floatVariablesFF[0],
		   ratio_ff_0, ratio_ff_1, ratio_ff_2,
		   diemptWeight, diemptWeightErr);
    /*
      evaluateTrialWeight(floatVariablesFF[4],
                          trialWeights_ff,
                          diemptWeight, diemptWeightErr);
    */
    ffNewTree->Fill();
  }

  ffTree->ResetBranchAddresses();
  ffNewTree->ResetBranchAddresses();

  gfNewTree->Branch("weight", &diemptWeight);
  gfNewTree->Branch("weightError", &diemptWeightErr);

  for(int j = 0; j < gfTree->GetEntries(); j++) {
    gfTree->GetEntry(j);
    evaluateWeight(intVariablesGF[0], floatVariablesGF[0],
		   ratio_gf_0, ratio_gf_1, ratio_gf_2,
		   diemptWeight, diemptWeightErr);
    /*
      evaluateTrialWeight(floatVariablesGF[4],
                          trialWeights_gf,
                          diemptWeight_gf, diemptWeightErr_gf);
    */
    gfNewTree->Fill();
  }

  gfTree->ResetBranchAddresses();
  gfNewTree->ResetBranchAddresses();

  egNewTree->Branch("weight", &diemptWeight);
  egNewTree->Branch("weightError", &diemptWeightErr);

  for(int j = 0; j < egTree->GetEntries(); j++) {
    egTree->GetEntry(j);
    diemptWeight = 1.0;
    diemptWeightErr = 0.0;
    egNewTree->Fill();
  }

  egTree->ResetBranchAddresses();
  egNewTree->ResetBranchAddresses();

  ggNewTree->Branch("weight", &diemptWeight);
  ggNewTree->Branch("weightError", &diemptWeightErr);

  for(int j = 0; j < ggTree->GetEntries(); j++) {
    ggTree->GetEntry(j);
    diemptWeight = 1.0;
    diemptWeightErr = 0.0;
    ggNewTree->Fill();
  }

  ggTree->ResetBranchAddresses();
  ggNewTree->ResetBranchAddresses();

  out->Write();
  out->Close();
  in->Close();
}

void analyze(TString input, bool addMC, int channel, TString intLumi, int intLumi_int, bool useFF, bool useDifferenceSystematic, double metCut, bool displayKStest) {

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

  TFile * in = new TFile("mvaTree_"+channels[channel]+".root", "READ");

  // Start grabbing objects from input file
  TTree * ggTree = (TTree*)in->Get("gg_"+channels[channel]+"_EvtTree");
  TTree * gfTree = (TTree*)in->Get("gf_"+channels[channel]+"_EvtTree");
  TTree * ffTree = (TTree*)in->Get("ff_"+channels[channel]+"_EvtTree");
  TTree * egTree = (TTree*)in->Get("eg_"+channels[channel]+"_EvtTree");

  TFile * fSig460 = new TFile("../acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSig460->Get("gg_"+channels[channel]+"_EvtTree_mst_460_m1_175");

  TFile * fSig560 = new TFile("../acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSig560->Get("gg_"+channels[channel]+"_EvtTree_mst_560_m1_325");
  
  // pixel veto
  Float_t fakeRate = 0.019;
  Float_t fakeRate_err = 0.0001;
  Float_t fakeRate_sys = 0.002;

  // electron conversion veto
  //Float_t fakeRate = 0.08;
  //Float_t fakeRate_err = 0.02;

  // Calculate background scaling
  Float_t egScale = fakeRate/(1. - fakeRate);
  Float_t egScaleErr = fakeRate_err/(1. - fakeRate)/(1. - fakeRate);

  Float_t gfScale, gfScaleErr;
  bool gfWorks = calculateScaling(ggTree, egTree, gfTree, egScale, egScaleErr, gfScale, gfScaleErr);

  Float_t ffScale, ffScaleErr;
  bool ffWorks = calculateScaling(ggTree, egTree, ffTree, egScale, egScaleErr, ffScale, ffScaleErr);

  /*
  formatTable(met_gg,
	      met_eg_noNorm, ewk_normErr2, fakeRate, fakeRate_sys, egScale,
	      met_ff_Rew_noNorm, ff_norm2, ffScale,
	      met_gf_Rew_noNorm, gf_norm2, gfScale,
	      met_ee_Rew_noNorm, ee_norm2, eeScale,
	      met_ttgjets, useTTGJets,
	      channels[channel]);
  */

  if(channel == 0) {

    TFile * biggerInput = new TFile(input, "READ");

    TCanvas * can = new TCanvas("canvas_a", "Plot", 10, 10, 2000, 2000);
    
    // Make the correlation plot for MET filters
    TH2D * metFilter = (TH2D*)biggerInput->Get("metFilter");
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
    
    TH2D * DR_jet_gg = (TH2D*)biggerInput->Get("DR_jet_gg");
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
    
    biggerInput->Close();
  }

  PlotMaker * pMaker = new PlotMaker(intLumi_int, 
				     egScale, egScaleErr, 
				     ffScale, ffScaleErr, ffWorks,
				     gfScale, gfScaleErr, gfWorks,
				     useDifferenceSystematic, useFF,
				     channels[channel]);

  pMaker->SetTrees(ggTree, egTree,
		   ffTree, gfTree,
		   sigaTree, sigbTree);

  pMaker->SetDisplayKStest(displayKStest);

  TFile * limitOutput = new TFile("met_reweighted_"+channels[channel]+".root", "RECREATE");
  pMaker->SaveLimitOutput(limitOutput);
  limitOutput->Close();

  TFile * out = new TFile("plots_"+channels[channel]+".root", "RECREATE");

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
		     "min(#DeltaR_{#gamma_{lead}, jets}", "Number of Events",
		     0.5, 5., 
		     2.e-2, 3.e5,
		     0., 2.1,
		     false, false, false,
		     out, metCut);

  pMaker->CreatePlot("minDR_trailPhoton_jets", true,
		     50, 0., 5.,
		     "min(#DeltaR_{#gamma_{trail}, jets}", "Number of Events",
		     0.5, 5., 
		     2.e-2, 3.e5,
		     0., 2.1,
		     false, false, false,
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

  pMaker->PlotKolmogorovValues();

  delete pMaker;
    
  out->Close();

  in->Close();
  fSig460->Close();
  fSig560->Close();

}
