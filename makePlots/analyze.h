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

#include "rootRoutines.h"

using namespace std;

const TString gifOrPdf = ".pdf";

const int nChannels = 4;
TString channels[nChannels] = {"ele_jjj", "ele_bjj",
			       "muon_jjj", "muon_bjj"};

TH1D * HistoFromTree(TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double metCut, int nPhotons_req) {
  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();
  FillHistoFromTree(h, tree, variable, metCut, nPhotons_req);
  return h;
}
  
TH1D * HistoFromTree(TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t* customBins, double metCut, int nPhotons_req) {
  TH1D * h = new TH1D(name, title, nBins, customBins);
  h->Sumw2();
  FillHistoFromTree(h, tree, variable, metCut, nPhotons_req);
  return h;
}

TH1D * SignalHistoFromTree(double scale, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double metCut, int nPhotons_req) {
  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();
  FillSignalHistoFromTree(h, tree, variable, metCut, nPhotons_req, scale);
  return h;
}

TH1D * SignalHistoFromTree(double scale, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t* customBins, double metCut, int nPhotons_req) {
  TH1D * h = new TH1D(name, title, nBins, customBins);
  h->Sumw2();
  FillSignalHistoFromTree(h, tree, variable, metCut, nPhotons_req, scale);
  return h;
}

void calculateROC(TH1D * sig_a, TH1D * sig_b, TH1D * bkg, TString req, TString title) {

  TH1D * srootb_a = (TH1D*)sig_a->Clone("roc_a"); srootb_a->Reset();
  TH1D * srootb_b = (TH1D*)sig_b->Clone("roc_b"); srootb_b->Reset();

  int nbins = sig_a->GetNbinsX();
  Double_t x_a[nbins], x_b[nbins], y_a[nbins], y_b[nbins];

  Double_t s, b, serr, berr;

  for(int i = 0; i < nbins; i++) {
    s = sig_a->IntegralAndError(i+1, -1, serr, "");
    b = bkg->IntegralAndError(i+1, -1, berr, "");
    
    if(b == 0.) continue;

    x_a[i] = s / sig_a->Integral();
    y_a[i] = b / bkg->Integral();
    srootb_a->SetBinContent(i+1, x_a[i] / sqrt(y_a[i]));

    s = sig_b->IntegralAndError(i+1, -1, serr, "");
    x_b[i] = s / sig_b->Integral();
    y_b[i] = b / bkg->Integral();
    srootb_b->SetBinContent(i+1, x_b[i] / sqrt(y_b[i]));

  }

  TGraph * roc_a = new TGraph(nbins, x_a, y_a);
  TGraph * roc_b = new TGraph(nbins, x_b, y_b);

  roc_a->SetLineColor(kMagenta);
  srootb_a->SetLineColor(kMagenta);
  roc_b->SetLineColor(kBlue);
  srootb_b->SetLineColor(kBlue);

  TCanvas * canv = new TCanvas("roc_can_"+title+"_"+req, "ROC Plot", 10, 10, 2000, 2000);

  TPad * padhi = new TPad("padhi", "padhi", 0, 0.5, 1, 1);
  TPad * padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.5);

  padhi->Draw();
  padlo->Draw();
  padhi->cd();

  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetGridx(true);
  padhi->SetGridy(true);

  TH2D * blank = new TH2D("blank_"+title+"_"+req, "blank;#epsilon_{S};#epsilon_{B}", 1, 0, 1, 1, 0, 1);
  blank->Draw();
  roc_a->Draw("same L");
  roc_b->Draw("same L");

  padlo->cd();
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  Int_t lastBin = 0;
  Int_t bestCutBin_a = 0;
  Int_t bestCutBin_b = 0;
  float bestDiscrim_a = 0;
  float bestDiscrim_b = 0;
  for(int i = 0; i < srootb_a->GetNbinsX(); i++) {
    if(srootb_a->GetBinContent(i+1) > 0 || srootb_b->GetBinContent(i+1) > 0) lastBin = i+1;
    if(srootb_a->GetBinContent(i+1) > bestDiscrim_a) {
      bestDiscrim_a = srootb_a->GetBinContent(i+1);
      bestCutBin_a = i+1;
    }
    if(srootb_b->GetBinContent(i+1) > bestDiscrim_b) {
      bestDiscrim_b = srootb_b->GetBinContent(i+1);
      bestCutBin_b = i+1;
    }
  }
  
  if(lastBin < srootb_a->GetNbinsX()) srootb_a->GetXaxis()->SetRangeUser(srootb_a->GetBinLowEdge(1), srootb_a->GetBinLowEdge(lastBin+1) * 1.1);
  if(srootb_b->GetMaximum() > srootb_a->GetMaximum()) srootb_a->GetYaxis()->SetRangeUser(0, 1.1 * srootb_b->GetMaximum());
  else srootb_a->GetYaxis()->SetRangeUser(0, 1.1 * srootb_a->GetMaximum());
  srootb_a->GetXaxis()->SetTitle("Lower cut on "+title);
  srootb_a->GetYaxis()->SetTitle("S / #sqrt{B}");
  srootb_a->Draw("hist");
  srootb_b->Draw("hist same");

  float lineMaxY = (srootb_b->GetMaximum() > srootb_a->GetMaximum()) ? 1.1 * srootb_b->GetMaximum() : 1.1 * srootb_a->GetMaximum();

  TLine * bestCutLine_a = new TLine(srootb_a->GetBinLowEdge(bestCutBin_a), 0, srootb_a->GetBinLowEdge(bestCutBin_a), lineMaxY);
  bestCutLine_a->SetLineColor(kMagenta);
  bestCutLine_a->SetLineWidth(2);
  bestCutLine_a->Draw("same");
  TLine * bestCutLine_b = new TLine(srootb_b->GetBinLowEdge(bestCutBin_b), 0, srootb_b->GetBinLowEdge(bestCutBin_b), lineMaxY);
  bestCutLine_b->SetLineColor(kBlue);
  bestCutLine_b->SetLineWidth(2);
  bestCutLine_b->Draw("same");

  canv->SaveAs("roc_"+title+"_"+req+".pdf");

}

class PlotMaker : public TObject {
  
  ClassDef(PlotMaker, 1);

 public:
  PlotMaker(Int_t lumi,
	    TString requirement,
	    bool blind);

  ~PlotMaker();
  
  bool LoadMCBackground(TString fileName, TString scanName,
			Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfError,
			int channel, int layer, int color, TString legendEntry);
  
  void SetTrees(TTree * gg, TTree * qcd,
		TTree * sig_a, TTree * sig_b);
  
  void SetDisplayKStest(bool v) { displayKStest = v; }

  void BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi);
  void BookHistogram(TString variable, Int_t nBins, Double_t* customBins);

  void FillHistograms(double metCut, int nPhotons_req, int nBtagReq);
  void SubtractMCFromQCD();

  void CreatePlot(TString variable,
		  bool divideByWidth,
		  TString xaxisTitle, TString yaxisTitle,
		  Float_t xmin, Float_t xmax,
		  Float_t ymin, Float_t ymax,
		  Float_t ratiomin, Float_t ratiomax,
		  bool drawSignal, bool drawLegend, bool drawPrelim,
		  TFile*& out);

  void DrawPlot(int variableNumber, TString variable,
		TString xaxisTitle, TString yaxisTitle,
		Float_t xmin, Float_t xmax,
		Float_t ymin, Float_t ymax,
		Float_t ratiomin, Float_t ratiomax,
		bool drawSignal, bool drawLegend, bool drawPrelim,
		TFile*& out);

  void CreateTable();

  void PlotKolmogorovValues();

 private:

  vector<TFile*> mcFiles;
  vector<TTree*> mcTrees;
  vector<TTree*> mcQCDTrees;
  vector<Double_t> mcNGen;
  vector<Double_t> crossSections;
  vector<Double_t> scaleErrUp;
  vector<Double_t> scaleErrDown;
  vector<Double_t> pdfErr;
  vector<int> mcLayerNumbers;
  vector<int> mcLayerColors;
  vector<TString> mcNames;
  vector<TString> legendNames;
  vector< vector<TH1D*> > mcHistograms;
  vector< vector<TH1D*> > mcQCDHistograms;

  TTree * ggTree;
  TTree * qcdTree;

  TTree * sigaTree;
  TTree * sigbTree;

  vector<TString> variables;

  vector<TH1D*> h_gg;
  vector<TH1D*> h_qcd;

  vector<TH1D*> h_siga;
  vector<TH1D*> h_sigb;

  vector<pair<TString, double> > KSscores;

  Int_t intLumi_int;
  TString intLumi;
  TString req;

  bool displayKStest;
  bool blinded;

};

PlotMaker::PlotMaker(Int_t lumi, TString requirement, bool blind) :
  intLumi_int(lumi),
  req(requirement),
  blinded(blind)
{
  char buffer[50];
  sprintf(buffer, "%.3f", (float)intLumi_int / 1000.);
  intLumi = buffer;

  displayKStest = false;

  KSscores.clear();

  variables.clear();

  h_gg.clear();
  h_qcd.clear();

  h_siga.clear();
  h_sigb.clear();

  mcFiles.clear();
  mcTrees.clear();
  mcQCDTrees.clear();
  mcNGen.clear();
  crossSections.clear();
  scaleErrUp.clear();
  scaleErrDown.clear();
  pdfErr.clear();
  mcLayerNumbers.clear();
  mcLayerColors.clear();
  mcNames.clear();
  legendNames.clear();
  mcHistograms.clear();
  mcQCDHistograms.clear();

}

PlotMaker::~PlotMaker() { 

    KSscores.clear();
    variables.clear();

    delete ggTree;
    delete qcdTree;

    mcHistograms.clear();
    mcQCDHistograms.clear();
    mcTrees.clear();
    mcQCDTrees.clear();
    mcFiles.clear();
    mcNGen.clear();
    crossSections.clear();
    scaleErrUp.clear();
    scaleErrDown.clear();
    pdfErr.clear();
    mcLayerNumbers.clear();
    mcLayerColors.clear();
    mcNames.clear();
    legendNames.clear();

    delete sigaTree;
    delete sigbTree;
    
    h_gg.clear();
    h_qcd.clear();

    h_siga.clear();
    h_sigb.clear();

}

void PlotMaker::SetTrees(TTree * gg, TTree * qcd,
			 TTree * sig_a, TTree * sig_b) {

  ggTree = gg;
  qcdTree = qcd;
  
  sigaTree = sig_a;
  sigbTree = sig_b;

}

bool PlotMaker::LoadMCBackground(TString fileName, TString scanName,
				 Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfError,
				 int channel, int layer, int color, TString legendEntry) {

  mcFiles.push_back(new TFile(fileName, "READ"));
  if(!mcFiles.back()) {
    cout << "Could not load TFile " << fileName << endl;
    return false;
  }

  mcTrees.push_back((TTree*)mcFiles.back()->Get(channels[channel]+"_signalTree"));
  if(!mcTrees.back()) {
    cout << "Could not load TTree " << channels[channel] << "_signalTree from TFile " << fileName << endl;
    return false;
  }

  mcQCDTrees.push_back((TTree*)mcFiles.back()->Get(channels[channel]+"_eQCDTree"));
  if(!mcQCDTrees.back()) {
    cout << "Could not load TTree " << channels[channel] << "_eQCDTree from TFile " << fileName << endl;
    return false;
  }

  TH1D * h_nGen = (TH1D*)mcFiles.back()->Get("nEvents_"+scanName);
  if(!h_nGen) {
    cout << "Could not load histogram " << "nEvents_" << scanName << " from TFile " << fileName << endl;
    return false;
  }
  mcNGen.push_back(h_nGen->Integral());

  crossSections.push_back(xsec);
  scaleErrUp.push_back(scaleErrorUp);
  scaleErrDown.push_back(scaleErrorDown);
  pdfErr.push_back(pdfError);
  mcLayerNumbers.push_back(layer);
  mcLayerColors.push_back(color);
  mcNames.push_back(scanName);
  legendNames.push_back(legendEntry);

  mcHistograms.resize(mcHistograms.size() + 1);
  mcQCDHistograms.resize(mcQCDHistograms.size() + 1);
  
  return true;
}

void PlotMaker::BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi) {
  
  variables.push_back(variable);

  TH1D * gg = new TH1D(variable+"_gg_"+req, variable, nBins, xlo, xhi);
  gg->Sumw2();
  h_gg.push_back(gg);

  TH1D * qcd = new TH1D(variable+"_qcd_"+req, variable, nBins, xlo, xhi);
  qcd->Sumw2();
  h_qcd.push_back(qcd);
  
  TH1D * h_bkg;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_"+mcNames[i]+"_"+req, variable, nBins, xlo, xhi);
    h_bkg->Sumw2();
    mcHistograms[i].push_back(h_bkg);
  }

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+req, variable, nBins, xlo, xhi);
    h_bkg->Sumw2();
    mcQCDHistograms[i].push_back(h_bkg);
  }

  TH1D * sig_a = new TH1D(variable+"_a_"+req, variable, nBins, xlo, xhi);
  sig_a->Sumw2();
  h_siga.push_back(sig_a);
  
  TH1D * sig_b = new TH1D(variable+"_b_"+req, variable, nBins, xlo, xhi);
  sig_b->Sumw2();
  h_sigb.push_back(sig_b);
  
}

void PlotMaker::BookHistogram(TString variable, Int_t nBins, Double_t* customBins) {

  variables.push_back(variable);

  TH1D * gg = new TH1D(variable+"_gg_"+req, variable, nBins, customBins);
  gg->Sumw2();
  h_gg.push_back(gg);
  
  TH1D * qcd = new TH1D(variable+"_qcd_"+req, variable, nBins, customBins);
  qcd->Sumw2();
  h_qcd.push_back(qcd);

  TH1D * h_bkg;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_"+mcNames[i]+"_"+req, variable, nBins, customBins);
    h_bkg->Sumw2();
    mcHistograms[i].push_back(h_bkg);
  }
  
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+req, variable, nBins, customBins);
    h_bkg->Sumw2();
    mcQCDHistograms[i].push_back(h_bkg);
  }

  TH1D * sig_a = new TH1D(variable+"_a_"+req, variable, nBins, customBins);
  sig_a->Sumw2();
  h_siga.push_back(sig_a);
  
  TH1D * sig_b = new TH1D(variable+"_b_"+req, variable, nBins, customBins);
  sig_b->Sumw2();
  h_sigb.push_back(sig_b);
  
}

// expects BookHistogram on nphotons, then met, then others
void PlotMaker::FillHistograms(double metCut, int nPhotons_req, int nBtagReq) {

  vector<Float_t> vars;
  vars.resize(variables.size());

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr, btagWeightUp, btagWeightDown;

  for(unsigned int i = 0; i < variables.size(); i++) {

    ggTree->SetBranchAddress(variables[i], &(vars[i]));
    qcdTree->SetBranchAddress(variables[i], &(vars[i]));

    for(unsigned int j = 0; j < mcTrees.size(); j++) {
      mcTrees[j]->SetBranchAddress(variables[i], &(vars[i]));
      mcQCDTrees[j]->SetBranchAddress(variables[i], &(vars[i]));
    }

    sigaTree->SetBranchAddress(variables[i], &(vars[i]));
    sigbTree->SetBranchAddress(variables[i], &(vars[i]));

  }

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    mcTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);
    
    mcQCDTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcQCDTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcQCDTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);
  }

  sigaTree->SetBranchAddress("pileupWeight", &puWeight);
  sigaTree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigaTree->SetBranchAddress("btagWeight", &btagWeight);
  sigaTree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigaTree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigaTree->SetBranchAddress("btagWeightDown", &btagWeightDown);

  sigbTree->SetBranchAddress("pileupWeight", &puWeight);
  sigbTree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigbTree->SetBranchAddress("btagWeight", &btagWeight);
  sigbTree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigbTree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigbTree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  
  for(int i = 0; i < ggTree->GetEntries(); i++) {
    ggTree->GetEntry(i);

    if(metCut > 0. && vars[1] >= metCut) continue;

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      if(blinded && vars[0] == 2) continue;
      if(blinded && vars[0] == 1 && variables[j] == "pfMET" && vars[1] > 50.) continue;
      if(blinded && vars[0] == 1 && variables[j] == "HT" && vars[2] > 400.) continue;

      h_gg[j]->Fill(vars[j]);
    }

  }

  for(int i = 0; i < qcdTree->GetEntries(); i++) {
    qcdTree->GetEntry(i);

    if(metCut > 0. && vars[1] >= metCut) continue;

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      h_qcd[j]->Fill(vars[j]);
    }

  }

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    
    for(int j = 0; j < mcTrees[i]->GetEntries(); j++) {
      mcTrees[i]->GetEntry(j);

      if(nBtagReq == 0) {
	btagWeight = 1.;
	btagWeightErr = 0.;
	btagWeightUp = 1.;
	btagWeightDown = 1.;
      }

      if(btagWeight != btagWeight) continue;
      if(metCut > 0. && vars[1] >= metCut) continue;

      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
      Float_t btagSFsys = (fabs(btagWeight - btagWeightUp) + fabs(btagWeight - btagWeightDown))/2.;
      Float_t btag_toterr = sqrt(btagWeightErr*btagWeightErr + btagSFsys*btagSFsys);
      Float_t addError2 = puWeight*puWeight*btag_toterr*btag_toterr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      for(unsigned int k = 0; k < vars.size(); k++) {
	if(variables[k] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

	Float_t oldError = mcHistograms[i][k]->GetBinError(mcHistograms[i][k]->FindBin(vars[k]));
	Float_t newerror = sqrt(oldError*oldError + addError2);
	mcHistograms[i][k]->Fill(vars[k]);
	mcHistograms[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);
      }

    }

    for(unsigned int j = 0; j < vars.size(); j++) mcHistograms[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
  }

  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) {
    
    for(int j = 0; j < mcQCDTrees[i]->GetEntries(); j++) {
      mcQCDTrees[i]->GetEntry(j);

      if(nBtagReq == 0) {
	btagWeight = 1.;
	btagWeightErr = 0.;
	btagWeightUp = 1.;
	btagWeightDown = 1.;
      }

      if(btagWeight != btagWeight) continue;
      if(metCut > 0. && vars[1] >= metCut) continue;

      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
      Float_t btagSFsys = (fabs(btagWeight - btagWeightUp) + fabs(btagWeight - btagWeightDown))/2.;
      Float_t btag_toterr = sqrt(btagWeightErr*btagWeightErr + btagSFsys*btagSFsys);
      Float_t addError2 = puWeight*puWeight*btag_toterr*btag_toterr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      for(unsigned int k = 0; k < vars.size(); k++) {
	if(variables[k] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

	Float_t oldError = mcQCDHistograms[i][k]->GetBinError(mcQCDHistograms[i][k]->FindBin(vars[k]));
	Float_t newerror = sqrt(oldError*oldError + addError2);
	mcQCDHistograms[i][k]->Fill(vars[k]);
	mcQCDHistograms[i][k]->SetBinError(mcQCDHistograms[i][k]->FindBin(vars[k]), newerror);
      }

    }

    for(unsigned int j = 0; j < vars.size(); j++) mcQCDHistograms[i][j]->Scale(intLumi_int * crossSections[j] / mcNGen[j]);
  }

  for(int i = 0; i < sigaTree->GetEntries(); i++) {
    sigaTree->GetEntry(i);

    if(nBtagReq == 0) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeight != btagWeight) continue;
    if(metCut > 0. && vars[1] >= metCut) continue;

    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    Float_t btagSFsys = (fabs(btagWeight - btagWeightUp) + fabs(btagWeight - btagWeightDown))/2.;
    Float_t btag_toterr = sqrt(btagWeightErr*btagWeightErr + btagSFsys*btagSFsys);
    Float_t addError2 = puWeight*puWeight*btag_toterr*btag_toterr + btagWeight*btagWeight*puWeightErr*puWeightErr;

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;
      Float_t olderror = h_siga[j]->GetBinError(h_siga[j]->FindBin(vars[j]));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_siga[j]->Fill(vars[j], puWeight * btagWeight);
      h_siga[j]->SetBinError(h_siga[j]->FindBin(vars[j]), newerror);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) h_siga[j]->Scale(intLumi_int * 0.147492 / 15000.);

  for(int i = 0; i < sigbTree->GetEntries(); i++) {
    sigbTree->GetEntry(i);

    if(nBtagReq == 0) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeight != btagWeight) continue;
    if(metCut > 0. && vars[1] >= metCut) continue;

    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    Float_t btagSFsys = (fabs(btagWeight - btagWeightUp) + fabs(btagWeight - btagWeightDown))/2.;
    Float_t btag_toterr = sqrt(btagWeightErr*btagWeightErr + btagSFsys*btagSFsys);
    Float_t addError2 = puWeight*puWeight*btag_toterr*btag_toterr + btagWeight*btagWeight*puWeightErr*puWeightErr;

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;
      Float_t olderror = h_sigb[j]->GetBinError(h_sigb[j]->FindBin(vars[j]));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_sigb[j]->Fill(vars[j], puWeight * btagWeight);
      h_sigb[j]->SetBinError(h_sigb[j]->FindBin(vars[j]), newerror);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) h_sigb[j]->Scale(intLumi_int * 0.0399591 / 15000.);

  ggTree->ResetBranchAddresses();
  qcdTree->ResetBranchAddresses();

  for(unsigned int i = 0; i < mcTrees.size(); i++) mcTrees[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) mcQCDTrees[i]->ResetBranchAddresses();

  sigaTree->ResetBranchAddresses();
  sigbTree->ResetBranchAddresses();

}

void PlotMaker::SubtractMCFromQCD() {

  for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms[i].size(); j++) {
      h_qcd[j]->Add(mcQCDHistograms[i][j], -1.);
    }
  }

}

void PlotMaker::CreatePlot(TString variable,
			   bool divideByWidth,
			   TString xaxisTitle, TString yaxisTitle,
			   Float_t xmin, Float_t xmax,
			   Float_t ymin, Float_t ymax,
			   Float_t ratiomin, Float_t ratiomax,
			   bool drawSignal, bool drawLegend, bool drawPrelim,
			   TFile*& out) {

  unsigned int var_num = 10000;

  for(unsigned int i = 0; i < variables.size(); i++) {
    if(variable == variables[i]) {
      var_num = i;
      break;
    }
  }

  if(var_num > variables.size()) return;

  if(divideByWidth) {
    yaxisTitle = "Number of Events / GeV";

    h_gg[var_num] = (TH1D*)DivideByBinWidth(h_gg[var_num]);
    h_qcd[var_num] = (TH1D*)DivideByBinWidth(h_qcd[var_num]);

    for(unsigned int i = 0; i < mcHistograms.size(); i++) mcHistograms[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms[i][var_num]);
    for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) mcQCDHistograms[i][var_num] = (TH1D*)DivideByBinWidth(mcQCDHistograms[i][var_num]);

    h_siga[var_num] = (TH1D*)DivideByBinWidth(h_siga[var_num]);
    h_sigb[var_num] = (TH1D*)DivideByBinWidth(h_sigb[var_num]);
  }

  DrawPlot(var_num, variable,
	   xaxisTitle, yaxisTitle,
	   xmin, xmax,
	   ymin, ymax,
	   ratiomin, ratiomax,
	   drawSignal, drawLegend, drawPrelim,
	   out);

}

void PlotMaker::DrawPlot(int variableNumber, TString variable,
			 TString xaxisTitle, TString yaxisTitle,
			 Float_t xmin, Float_t xmax,
			 Float_t ymin, Float_t ymax,
			 Float_t ratiomin, Float_t ratiomax,
			 bool drawSignal, bool drawLegend, bool drawPrelim,
			 TFile*& out) {

  out->cd();
  h_gg[variableNumber]->Write();
  h_qcd[variableNumber]->Write();
  for(unsigned int i = 0; i < mcHistograms.size(); i++) mcHistograms[i][variableNumber]->Write();

  // Stack histograms; qcd is on top
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_qcd[variableNumber]->Add(mcHistograms[i][variableNumber]);
    for(unsigned int j = i + 1; j < mcHistograms.size(); j++) {
      mcHistograms[i][variableNumber]->Add(mcHistograms[j][variableNumber]);
    }
  }

  TH1D * bkg = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req);
  bkg->Write();

  Double_t kolm = h_gg[variableNumber]->KolmogorovTest(bkg);
  TString kolmText = Form("KS test probability = %5.3g", kolm);
  TText * tt = new TText(0.92, 0.5, kolmText);
  tt->SetTextAngle(90.);
  tt->SetNDC(); tt->SetTextSize( 0.032 );

  TH1D * errors = (TH1D*)bkg->Clone("errors");

  if(drawSignal) calculateROC(h_siga[variableNumber], h_sigb[variableNumber], bkg, req, variable);

  TLegend * leg = new TLegend(0.50, 0.65, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(h_gg[variableNumber], "#gamma#gamma Candidate Sample", "LP");
  leg->AddEntry(errors, "Total Background Uncertainty", "F");
  leg->AddEntry(bkg, "QCD", "F");
  leg->AddEntry(mcHistograms[0][variableNumber], legendNames[0], "F");
  for(unsigned int i = 1; i < mcHistograms.size(); i++) {
    if(mcLayerNumbers[i] != mcLayerNumbers[i-1]) leg->AddEntry(mcHistograms[i][variableNumber], legendNames[i], "F");
  }
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  TPaveText * prelim = new TPaveText(0.50, 0.42, 0.85, 0.62, "NDC");
  prelim->SetFillColor(0);
  prelim->SetFillStyle(0);
  prelim->SetLineColor(0);
  prelim->AddText("CMS Preliminary 2013");
  prelim->AddText(" ");
  prelim->AddText("#sqrt{s} = 8 TeV, #intL = "+intLumi+" fb^{-1}");
  prelim->AddText(req+" Requirement");

  h_gg[variableNumber]->SetMarkerStyle(20); 
  h_gg[variableNumber]->SetMarkerSize(1.5);

  errors->SetFillColor(kOrange+10);
  errors->SetFillStyle(3154);
  errors->SetMarkerSize(0);

  bkg->SetFillColor(kSpring);
  bkg->SetMarkerSize(0);
  bkg->SetLineColor(1);

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    mcHistograms[i][variableNumber]->SetFillColor(mcLayerColors[i]);
    mcHistograms[i][variableNumber]->SetMarkerSize(0);
    mcHistograms[i][variableNumber]->SetLineColor(1);
  }

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);

  TPad * padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  TPad * padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);

  padhi->Draw();
  padlo->Draw();
  padhi->cd();

  padhi->SetLogy(false);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  //padhi->SetGridx(true);
  //padhi->SetGridy(true);
  padhi->SetBottomMargin(0);

  bkg->SetTitle(variable);
  bkg->GetXaxis()->SetTitle(xaxisTitle);
  bkg->GetYaxis()->SetTitle(yaxisTitle);

  if(xmax > xmin) bkg->GetXaxis()->SetRangeUser(xmin, xmax);
  bkg->GetYaxis()->SetRangeUser(ymin, ymax);

  bkg->Draw("hist");
  mcHistograms[0][variableNumber]->Draw("same hist");
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    if(i != 0 && mcLayerNumbers[i] != mcLayerNumbers[i-1]) mcHistograms[i][variableNumber]->Draw("same hist");
  }
  errors->Draw("same e2");
  h_gg[variableNumber]->Draw("same e1");
  bkg->Draw("same axis");

  if(drawSignal) {
    h_siga[variableNumber]->SetLineColor(kMagenta);
    h_siga[variableNumber]->SetLineWidth(3);
    leg->AddEntry(h_siga[variableNumber], "GGM #gamma#gamma (460_175)", "L");
    h_siga[variableNumber]->Draw("same hist");
    
    h_sigb[variableNumber]->SetLineColor(kBlue);
    h_sigb[variableNumber]->SetLineWidth(3);
    leg->AddEntry(h_sigb[variableNumber], "GGM #gamma#gamma (560_325)", "L");
    h_sigb[variableNumber]->Draw("same hist");
  }

  if(drawLegend) leg->Draw("same");
  if(drawPrelim && drawLegend) prelim->Draw("same");
  if(displayKStest) {
    tt->AppendPad();
    KSscores.push_back(make_pair(variable, kolm));
  }

  padlo->cd();
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  TH1D * ratio = (TH1D*)h_gg[variableNumber]->Clone("ratio");
  ratio->Reset();
  ratio->SetTitle("Data / Background");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(bkg->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, h_gg[variableNumber]->GetBinContent(i+1) / bkg->GetBinContent(i+1));
    ratio->SetBinError(i+1, h_gg[variableNumber]->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  TH1D * ratio_sys;
  ratio_sys = (TH1D*)bkg->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    ratio_sys->SetBinContent(i+1, 1.);
    if(bkg->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, ratio_sys->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  if(xmax > xmin) ratio->GetXaxis()->SetRangeUser(xmin, xmax);

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  ratio->GetXaxis()->SetTitle(xaxisTitle);
  ratio->GetXaxis()->SetLabelFont(63);
  ratio->GetXaxis()->SetLabelSize(48);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(0.6);
  ratio->GetYaxis()->SetTitle("Data / Background");
  ratio->GetYaxis()->SetLabelFont(63);
  ratio->GetYaxis()->SetLabelSize(48);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetRangeUser(ratiomin, ratiomax);
  ratio->GetYaxis()->SetNdivisions(508);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");

  TLine * oneLine = new TLine(xmin, 1, xmax, 1);
  oneLine->SetLineStyle(2);
  oneLine->Draw();  

  padhi->cd();
  padhi->SetLogy(true);
  can->SaveAs(variable+"_"+req+".pdf");

  delete can;

}

void PlotMaker::CreateTable() {
  /*
    const int nBins = 5;
    Double_t xbins[nBins+1] = {0, 20, 50, 80, 100, 1000};

    Double_t rangeLow[nBins] = {0, 0, 50, 80, 100};
    Double_t rangeHigh[nBins] = {20, 50, -1, -1, -1};

    TH1D * h_gg = HistoFromTree("pfMET", ggTree, "pfMet2_gg_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * h_ewk = HistoFromTree("pfMET", egTree, "pfMet2_eg_"+req, "pfMet2", nBins, xbins, -1.);
  
    TH1D * ewk_noNorm = (TH1D*)h_ewk->Clone();
    h_ewk->Scale(egScale);
    for(int i = 0; i < h_ewk->GetNbinsX(); i++) {
    Float_t normerr = egScaleErr*(ewk_noNorm->GetBinContent(i+1));
    Float_t staterr = h_ewk->GetBinError(i+1);
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    h_ewk->SetBinError(i+1, new_err);
    }

    TH1D * qcd30to40 = SignalHistoFromTree(intLumi_int * 5.195E7 * 2.35E-4 * 1.019 * 1.019 / 6061407., true, "pfMET", qcd30to40Tree, "pfMet2_qcd30to40_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * qcd40 = SignalHistoFromTree(intLumi_int * 5.195E7 * 0.002175 * 1.019 * 1.019 / 9782735., true, "pfMET", qcd40Tree, "pfMet2_qcd40_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * h_qcd = (TH1D*)qcd30to40->Clone("pfMet2_qcd_"+req);
    h_qcd->Add(qcd40);

    TH1D * gjet20to40 = SignalHistoFromTree(intLumi_int * 81930.0 * 0.001835 * 1.019 * 1.019 / 5907942., true, "pfMET", gjet20to40Tree, "pfMet2_gjet20to40_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * gjet40 = SignalHistoFromTree(intLumi_int * 8884.0 * 0.05387 * 1.019 * 1.019 / 5956149., true, "pfMET", gjet40Tree, "pfMet2_gjet40_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * h_gjet = (TH1D*)gjet20to40->Clone("pfMet2_gjet_"+req);
    h_gjet->Add(gjet40);

    TH1D * diphotonjets = SignalHistoFromTree(intLumi_int * 75.39 * 1.019 * 1.019 / 1156030., true, "pfMET", diphotonjetsTree, "pfMet2_diphotonjets_"+req, "pfMet2", nBins, xbins, -1.);

    TH1D * ttHadronic = SignalHistoFromTree(intLumi_int * 53.4 * 1.019 * 1.019 / 10537444., true, "pfMET", ttHadronicTree, "pfMet2_ttHadronic_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * ttSemiLep = SignalHistoFromTree(intLumi_int * 53.2 * 1.019 * 1.019 / 25424818., true, "pfMET", ttSemiLepTree, "pfMet2_ttSemiLep_"+req, "pfMet2", nBins, xbins, -1.);
    TH1D * ttbar = (TH1D*)ttHadronic->Clone("pfMet2_ttbar_"+req);
    ttbar->Add(ttSemiLep);

    TH1D * h_ttg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 14.0 / 1719954., true, "pfMET", ttgjetsTree, "pfMet2_ttgjets_"+req, "pfMet2", nBins, xbins, -1.);

    TH1D * bkg = (TH1D*)h_qcd->Clone("pfMet2_bkg_"+req);

    bkg->Add(h_gjet);
    bkg->Add(diphotonjets);
    bkg->Add(h_ewk);
    bkg->Add(h_ttg);

    // Calculate entries

    FILE * tableFile = fopen("errorTable_"+req+".temp", "w");

    Double_t binLow[nBins], binHigh[nBins];
    for(int i = 0; i < nBins; i++) {
    binLow[i] = h_gg->GetXaxis()->FindBin(rangeLow[i]);
    binHigh[i] = (rangeHigh[i] == -1) ? -1 : h_gg->GetXaxis()->FindBin(rangeHigh[i]) - 1;
    }

    for(int i = 0; i < nBins; i++) {

    Double_t gg, ggerr;
    gg = h_gg->IntegralAndError(binLow[i], binHigh[i], ggerr);
    fprintf(tableFile, "ggval%dx:%.0f\nggstat%dx:%.1f\n", i+1, gg, i+1, ggerr);

    Double_t qcd, qcderr;
    qcd = h_qcd->IntegralAndError(binLow[i], binHigh[i], qcderr);
    fprintf(tableFile, "qcdval%dx:%.0f\nqcdstat%dx:%.1f\n", i+1, qcd, i+1, qcderr);

    Double_t gjet, gjeterr;
    gjet = h_gjet->IntegralAndError(binLow[i], binHigh[i], gjeterr);
    fprintf(tableFile, "gjetval%dx:%.0f\ngjetstat%dx:%.1f\n", i+1, gjet, i+1, gjeterr);

    Double_t dipho, diphoerr;
    dipho = diphotonjets->IntegralAndError(binLow[i], binHigh[i], diphoerr);
    fprintf(tableFile, "diphoval%dx:%.0f\ndiphostat%dx:%.1f\n", i+1, dipho, i+1, diphoerr);

    Double_t ewk, ewkerr;
    ewk = h_ewk->IntegralAndError(binLow[i], binHigh[i], ewkerr);
    fprintf(tableFile, "ewkval%dx:%.1f\newkstat%dx:%.2f\n", i+1, ewk, i+1, ewkerr);
    
    Double_t ttgg, ttggerr;
    ttgg = h_ttg->IntegralAndError(binLow[i], binHigh[i], ttggerr);
    fprintf(tableFile, "ttgval%dx:%.1f\nttgstat%dx:%.2f\n", i+1, ttgg, i+1, ttggerr);

    Double_t total, totalerr;
    total = bkg->IntegralAndError(binLow[i], binHigh[i], totalerr);
    fprintf(tableFile, "bkgval%dx:%.1f\nbkgstat%dx:%.2f\n", i+1, total, i+1, totalerr);

    }

    fclose(tableFile);
  */
}

void PlotMaker::PlotKolmogorovValues() {

  if(!displayKStest) return;

  TH1D * h_ks = new TH1D("h_ks_"+req, "h_ks_"+req, (int)KSscores.size(), 0, (int)KSscores.size());

  for(unsigned int i = 0; i < KSscores.size(); i++) {
    h_ks->SetBinContent(i+1, KSscores[i].second);
    h_ks->GetXaxis()->SetBinLabel(i+1, KSscores[i].first);
  }

  TCanvas * can = new TCanvas("ks_can_"+req, "Plot", 10, 10, 2000, 2000);
  can->SetLogy(true);

  h_ks->Draw("hist");
  can->SaveAs("ksScores_"+req+".pdf");

  delete can;

}

void prep_signal(TString req, int nPhotons_req) {

  Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst[i] + mst[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino[i] + mBino[i-1])/2.;
  ybins[32] = 2175;

  char code[100];
  int index1, index2;

  TH2D * h_acc = new TH2D("acc_"+req, "acc_"+req, 30, xbins, 32, ybins);
  
  TFile * out = new TFile("signal_"+req+".root", "RECREATE");

  for(int i = 0; i < 899; i++) {

    index1 = mst[int(i)/31];
    index2 = mBino[int(i)%31];
    sprintf(code, "_mst_%d_m1_%d", index1, index2);
    TString code_t = code;

    TFile * f = new TFile("../acceptance/signal_contamination"+code_t+".root", "READ");
    if(f->IsZombie()) {
      f->Close();
      continue;
    }
    
    TTree * ggTree = (TTree*)f->Get(req+"_signalTree");
    
    TH1D * gg;
    
    if(ggTree->GetEntries() > 0) {
      gg = (TH1D*)SignalHistoFromTree(1.0, "pfMET", ggTree, "met_gg_"+req+code_t, "met_gg_"+req+code_t, 400, 0., 2000., -1, nPhotons_req);
      
      out->cd();
      gg->Write();
    }
    else {
      f->Close();
      continue;
    }

    double n = 15000.;

    double acceptance = gg->Integral();
    if(n > 0) h_acc->Fill(index1, index2, acceptance / n);

    f->Close();
  }

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  h_acc->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_acc->GetXaxis()->SetRangeUser(0, 1600);
  h_acc->GetXaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_acc->GetYaxis()->SetTitleOffset(1.3);
  h_acc->GetYaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetRangeUser(0, 1600);
  h_acc->GetZaxis()->SetLabelSize(0.02);
  h_acc->Draw("colz");
  can->SaveAs("acceptance_"+req+".pdf");
  
  out->Write();
  out->Close();

}

void plotReducedChi2(vector<TH1D*> gg, vector<TH1D*> gf, vector<TH1D*> ff,
		     TH2D*& gf_chi2,
		     TH2D*& ff_chi2,
		     Int_t binx) {

  for(unsigned int i = 0; i < gg.size(); i++) {
    if(gg[i]->Integral() >= 1.) gg[i]->Scale(1./gg[i]->Integral());
  }
  for(unsigned int i = 0; i < gf.size(); i++) {
    if(gf[i]->Integral() >= 1.) gf[i]->Scale(1./gf[i]->Integral());
  }
  for(unsigned int i = 0; i < ff.size(); i++) {
    if(ff[i]->Integral() >= 1.) ff[i]->Scale(1./ff[i]->Integral());
  }

  for(int i = 0; i < 10; i++) {

    Float_t chi2 = 0.;
    Int_t nBins = 0;
    
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - gf[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + gf[i]->GetBinError(j+1)*gf[i]->GetBinError(j+1);

      if(val_den == 0.) continue;

      chi2 += val_num * val_num / val_den;
      nBins++;
    }
    chi2 /= nBins;
    gf_chi2->SetBinContent(gf_chi2->FindBin(binx, i), chi2);

    chi2 = 0.;
    nBins = 0;
    
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - ff[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + ff[i]->GetBinError(j+1)*ff[i]->GetBinError(j+1);

      if(val_den == 0.) continue;

      chi2 += val_num * val_num / val_den;
      nBins++;
    }
    chi2 /= nBins;
    ff_chi2->SetBinContent(gf_chi2->FindBin(binx, i), chi2);
  }

  // now fill bin 10, chi2 across all variables
  Float_t chi2_all = 0.;
  Float_t nBins_all = 0;
  for(unsigned int i = 0; i < gg.size(); i++) {
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - gf[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + gf[i]->GetBinError(j+1)*gf[i]->GetBinError(j+1);

      if(val_den == 0.) continue;

      chi2_all += val_num * val_num / val_den;
      nBins_all++;
    }
  }
  chi2_all /= nBins_all;
  gf_chi2->SetBinContent(gf_chi2->FindBin(binx, 10), chi2_all);

  chi2_all = 0.;
  nBins_all = 0;
  for(unsigned int i = 0; i < gg.size(); i++) {
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - ff[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + ff[i]->GetBinError(j+1)*ff[i]->GetBinError(j+1);
    
      if(val_den == 0.) continue;
  
      chi2_all += val_num * val_num / val_den;
      nBins_all++;
    }
  }
  chi2_all /= nBins_all;
  ff_chi2->SetBinContent(ff_chi2->FindBin(binx, 10), chi2_all);

}
