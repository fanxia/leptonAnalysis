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
    srootb_a->SetBinContent(i+1, s / sqrt(b));
    Double_t error = (s / sqrt(b)) * sqrt(serr*serr/(s*s) + berr*berr/(4*b*b));
    srootb_a->SetBinError(i+1, error);

    s = sig_b->IntegralAndError(i+1, -1, serr, "");
    x_b[i] = s / sig_b->Integral();
    y_b[i] = b / bkg->Integral();
    srootb_b->SetBinContent(i+1, s / sqrt(b));
    error = (s / sqrt(b)) * sqrt(serr*serr/(s*s) + berr*berr/(4*b*b));
    srootb_b->SetBinContent(i+1, error);

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
  srootb_a->Draw("e1");
  srootb_b->Draw("e1 same");

  float lineMaxY = (srootb_b->GetMaximum() > srootb_a->GetMaximum()) ? 1.1 * srootb_b->GetMaximum() : 1.1 * srootb_a->GetMaximum();

  TLine * bestCutLine_a = new TLine(srootb_a->GetBinLowEdge(bestCutBin_a), 0, srootb_a->GetBinLowEdge(bestCutBin_a), lineMaxY);
  bestCutLine_a->SetLineColor(kMagenta);
  bestCutLine_a->SetLineWidth(2);
  //bestCutLine_a->Draw("same");
  TLine * bestCutLine_b = new TLine(srootb_b->GetBinLowEdge(bestCutBin_b), 0, srootb_b->GetBinLowEdge(bestCutBin_b), lineMaxY);
  bestCutLine_b->SetLineColor(kBlue);
  bestCutLine_b->SetLineWidth(2);
  //bestCutLine_b->Draw("same");

  canv->SaveAs("roc_"+title+"_"+req+".pdf");

}

class PlotMaker : public TObject {
  
  ClassDef(PlotMaker, 1);

 public:
  PlotMaker(Int_t lumi,
	    TString requirement,
	    bool blind);

  ~PlotMaker();
  
  void LoadLeptonSFs(TString fileName);
  void LoadPhotonSFs(TString fileName);

  bool LoadMCBackground(TString fileName, TString scanName,
			Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfErrorUp, Double_t pdfErrorDown,
			bool removeTTA, bool reweightTop,
			int channel, int layer, int color, TString legendEntry, TString tableEntry);
  
  void SetTrees(TTree * gg, TTree * qcd,
		TTree * sig_a, TTree * sig_a_JECup, TTree * sig_a_JECdown,
		TTree * sig_b, TTree * sig_b_JECup, TTree * sig_b_JECdown);

  void SetDisplayKStest(bool v) { displayKStest = v; }
  void SetUseWHIZARD(bool v) { useWHIZARD = v; }

  void BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi);
  void BookHistogram(TString variable, Int_t nBins, Double_t* customBins);

  void FillHistograms(double metCut, int nPhotons_req, int nBtagReq, int chan);
  void SubtractMCFromQCD();
  void NormalizeQCD();

  void CreatePlot(TString variable, bool divideByWidth, bool needsQCD, TString xaxisTitle, TString yaxisTitle,
		  Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
		  Float_t ratiomin, Float_t ratiomax,
		  bool drawSignal, bool drawLegend, bool drawPrelim,
		  TFile*& out);

  void DrawPlot(int variableNumber, TString variable, bool needsQCD,
		TString xaxisTitle, TString yaxisTitle,
		Float_t xmin, Float_t xmax,
		Float_t ymin, Float_t ymax,
		Float_t ratiomin, Float_t ratiomax,
		bool drawSignal, bool drawLegend, bool drawPrelim,
		TFile*& out);

  void CreateTable();

  void PlotKolmogorovValues();

  void GetLeptonSF(vector<Float_t> vars, int chan, Float_t& central, Float_t& up, Float_t& down);
  void GetPhotonSF(vector<Float_t> vars, Float_t& central, Float_t& up, Float_t& down);

 private:

  vector<TFile*> mcFiles;
  vector<TTree*> mcTrees;
  vector<TTree*> mcTrees_JECup;
  vector<TTree*> mcTrees_JECdown;
  vector<TTree*> mcQCDTrees;
  vector<Double_t> mcNGen;
  vector<Double_t> crossSections;
  vector<Double_t> scaleErrUp;
  vector<Double_t> scaleErrDown;
  vector<Double_t> pdfErrUp;
  vector<Double_t> pdfErrDown;
  vector<int> mcLayerNumbers;
  vector<int> mcLayerColors;
  vector<TString> mcNames;
  vector<TString> legendNames;
  vector<TString> tableNames;
  vector<bool> removeTTAoverlap;
  vector<bool> reweightTopPt;

  vector< vector<TH1D*> > mcHistograms;
  vector< vector<TH1D*> > mcHistograms_btagWeightUp;
  vector< vector<TH1D*> > mcHistograms_btagWeightDown;
  vector< vector<TH1D*> > mcHistograms_scaleUp;
  vector< vector<TH1D*> > mcHistograms_scaleDown;
  vector< vector<TH1D*> > mcHistograms_pdfUp;
  vector< vector<TH1D*> > mcHistograms_pdfDown;
  vector< vector<TH1D*> > mcHistograms_topPtUp;
  vector< vector<TH1D*> > mcHistograms_topPtDown;
  vector< vector<TH1D*> > mcHistograms_JECup;
  vector< vector<TH1D*> > mcHistograms_JECdown;
  vector< vector<TH1D*> > mcHistograms_leptonSFup;
  vector< vector<TH1D*> > mcHistograms_leptonSFdown;
  vector< vector<TH1D*> > mcHistograms_photonSFup;
  vector< vector<TH1D*> > mcHistograms_photonSFdown;

  vector< vector<TH1D*> > mcQCDHistograms;

  TTree * ggTree;
  TTree * qcdTree;

  TTree * sigaTree;
  TTree * sigaTree_JECup;
  TTree * sigaTree_JECdown;

  TTree * sigbTree;
  TTree * sigbTree_JECup;
  TTree * sigbTree_JECdown;

  vector<TString> variables;

  vector<TH1D*> h_gg;
  vector<TH1D*> h_qcd;

  vector<TH1D*> h_siga;
  vector<TH1D*> h_siga_btagWeightUp;
  vector<TH1D*> h_siga_btagWeightDown;
  vector<TH1D*> h_siga_topPtUp;
  vector<TH1D*> h_siga_topPtDown;
  vector<TH1D*> h_siga_JECup;
  vector<TH1D*> h_siga_JECdown;
  vector<TH1D*> h_siga_leptonSFup;
  vector<TH1D*> h_siga_leptonSFdown;
  vector<TH1D*> h_siga_photonSFup;
  vector<TH1D*> h_siga_photonSFdown;

  vector<TH1D*> h_sigb;
  vector<TH1D*> h_sigb_btagWeightUp;
  vector<TH1D*> h_sigb_btagWeightDown;
  vector<TH1D*> h_sigb_topPtUp;
  vector<TH1D*> h_sigb_topPtDown;
  vector<TH1D*> h_sigb_JECup;
  vector<TH1D*> h_sigb_JECdown;
  vector<TH1D*> h_sigb_leptonSFup;
  vector<TH1D*> h_sigb_leptonSFdown;
  vector<TH1D*> h_sigb_photonSFup;
  vector<TH1D*> h_sigb_photonSFdown;

  vector<pair<TString, double> > KSscores;

  TFile * fLeptonSF;
  TH2D * sf_muon;
  TH2D * sf_electron;
  TH2D * sf_SingleElectronTrigger;

  TFile * fPhotonSF;
  TH2D * sf_photon_id;
  TH2D * sf_photon_veto;

  Int_t intLumi_int;
  TString intLumi;
  TString req;

  bool displayKStest;
  bool blinded;

  bool useWHIZARD;

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
  useWHIZARD = false;

  KSscores.clear();

  variables.clear();

  h_gg.clear();
  h_qcd.clear();

  h_siga.clear();
  h_siga_btagWeightUp.clear();
  h_siga_btagWeightDown.clear();
  h_siga_topPtUp.clear();
  h_siga_topPtDown.clear();
  h_siga_JECup.clear();
  h_siga_JECdown.clear();
  h_siga_leptonSFup.clear();
  h_siga_leptonSFdown.clear();
  h_siga_photonSFup.clear();
  h_siga_photonSFdown.clear();

  h_sigb.clear();
  h_sigb_btagWeightUp.clear();
  h_sigb_btagWeightDown.clear();
  h_sigb_topPtUp.clear();
  h_sigb_topPtDown.clear();
  h_sigb_JECup.clear();
  h_sigb_JECdown.clear();
  h_sigb_leptonSFup.clear();
  h_sigb_leptonSFdown.clear();
  h_sigb_photonSFup.clear();
  h_sigb_photonSFdown.clear();

  mcFiles.clear();
  mcTrees.clear();
  mcTrees_JECup.clear();
  mcTrees_JECdown.clear();
  mcQCDTrees.clear();
  mcNGen.clear();
  crossSections.clear();
  scaleErrUp.clear();
  scaleErrDown.clear();
  pdfErrUp.clear();
  pdfErrDown.clear();
  mcLayerNumbers.clear();
  mcLayerColors.clear();
  mcNames.clear();
  legendNames.clear();
  tableNames.clear();
  removeTTAoverlap.clear();
  reweightTopPt.clear();

  mcHistograms.clear();
  mcHistograms_btagWeightUp.clear();
  mcHistograms_btagWeightDown.clear();
  mcHistograms_scaleUp.clear();
  mcHistograms_scaleDown.clear();
  mcHistograms_pdfUp.clear();
  mcHistograms_pdfDown.clear();
  mcHistograms_topPtUp.clear();
  mcHistograms_topPtDown.clear();
  mcHistograms_JECup.clear();
  mcHistograms_JECdown.clear();
  mcHistograms_leptonSFup.clear();
  mcHistograms_leptonSFdown.clear();
  mcHistograms_photonSFup.clear();
  mcHistograms_photonSFdown.clear();

  mcQCDHistograms.clear();

}

PlotMaker::~PlotMaker() { 

    KSscores.clear();
    variables.clear();

    delete ggTree;
    delete qcdTree;

    mcHistograms.clear();
    mcHistograms_btagWeightUp.clear();
    mcHistograms_btagWeightDown.clear();
    mcHistograms_scaleUp.clear();
    mcHistograms_scaleDown.clear();
    mcHistograms_pdfUp.clear();
    mcHistograms_pdfDown.clear();
    mcHistograms_topPtUp.clear();
    mcHistograms_topPtDown.clear();
    mcHistograms_JECup.clear();
    mcHistograms_JECdown.clear();
    mcHistograms_leptonSFup.clear();
    mcHistograms_leptonSFdown.clear();
    mcHistograms_photonSFup.clear();
    mcHistograms_photonSFdown.clear();
    mcQCDHistograms.clear();
    mcTrees.clear();
    mcTrees_JECup.clear();
    mcTrees_JECdown.clear();
    mcQCDTrees.clear();
    mcFiles.clear();
    mcNGen.clear();
    crossSections.clear();
    scaleErrUp.clear();
    scaleErrDown.clear();
    pdfErrUp.clear();
    pdfErrDown.clear();
    mcLayerNumbers.clear();
    mcLayerColors.clear();
    mcNames.clear();
    legendNames.clear();
    tableNames.clear();
    removeTTAoverlap.clear();
    reweightTopPt.clear();

    fLeptonSF->Close();
    fPhotonSF->Close();

    delete sigaTree;
    delete sigaTree_JECup;
    delete sigaTree_JECdown;

    delete sigbTree;
    delete sigbTree_JECup;
    delete sigbTree_JECdown;
    
    h_gg.clear();
    h_qcd.clear();

    h_siga.clear();
    h_siga_btagWeightUp.clear();
    h_siga_btagWeightDown.clear();
    h_siga_topPtUp.clear();
    h_siga_topPtDown.clear();
    h_siga_JECup.clear();
    h_siga_JECdown.clear();
    h_siga_leptonSFup.clear();
    h_siga_leptonSFdown.clear();
    h_siga_photonSFup.clear();
    h_siga_photonSFdown.clear();
    
    h_sigb.clear();
    h_sigb_btagWeightUp.clear();
    h_sigb_btagWeightDown.clear();
    h_sigb_topPtUp.clear();
    h_sigb_topPtDown.clear();
    h_sigb_JECup.clear();
    h_sigb_JECdown.clear();
    h_sigb_leptonSFup.clear();
    h_sigb_leptonSFdown.clear();
    h_sigb_photonSFup.clear();
    h_sigb_photonSFdown.clear();

}

void PlotMaker::SetTrees(TTree * gg, TTree * qcd,
			 TTree * sig_a, TTree * sig_a_JECup, TTree * sig_a_JECdown,
			 TTree * sig_b, TTree * sig_b_JECup, TTree * sig_b_JECdown) {

  ggTree = gg;
  qcdTree = qcd;
  
  sigaTree = sig_a;
  sigaTree_JECup = sig_a_JECup;
  sigaTree_JECdown = sig_a_JECdown;

  sigbTree = sig_b;
  sigbTree_JECup = sig_b_JECup;
  sigbTree_JECdown = sig_b_JECdown;

}

void PlotMaker::LoadLeptonSFs(TString fileName) {
  
  fLeptonSF = new TFile(fileName, "READ");
  sf_electron = (TH2D*)fLeptonSF->Get("TightEleIdIsoSF");
  sf_SingleElectronTrigger = (TH2D*)fLeptonSF->Get("TightEleTriggerSF");

  sf_muon = (TH2D*)fLeptonSF->Get("mu_pt_eta_full_id_iso_hlt_8TeV");

}

void PlotMaker::LoadPhotonSFs(TString fileName) {
  
  fPhotonSF = new TFile(fileName, "READ");
  sf_photon_id = (TH2D*)fPhotonSF->Get("PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
  sf_photon_veto = (TH2D*)fPhotonSF->Get("PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");

}

bool PlotMaker::LoadMCBackground(TString fileName, TString scanName,
				 Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfErrorUp, Double_t pdfErrorDown,
				 bool removeTTA, bool reweightTop,
				 int channel, int layer, int color, TString legendEntry, TString tableEntry) {

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

  mcTrees_JECup.push_back((TTree*)mcFiles.back()->Get(channels[channel]+"_signalTree_JECup"));
  if(!mcTrees_JECup.back()) {
    cout << "Could not load TTree " << channels[channel] << "_signalTree_JECup from TFile " << fileName << endl;
    return false;
  }

  mcTrees_JECdown.push_back((TTree*)mcFiles.back()->Get(channels[channel]+"_signalTree_JECdown"));
  if(!mcTrees_JECdown.back()) {
    cout << "Could not load TTree " << channels[channel] << "_signalTree_JECdown from TFile " << fileName << endl;
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
  pdfErrUp.push_back(pdfErrorUp);
  pdfErrDown.push_back(pdfErrorDown);  
  mcLayerNumbers.push_back(layer);
  mcLayerColors.push_back(color);
  mcNames.push_back(scanName);
  legendNames.push_back(legendEntry);
  tableNames.push_back(tableEntry);
  removeTTAoverlap.push_back(removeTTA);
  reweightTopPt.push_back(reweightTop);

  mcHistograms.resize(mcHistograms.size() + 1);
  mcHistograms_btagWeightUp.resize(mcHistograms_btagWeightUp.size() + 1);
  mcHistograms_btagWeightDown.resize(mcHistograms_btagWeightDown.size() + 1);
  mcHistograms_scaleUp.resize(mcHistograms_scaleUp.size() + 1);
  mcHistograms_scaleDown.resize(mcHistograms_scaleDown.size() + 1);
  mcHistograms_pdfUp.resize(mcHistograms_pdfUp.size() + 1);
  mcHistograms_pdfDown.resize(mcHistograms_pdfDown.size() + 1);
  mcHistograms_topPtUp.resize(mcHistograms_topPtUp.size() + 1);
  mcHistograms_topPtDown.resize(mcHistograms_topPtDown.size() + 1);
  mcHistograms_JECup.resize(mcHistograms_JECup.size() + 1);
  mcHistograms_JECdown.resize(mcHistograms_JECdown.size() + 1);
  mcHistograms_leptonSFup.resize(mcHistograms_leptonSFup.size() + 1);
  mcHistograms_leptonSFdown.resize(mcHistograms_leptonSFdown.size() + 1);
  mcHistograms_photonSFup.resize(mcHistograms_photonSFup.size() + 1);
  mcHistograms_photonSFdown.resize(mcHistograms_photonSFdown.size() + 1);
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

    TH1D * h_bkg_btagWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightUp");
    mcHistograms_btagWeightUp[i].push_back(h_bkg_btagWeightUp);

    TH1D * h_bkg_btagWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightDown");
    mcHistograms_btagWeightDown[i].push_back(h_bkg_btagWeightDown);

    TH1D * h_bkg_scaleUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleUp");
    mcHistograms_scaleUp[i].push_back(h_bkg_scaleUp);

    TH1D * h_bkg_scaleDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleDown");
    mcHistograms_scaleDown[i].push_back(h_bkg_scaleDown);

    TH1D * h_bkg_pdfUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfUp");
    mcHistograms_pdfUp[i].push_back(h_bkg_pdfUp);

    TH1D * h_bkg_pdfDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfDown");
    mcHistograms_pdfDown[i].push_back(h_bkg_pdfDown);

    TH1D * h_bkg_topPtUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtUp");
    mcHistograms_topPtUp[i].push_back(h_bkg_topPtUp);

    TH1D * h_bkg_topPtDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtDown");
    mcHistograms_topPtDown[i].push_back(h_bkg_topPtDown);

    TH1D * h_bkg_JECup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECup");
    mcHistograms_JECup[i].push_back(h_bkg_JECup);

    TH1D * h_bkg_JECdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECdown");
    mcHistograms_JECdown[i].push_back(h_bkg_JECdown);

    TH1D * h_bkg_leptonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFup");
    mcHistograms_leptonSFup[i].push_back(h_bkg_leptonSFup);

    TH1D * h_bkg_leptonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFdown");
    mcHistograms_leptonSFdown[i].push_back(h_bkg_leptonSFdown);

    TH1D * h_bkg_photonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFup");
    mcHistograms_photonSFup[i].push_back(h_bkg_photonSFup);

    TH1D * h_bkg_photonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFdown");
    mcHistograms_photonSFdown[i].push_back(h_bkg_photonSFdown);
  }

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+req, variable, nBins, xlo, xhi);
    h_bkg->Sumw2();
    mcQCDHistograms[i].push_back(h_bkg);
  }

  TH1D * sig_a = new TH1D(variable+"_a_"+req, variable, nBins, xlo, xhi);
  sig_a->Sumw2();
  h_siga.push_back(sig_a);

  TH1D * sig_a_btagWeightUp = new TH1D(variable+"_a_"+req+"_btagWeightUp", variable, nBins, xlo, xhi);
  sig_a_btagWeightUp->Sumw2();
  h_siga_btagWeightUp.push_back(sig_a_btagWeightUp);

  TH1D * sig_a_btagWeightDown = new TH1D(variable+"_a_"+req+"_btagWeightDown", variable, nBins, xlo, xhi);
  sig_a_btagWeightDown->Sumw2();
  h_siga_btagWeightDown.push_back(sig_a_btagWeightDown);

  TH1D * sig_a_topPtUp = new TH1D(variable+"_a_"+req+"_topPtUp", variable, nBins, xlo, xhi);
  sig_a_topPtUp->Sumw2();
  h_siga_topPtUp.push_back(sig_a_topPtUp);

  TH1D * sig_a_topPtDown = new TH1D(variable+"_a_"+req+"_topPtDown", variable, nBins, xlo, xhi);
  sig_a_topPtDown->Sumw2();
  h_siga_topPtDown.push_back(sig_a_topPtDown);
  
  TH1D * sig_a_JECup = new TH1D(variable+"_a_"+req+"_JECup", variable, nBins, xlo, xhi);
  sig_a_JECup->Sumw2();
  h_siga_JECup.push_back(sig_a_JECup);

  TH1D * sig_a_JECdown = new TH1D(variable+"_a_"+req+"_JECdown", variable, nBins, xlo, xhi);
  sig_a_JECdown->Sumw2();
  h_siga_JECdown.push_back(sig_a_JECdown);

  TH1D * sig_a_leptonSFup = new TH1D(variable+"_a_"+req+"_leptonSFup", variable, nBins, xlo, xhi);
  sig_a_leptonSFup->Sumw2();
  h_siga_leptonSFup.push_back(sig_a_leptonSFup);

  TH1D * sig_a_leptonSFdown = new TH1D(variable+"_a_"+req+"_leptonSFdown", variable, nBins, xlo, xhi);
  sig_a_leptonSFdown->Sumw2();
  h_siga_leptonSFdown.push_back(sig_a_leptonSFdown);

  TH1D * sig_a_photonSFup = new TH1D(variable+"_a_"+req+"_photonSFup", variable, nBins, xlo, xhi);
  sig_a_photonSFup->Sumw2();
  h_siga_photonSFup.push_back(sig_a_photonSFup);

  TH1D * sig_a_photonSFdown = new TH1D(variable+"_a_"+req+"_photonSFdown", variable, nBins, xlo, xhi);
  sig_a_photonSFdown->Sumw2();
  h_siga_photonSFdown.push_back(sig_a_photonSFdown);

  TH1D * sig_b = new TH1D(variable+"_b_"+req, variable, nBins, xlo, xhi);
  sig_b->Sumw2();
  h_sigb.push_back(sig_b);

  TH1D * sig_b_btagWeightUp = new TH1D(variable+"_b_"+req+"_btagWeightUp", variable, nBins, xlo, xhi);
  sig_b_btagWeightUp->Sumw2();
  h_sigb_btagWeightUp.push_back(sig_b_btagWeightUp);

  TH1D * sig_b_btagWeightDown = new TH1D(variable+"_b_"+req+"_btagWeightDown", variable, nBins, xlo, xhi);
  sig_b_btagWeightDown->Sumw2();
  h_sigb_btagWeightDown.push_back(sig_b_btagWeightDown);

  TH1D * sig_b_topPtUp = new TH1D(variable+"_b_"+req+"_topPtUp", variable, nBins, xlo, xhi);
  sig_b_topPtUp->Sumw2();
  h_sigb_topPtUp.push_back(sig_b_topPtUp);

  TH1D * sig_b_topPtDown = new TH1D(variable+"_b_"+req+"_topPtDown", variable, nBins, xlo, xhi);
  sig_b_topPtDown->Sumw2();
  h_sigb_topPtDown.push_back(sig_b_topPtDown);
  
  TH1D * sig_b_JECup = new TH1D(variable+"_b_"+req+"_JECup", variable, nBins, xlo, xhi);
  sig_b_JECup->Sumw2();
  h_sigb_JECup.push_back(sig_b_JECup);

  TH1D * sig_b_JECdown = new TH1D(variable+"_b_"+req+"_JECdown", variable, nBins, xlo, xhi);
  sig_b_JECdown->Sumw2();
  h_sigb_JECdown.push_back(sig_b_JECdown);

  TH1D * sig_b_leptonSFup = new TH1D(variable+"_b_"+req+"_leptonSFup", variable, nBins, xlo, xhi);
  sig_b_leptonSFup->Sumw2();
  h_sigb_leptonSFup.push_back(sig_b_leptonSFup);

  TH1D * sig_b_leptonSFdown = new TH1D(variable+"_b_"+req+"_leptonSFdown", variable, nBins, xlo, xhi);
  sig_b_leptonSFdown->Sumw2();
  h_sigb_leptonSFdown.push_back(sig_b_leptonSFdown);
  
  TH1D * sig_b_photonSFup = new TH1D(variable+"_b_"+req+"_photonSFup", variable, nBins, xlo, xhi);
  sig_b_photonSFup->Sumw2();
  h_sigb_photonSFup.push_back(sig_b_photonSFup);

  TH1D * sig_b_photonSFdown = new TH1D(variable+"_b_"+req+"_photonSFdown", variable, nBins, xlo, xhi);
  sig_b_photonSFdown->Sumw2();
  h_sigb_photonSFdown.push_back(sig_b_photonSFdown);

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

    TH1D * h_bkg_btagWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightUp");
    mcHistograms_btagWeightUp[i].push_back(h_bkg_btagWeightUp);

    TH1D * h_bkg_btagWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightDown");
    mcHistograms_btagWeightDown[i].push_back(h_bkg_btagWeightDown);

    TH1D * h_bkg_scaleUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleUp");
    mcHistograms_scaleUp[i].push_back(h_bkg_scaleUp);

    TH1D * h_bkg_scaleDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleDown");
    mcHistograms_scaleDown[i].push_back(h_bkg_scaleDown);

    TH1D * h_bkg_pdfUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfUp");
    mcHistograms_pdfUp[i].push_back(h_bkg_pdfUp);

    TH1D * h_bkg_pdfDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfDown");
    mcHistograms_pdfDown[i].push_back(h_bkg_pdfDown);
    
    TH1D * h_bkg_topPtUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtUp");
    mcHistograms_topPtUp[i].push_back(h_bkg_topPtUp);

    TH1D * h_bkg_topPtDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtDown");
    mcHistograms_topPtDown[i].push_back(h_bkg_topPtDown);

    TH1D * h_bkg_JECup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECup");
    mcHistograms_JECup[i].push_back(h_bkg_JECup);

    TH1D * h_bkg_JECdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECdown");
    mcHistograms_JECdown[i].push_back(h_bkg_JECdown);

    TH1D * h_bkg_leptonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFup");
    mcHistograms_leptonSFup[i].push_back(h_bkg_leptonSFup);

    TH1D * h_bkg_leptonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFdown");
    mcHistograms_leptonSFdown[i].push_back(h_bkg_leptonSFdown);

    TH1D * h_bkg_photonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFup");
    mcHistograms_photonSFup[i].push_back(h_bkg_photonSFup);

    TH1D * h_bkg_photonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFdown");
    mcHistograms_photonSFdown[i].push_back(h_bkg_photonSFdown);
  }
  
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+req, variable, nBins, customBins);
    h_bkg->Sumw2();
    mcQCDHistograms[i].push_back(h_bkg);
  }

  TH1D * sig_a = new TH1D(variable+"_a_"+req, variable, nBins, customBins);
  sig_a->Sumw2();
  h_siga.push_back(sig_a);

  TH1D * sig_a_btagWeightUp = new TH1D(variable+"_a_"+req+"_btagWeightUp", variable, nBins, customBins);
  sig_a_btagWeightUp->Sumw2();
  h_siga_btagWeightUp.push_back(sig_a_btagWeightUp);

  TH1D * sig_a_btagWeightDown = new TH1D(variable+"_a_"+req+"_btagWeightDown", variable, nBins, customBins);
  sig_a_btagWeightDown->Sumw2();
  h_siga_btagWeightDown.push_back(sig_a_btagWeightDown);

  TH1D * sig_a_topPtUp = new TH1D(variable+"_a_"+req+"_topPtUp", variable, nBins, customBins);
  sig_a_topPtUp->Sumw2();
  h_siga_topPtUp.push_back(sig_a_topPtUp);

  TH1D * sig_a_topPtDown = new TH1D(variable+"_a_"+req+"_topPtDown", variable, nBins, customBins);
  sig_a_topPtDown->Sumw2();
  h_siga_topPtDown.push_back(sig_a_topPtDown);
  
  TH1D * sig_a_JECup = new TH1D(variable+"_a_"+req+"_JECup", variable, nBins, customBins);
  sig_a_JECup->Sumw2();
  h_siga_JECup.push_back(sig_a_JECup);

  TH1D * sig_a_JECdown = new TH1D(variable+"_a_"+req+"_JECdown", variable, nBins, customBins);
  sig_a_JECdown->Sumw2();
  h_siga_JECdown.push_back(sig_a_JECdown);

  TH1D * sig_a_leptonSFup = new TH1D(variable+"_a_"+req+"_leptonSFup", variable, nBins, customBins);
  sig_a_leptonSFup->Sumw2();
  h_siga_leptonSFup.push_back(sig_a_leptonSFup);

  TH1D * sig_a_leptonSFdown = new TH1D(variable+"_a_"+req+"_leptonSFdown", variable, nBins, customBins);
  sig_a_leptonSFdown->Sumw2();
  h_siga_leptonSFdown.push_back(sig_a_leptonSFdown);

  TH1D * sig_a_photonSFup = new TH1D(variable+"_a_"+req+"_photonSFup", variable, nBins, customBins);
  sig_a_photonSFup->Sumw2();
  h_siga_photonSFup.push_back(sig_a_photonSFup);

  TH1D * sig_a_photonSFdown = new TH1D(variable+"_a_"+req+"_photonSFdown", variable, nBins, customBins);
  sig_a_photonSFdown->Sumw2();
  h_siga_photonSFdown.push_back(sig_a_photonSFdown);

  TH1D * sig_b = new TH1D(variable+"_b_"+req, variable, nBins, customBins);
  sig_b->Sumw2();
  h_sigb.push_back(sig_b);

  TH1D * sig_b_btagWeightUp = new TH1D(variable+"_b_"+req+"_btagWeightUp", variable, nBins, customBins);
  sig_b_btagWeightUp->Sumw2();
  h_sigb_btagWeightUp.push_back(sig_b_btagWeightUp);

  TH1D * sig_b_btagWeightDown = new TH1D(variable+"_b_"+req+"_btagWeightDown", variable, nBins, customBins);
  sig_b_btagWeightDown->Sumw2();
  h_sigb_btagWeightDown.push_back(sig_b_btagWeightDown);

  TH1D * sig_b_topPtUp = new TH1D(variable+"_b_"+req+"_topPtUp", variable, nBins, customBins);
  sig_b_topPtUp->Sumw2();
  h_sigb_topPtUp.push_back(sig_b_topPtUp);

  TH1D * sig_b_topPtDown = new TH1D(variable+"_b_"+req+"_topPtDown", variable, nBins, customBins);
  sig_b_topPtDown->Sumw2();
  h_sigb_topPtDown.push_back(sig_b_topPtDown);
  
  TH1D * sig_b_JECup = new TH1D(variable+"_b_"+req+"_JECup", variable, nBins, customBins);
  sig_b_JECup->Sumw2();
  h_sigb_JECup.push_back(sig_b_JECup);

  TH1D * sig_b_JECdown = new TH1D(variable+"_b_"+req+"_JECdown", variable, nBins, customBins);
  sig_b_JECdown->Sumw2();
  h_sigb_JECdown.push_back(sig_b_JECdown);

  TH1D * sig_b_leptonSFup = new TH1D(variable+"_b_"+req+"_leptonSFup", variable, nBins, customBins);
  sig_b_leptonSFup->Sumw2();
  h_sigb_leptonSFup.push_back(sig_b_leptonSFup);

  TH1D * sig_b_leptonSFdown = new TH1D(variable+"_b_"+req+"_leptonSFdown", variable, nBins, customBins);
  sig_b_leptonSFdown->Sumw2();
  h_sigb_leptonSFdown.push_back(sig_b_leptonSFdown);

  TH1D * sig_b_photonSFup = new TH1D(variable+"_b_"+req+"_photonSFup", variable, nBins, customBins);
  sig_b_photonSFup->Sumw2();
  h_sigb_photonSFup.push_back(sig_b_photonSFup);

  TH1D * sig_b_photonSFdown = new TH1D(variable+"_b_"+req+"_photonSFdown", variable, nBins, customBins);
  sig_b_photonSFdown->Sumw2();
  h_sigb_photonSFdown.push_back(sig_b_photonSFdown);
  
}

// expects BookHistogram on nphotons, then met, then others
void PlotMaker::FillHistograms(double metCut, int nPhotons_req, int nBtagReq, int chan) {

  vector<Float_t> vars;
  vars.resize(variables.size());

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr, btagWeightUp, btagWeightDown;
  Float_t overlaps_ttA;
  Float_t topPtReweighting;

  for(unsigned int i = 0; i < variables.size(); i++) {

    ggTree->SetBranchAddress(variables[i], &(vars[i]));
    qcdTree->SetBranchAddress(variables[i], &(vars[i]));

    for(unsigned int j = 0; j < mcTrees.size(); j++) {
      mcTrees[j]->SetBranchAddress(variables[i], &(vars[i]));
      mcTrees_JECup[j]->SetBranchAddress(variables[i], &(vars[i]));
      mcTrees_JECdown[j]->SetBranchAddress(variables[i], &(vars[i]));
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
    
    mcTrees_JECup[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees_JECup[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees_JECdown[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees_JECdown[i]->SetBranchAddress("btagWeight", &btagWeight);

    mcQCDTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcQCDTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcQCDTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);

    if(removeTTAoverlap[i]) {
      mcTrees[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
      mcTrees_JECup[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
      mcTrees_JECdown[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
      mcQCDTrees[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
    }

    if(reweightTopPt[i]) {
      mcTrees[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
      mcTrees_JECup[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
      mcTrees_JECdown[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
      mcQCDTrees[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    }

  }

  sigaTree->SetBranchAddress("pileupWeight", &puWeight);
  sigaTree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigaTree->SetBranchAddress("btagWeight", &btagWeight);
  sigaTree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigaTree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigaTree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigaTree->SetBranchAddress("TopPtReweighting", &topPtReweighting);
  
  sigaTree_JECup->SetBranchAddress("pileupWeight", &puWeight);
  sigaTree_JECup->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigaTree_JECup->SetBranchAddress("btagWeight", &btagWeight);
  sigaTree_JECup->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigaTree_JECup->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigaTree_JECup->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigaTree_JECup->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  sigaTree_JECdown->SetBranchAddress("pileupWeight", &puWeight);
  sigaTree_JECdown->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigaTree_JECdown->SetBranchAddress("btagWeight", &btagWeight);
  sigaTree_JECdown->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigaTree_JECdown->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigaTree_JECdown->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigaTree_JECdown->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  sigbTree->SetBranchAddress("pileupWeight", &puWeight);
  sigbTree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigbTree->SetBranchAddress("btagWeight", &btagWeight);
  sigbTree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigbTree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigbTree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigbTree->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  sigbTree_JECup->SetBranchAddress("pileupWeight", &puWeight);
  sigbTree_JECup->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigbTree_JECup->SetBranchAddress("btagWeight", &btagWeight);
  sigbTree_JECup->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigbTree_JECup->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigbTree_JECup->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigbTree_JECup->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  sigbTree_JECdown->SetBranchAddress("pileupWeight", &puWeight);
  sigbTree_JECdown->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigbTree_JECdown->SetBranchAddress("btagWeight", &btagWeight);
  sigbTree_JECdown->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigbTree_JECdown->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigbTree_JECdown->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigbTree_JECdown->SetBranchAddress("TopPtReweighting", &topPtReweighting);
  
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

  Float_t leptonSF, leptonSFup, leptonSFdown;
  Float_t photonSF, photonSFup, photonSFdown;

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    
    for(int j = 0; j < mcTrees[i]->GetEntries(); j++) {
      mcTrees[i]->GetEntry(j);

      if(nBtagReq == 0) {
	btagWeight = 1.;
	btagWeightErr = 0.;
	btagWeightUp = 1.;
	btagWeightDown = 1.;
      }

      if(useWHIZARD && removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;

      if(btagWeight != btagWeight) continue;
      if(metCut > 0. && vars[1] >= metCut) continue;

      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;

      if(topPtReweighting < 0) topPtReweighting = 1.;

      GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      Float_t addError2_puOnly = btagWeight*btagWeight*puWeightErr*puWeightErr;

      for(unsigned int k = 0; k < vars.size(); k++) {
	if(variables[k] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

	double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;

	Float_t oldError = mcHistograms[i][k]->GetBinError(mcHistograms[i][k]->FindBin(vars[k]));
	Float_t newerror = sqrt(oldError*oldError + addError2);
	mcHistograms[i][k]->Fill(vars[k], totalWeight);
	mcHistograms[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeightUp * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	oldError = mcHistograms_btagWeightUp[i][k]->GetBinError(mcHistograms_btagWeightUp[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2_puOnly);
	mcHistograms_btagWeightUp[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_btagWeightUp[i][k]->SetBinError(mcHistograms_btagWeightUp[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeightDown * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	oldError = mcHistograms_btagWeightDown[i][k]->GetBinError(mcHistograms_btagWeightDown[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2_puOnly);
	mcHistograms_btagWeightDown[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_btagWeightDown[i][k]->SetBinError(mcHistograms_btagWeightDown[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeight * leptonSFup * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	oldError = mcHistograms_leptonSFup[i][k]->GetBinError(mcHistograms_leptonSFup[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_leptonSFup[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_leptonSFup[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeight * leptonSFdown * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	oldError = mcHistograms_leptonSFdown[i][k]->GetBinError(mcHistograms_leptonSFdown[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_leptonSFdown[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_leptonSFdown[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSFup;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	oldError = mcHistograms_photonSFup[i][k]->GetBinError(mcHistograms_photonSFup[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_photonSFup[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_photonSFup[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSFdown;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	oldError = mcHistograms_photonSFdown[i][k]->GetBinError(mcHistograms_photonSFdown[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_photonSFdown[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_photonSFdown[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting * topPtReweighting;
	oldError = mcHistograms_topPtUp[i][k]->GetBinError(mcHistograms_topPtUp[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_topPtUp[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_topPtUp[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSF;
	oldError = mcHistograms_topPtDown[i][k]->GetBinError(mcHistograms_topPtDown[i][k]->FindBin(vars[k]));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_topPtDown[i][k]->Fill(vars[k], totalWeight);
	mcHistograms_topPtDown[i][k]->SetBinError(mcHistograms[i][k]->FindBin(vars[k]), newerror);
	
      }

    }

    for(unsigned int j = 0; j < vars.size(); j++) {
      for(int k = 0; k < mcHistograms[i][j]->GetNbinsX(); k++) {
	Double_t content = mcHistograms[i][j]->GetBinContent(k+1);
	Double_t error = mcHistograms[i][j]->GetBinError(k+1);
      
	mcHistograms_scaleUp[i][j]->SetBinContent(k+1, content);
	mcHistograms_scaleDown[i][j]->SetBinContent(k+1, content);
	mcHistograms_pdfUp[i][j]->SetBinContent(k+1, content);
	mcHistograms_pdfDown[i][j]->SetBinContent(k+1, content);

	mcHistograms_scaleUp[i][j]->SetBinError(k+1, error);
	mcHistograms_scaleDown[i][j]->SetBinError(k+1, error);
	mcHistograms_pdfUp[i][j]->SetBinError(k+1, error);
	mcHistograms_pdfDown[i][j]->SetBinError(k+1, error);
      }
    }

    for(unsigned int j = 0; j < vars.size(); j++) {
      mcHistograms[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_btagWeightUp[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_btagWeightDown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_scaleUp[i][j]->Scale(intLumi_int * (crossSections[i] + scaleErrUp[i]) / mcNGen[i]);
      mcHistograms_scaleDown[i][j]->Scale(intLumi_int * (crossSections[i] - scaleErrDown[i]) / mcNGen[i]);
      mcHistograms_pdfUp[i][j]->Scale(intLumi_int * (crossSections[i] + pdfErrUp[i]) / mcNGen[i]);
      mcHistograms_pdfDown[i][j]->Scale(intLumi_int * (crossSections[i] - pdfErrDown[i]) / mcNGen[i]);
      mcHistograms_topPtUp[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_topPtDown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_leptonSFup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_leptonSFdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_photonSFup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_photonSFdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    }

  }

  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECup[i]->GetEntries(); j++) {
      mcTrees_JECup[i]->GetEntry(j);

      if(nBtagReq == 0) btagWeight = 1.;

      if(useWHIZARD && removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;

      if(btagWeight != btagWeight) continue;
      if(metCut > 0. && vars[1] >= metCut) continue;

      if(topPtReweighting < 0) topPtReweighting = 1.;

      GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;

      for(unsigned int k = 0; k < vars.size(); k++) {
	if(variables[k] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;
	mcHistograms_JECup[i][k]->Fill(vars[k], totalWeight);
      }

    }

    for(unsigned int j = 0; j < vars.size(); j++) mcHistograms_JECup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);

  }

  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECdown[i]->GetEntries(); j++) {
      mcTrees_JECdown[i]->GetEntry(j);

      if(nBtagReq == 0) btagWeight = 1.;

      if(useWHIZARD && removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;

      if(btagWeight != btagWeight) continue;
      if(metCut > 0. && vars[1] >= metCut) continue;

      if(topPtReweighting < 0) topPtReweighting = 1.;

      GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;

      for(unsigned int k = 0; k < vars.size(); k++) {
	if(variables[k] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;
	mcHistograms_JECdown[i][k]->Fill(vars[k], totalWeight);
      }

    }

    for(unsigned int j = 0; j < vars.size(); j++) mcHistograms_JECdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);

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

      if(useWHIZARD && i < 3 && overlaps_ttA > 0.001) continue;

      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;

      GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;

      for(unsigned int k = 0; k < vars.size(); k++) {
	if(variables[k] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

	Float_t oldError = mcQCDHistograms[i][k]->GetBinError(mcQCDHistograms[i][k]->FindBin(vars[k]));
	Float_t newerror = sqrt(oldError*oldError + addError2);
	mcQCDHistograms[i][k]->Fill(vars[k], totalWeight);
	mcQCDHistograms[i][k]->SetBinError(mcQCDHistograms[i][k]->FindBin(vars[k]), newerror);
      }

    }

    for(unsigned int j = 0; j < vars.size(); j++) mcQCDHistograms[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
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

    if(topPtReweighting < 0) topPtReweighting = 1.;

    GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h_siga[j]->GetBinError(h_siga[j]->FindBin(vars[j]));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_siga[j]->Fill(vars[j], totalWeight);
      h_siga[j]->SetBinError(h_siga[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeightUp * leptonSF * photonSF * topPtReweighting;
      olderror = h_siga_btagWeightUp[j]->GetBinError(h_siga_btagWeightUp[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_btagWeightUp[j]->Fill(vars[j], totalWeight);
      h_siga_btagWeightUp[j]->SetBinError(h_siga_btagWeightUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeightDown * leptonSF * photonSF * topPtReweighting;
      olderror = h_siga_btagWeightDown[j]->GetBinError(h_siga_btagWeightDown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_btagWeightDown[j]->Fill(vars[j], totalWeight);
      h_siga_btagWeightDown[j]->SetBinError(h_siga_btagWeightDown[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting * topPtReweighting;
      olderror = h_siga_topPtUp[j]->GetBinError(h_siga_topPtUp[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_topPtUp[j]->Fill(vars[j], totalWeight);
      h_siga_topPtUp[j]->SetBinError(h_siga_topPtUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      olderror = h_siga_topPtDown[j]->GetBinError(h_siga_topPtDown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_topPtDown[j]->Fill(vars[j], totalWeight);
      h_siga_topPtDown[j]->SetBinError(h_siga_topPtDown[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSFup * photonSF * topPtReweighting;
      olderror = h_siga_leptonSFup[j]->GetBinError(h_siga_leptonSFup[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_leptonSFup[j]->Fill(vars[j], totalWeight);
      h_siga_leptonSFup[j]->SetBinError(h_siga_btagWeightUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSFdown * photonSF * topPtReweighting;
      olderror = h_siga_leptonSFdown[j]->GetBinError(h_siga_leptonSFdown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_leptonSFdown[j]->Fill(vars[j], totalWeight);
      h_siga_leptonSFdown[j]->SetBinError(h_siga_leptonSFdown[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFup * topPtReweighting;
      olderror = h_siga_photonSFup[j]->GetBinError(h_siga_photonSFup[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_photonSFup[j]->Fill(vars[j], totalWeight);
      h_siga_photonSFup[j]->SetBinError(h_siga_btagWeightUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFdown * topPtReweighting;
      olderror = h_siga_photonSFdown[j]->GetBinError(h_siga_photonSFdown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_siga_photonSFdown[j]->Fill(vars[j], totalWeight);
      h_siga_photonSFdown[j]->SetBinError(h_siga_photonSFdown[j]->FindBin(vars[j]), newerror);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) {
    h_siga[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_btagWeightUp[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_btagWeightDown[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_topPtUp[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_topPtDown[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_leptonSFup[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_leptonSFdown[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_photonSFup[j]->Scale(intLumi_int * 0.147492 / 15000.);
    h_siga_photonSFdown[j]->Scale(intLumi_int * 0.147492 / 15000.);
  }

  for(int i = 0; i < sigaTree_JECup->GetEntries(); i++) {
    sigaTree_JECup->GetEntry(i);

    if(nBtagReq == 0) btagWeight = 1.;

    if(btagWeight != btagWeight) continue;
    if(metCut > 0. && vars[1] >= metCut) continue;

    if(topPtReweighting < 0) topPtReweighting = 1.;

    GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      h_siga_JECup[j]->Fill(vars[j], totalWeight);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) h_siga_JECup[j]->Scale(intLumi_int * 0.147492 / 15000.);

  for(int i = 0; i < sigaTree_JECdown->GetEntries(); i++) {
    sigaTree_JECdown->GetEntry(i);

    if(nBtagReq == 0) btagWeight = 1.;

    if(btagWeight != btagWeight) continue;
    if(metCut > 0. && vars[1] >= metCut) continue;

    if(topPtReweighting < 0) topPtReweighting = 1.;

    GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      h_siga_JECdown[j]->Fill(vars[j], totalWeight);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) h_siga_JECdown[j]->Scale(intLumi_int * 0.147492 / 15000.);

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

    if(topPtReweighting < 0) topPtReweighting = 1.;

    GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h_sigb[j]->GetBinError(h_sigb[j]->FindBin(vars[j]));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_sigb[j]->Fill(vars[j], totalWeight);
      h_sigb[j]->SetBinError(h_sigb[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeightUp * leptonSF * photonSF * topPtReweighting;
      olderror = h_sigb_btagWeightUp[j]->GetBinError(h_sigb_btagWeightUp[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_btagWeightUp[j]->Fill(vars[j], totalWeight);
      h_sigb_btagWeightUp[j]->SetBinError(h_sigb_btagWeightUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeightDown * leptonSF * photonSF * topPtReweighting;
      olderror = h_sigb_btagWeightDown[j]->GetBinError(h_sigb_btagWeightDown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_btagWeightDown[j]->Fill(vars[j], totalWeight);
      h_sigb_btagWeightDown[j]->SetBinError(h_sigb_btagWeightDown[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting * topPtReweighting;
      olderror = h_sigb_topPtUp[j]->GetBinError(h_sigb_topPtUp[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_topPtUp[j]->Fill(vars[j], totalWeight);
      h_sigb_topPtUp[j]->SetBinError(h_sigb_topPtUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      olderror = h_sigb_topPtDown[j]->GetBinError(h_sigb_topPtDown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_topPtDown[j]->Fill(vars[j], totalWeight);
      h_sigb_topPtDown[j]->SetBinError(h_sigb_topPtDown[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSFup * photonSF * topPtReweighting;
      olderror = h_sigb_leptonSFup[j]->GetBinError(h_sigb_leptonSFup[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_leptonSFup[j]->Fill(vars[j], totalWeight);
      h_sigb_leptonSFup[j]->SetBinError(h_sigb_btagWeightUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSFdown * photonSF * topPtReweighting;
      olderror = h_sigb_leptonSFdown[j]->GetBinError(h_sigb_leptonSFdown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_leptonSFdown[j]->Fill(vars[j], totalWeight);
      h_sigb_leptonSFdown[j]->SetBinError(h_sigb_leptonSFdown[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFup * topPtReweighting;
      olderror = h_sigb_photonSFup[j]->GetBinError(h_sigb_photonSFup[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_photonSFup[j]->Fill(vars[j], totalWeight);
      h_sigb_photonSFup[j]->SetBinError(h_sigb_btagWeightUp[j]->FindBin(vars[j]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFdown * topPtReweighting;
      olderror = h_sigb_photonSFdown[j]->GetBinError(h_sigb_photonSFdown[j]->FindBin(vars[j]));
      newerror = sqrt(olderror*olderror + addError2);
      h_sigb_photonSFdown[j]->Fill(vars[j], totalWeight);
      h_sigb_photonSFdown[j]->SetBinError(h_sigb_photonSFdown[j]->FindBin(vars[j]), newerror);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) {
    h_sigb[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_btagWeightUp[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_btagWeightDown[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_topPtUp[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_topPtDown[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_leptonSFup[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_leptonSFdown[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_photonSFup[j]->Scale(intLumi_int * 0.0399591 / 15000.);
    h_sigb_photonSFdown[j]->Scale(intLumi_int * 0.0399591 / 15000.);
  }

  for(int i = 0; i < sigbTree_JECup->GetEntries(); i++) {
    sigbTree_JECup->GetEntry(i);

    if(nBtagReq == 0) btagWeight = 1.;

    if(btagWeight != btagWeight) continue;
    if(metCut > 0. && vars[1] >= metCut) continue;

    if(topPtReweighting < 0) topPtReweighting = 1.;

    GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      h_sigb_JECup[j]->Fill(vars[j], totalWeight);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) h_sigb_JECup[j]->Scale(intLumi_int * 0.0399591 / 15000.);

  for(int i = 0; i < sigbTree_JECdown->GetEntries(); i++) {
    sigbTree_JECdown->GetEntry(i);

    if(nBtagReq == 0) btagWeight = 1.;

    if(btagWeight != btagWeight) continue;
    if(metCut > 0. && vars[1] >= metCut) continue;

    if(topPtReweighting < 0) topPtReweighting = 1.;

    GetLeptonSF(vars, chan, leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(vars, photonSF, photonSFup, photonSFdown);

    for(unsigned int j = 0; j < vars.size(); j++) {
      if(variables[j] != "Nphotons" && (int)vars[0] != nPhotons_req) continue;

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      h_sigb_JECdown[j]->Fill(vars[j], totalWeight);
    }

  }
  for(unsigned int j = 0; j < vars.size(); j++) h_sigb_JECdown[j]->Scale(intLumi_int * 0.0399591 / 15000.);

  ggTree->ResetBranchAddresses();
  qcdTree->ResetBranchAddresses();

  for(unsigned int i = 0; i < mcTrees.size(); i++) mcTrees[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) mcTrees_JECup[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) mcTrees_JECdown[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) mcQCDTrees[i]->ResetBranchAddresses();

  sigaTree->ResetBranchAddresses();
  sigaTree_JECup->ResetBranchAddresses();
  sigaTree_JECdown->ResetBranchAddresses();

  sigbTree->ResetBranchAddresses();
  sigbTree_JECup->ResetBranchAddresses();
  sigbTree_JECdown->ResetBranchAddresses();

}

void PlotMaker::SubtractMCFromQCD() {

  for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms[i].size(); j++) {
      h_qcd[j]->Add(mcQCDHistograms[i][j], -1.);
     }
  }

  for(unsigned int i = 0; i < h_qcd.size(); i++) {
    for(Int_t j = 0; j < h_qcd[i]->GetNbinsX(); j++) {
      if(h_qcd[i]->GetBinContent(j+1) < 0) {
	h_qcd[i]->SetBinContent(j+1, 0.);
	h_qcd[i]->SetBinError(j+1, 0.);
      }
    }
  }

}

void PlotMaker::NormalizeQCD() {

  double n_qcd_before = h_qcd[1]->Integral();

  const int endBin = h_gg[1]->GetXaxis()->FindBin(20) - 1;

  double n_sig = h_gg[1]->Integral(0, endBin);
  double n_qcd = h_qcd[1]->Integral(0, endBin);

  if(n_qcd < 1) return;

  double n_mc = 0;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) n_mc += mcHistograms[i][1]->Integral(0, mcHistograms[i][1]->GetXaxis()->FindBin(20) - 1);

  double sigma_sig = 0;
  double sigma_qcd = 0;
  double sigma_mc = 0;
  
  for(int i = 0; i < endBin; i++) {
    sigma_sig += h_gg[1]->GetBinError(i+1) * h_gg[1]->GetBinError(i+1);
    sigma_qcd += h_qcd[1]->GetBinError(i+1) * h_qcd[1]->GetBinError(i+1);

    for(unsigned int j = 0; j < mcHistograms.size(); j++) sigma_mc += mcHistograms[j][1]->GetBinError(i+1) * mcHistograms[j][1]->GetBinError(i+1);
  }

  sigma_sig = sqrt(sigma_sig);
  sigma_qcd = sqrt(sigma_qcd);
  sigma_mc = sqrt(sigma_mc);

  double scale = (n_sig - n_mc) / n_qcd;
  
  for(unsigned int i = 0; i < h_qcd.size(); i++) {
    double newError_lowBins[endBin];

    if(i == 1) {
      for(int j = 0; j < endBin; j++) {
	newError_lowBins[j] = h_qcd[1]->GetBinContent(j+1) * h_qcd[1]->GetBinContent(j+1);
	newError_lowBins[j] *= sigma_sig*sigma_sig + sigma_mc*sigma_mc;
	newError_lowBins[j] += scale*scale * h_qcd[1]->GetBinError(j+1)*h_qcd[1]->GetBinError(j+1) * (n_qcd - h_qcd[1]->GetBinContent(j+1))*(n_qcd - h_qcd[1]->GetBinContent(j+1));
	for(int k = 0; k < endBin; k++) {
	  if(k == j) continue;
	  newError_lowBins[j] += scale*scale * h_qcd[1]->GetBinError(k+1)*h_qcd[1]->GetBinError(k+1) * h_qcd[1]->GetBinContent(k+1)*h_qcd[1]->GetBinContent(k+1);
	}
	newError_lowBins[j] /= n_qcd*n_qcd;
	newError_lowBins[j] = sqrt(newError_lowBins[j]);
      }
    }

    for(int j = 0; j < h_qcd[i]->GetNbinsX(); j++) {
      // do something different for pfMET < 20
      if(j < endBin+1 && i == 1) {
	h_qcd[i]->SetBinContent(j+1, h_qcd[i]->GetBinContent(j+1) * scale);
	h_qcd[i]->SetBinError(j+1, newError_lowBins[j]);
      }
      else {
	double newError = sigma_sig*sigma_sig + sigma_mc*sigma_mc + scale*scale*sigma_qcd*sigma_qcd;
	newError *= h_qcd[i]->GetBinContent(j+1)*h_qcd[i]->GetBinContent(j+1) / (n_qcd*n_qcd);
	newError += h_qcd[i]->GetBinError(j+1)*h_qcd[i]->GetBinError(j+1) * scale*scale;
	newError = sqrt(newError);

	h_qcd[i]->SetBinContent(j+1, h_qcd[i]->GetBinContent(j+1) * scale);
	h_qcd[i]->SetBinError(j+1, newError);
      }
    }

  }

  double n_qcd_after = h_qcd[1]->Integral();

  cout << endl << "NormalizeQCD(): Overall scaling factor of " << n_qcd_after / n_qcd_before << " applied." << endl << endl;

}
    
void PlotMaker::CreatePlot(TString variable,
			   bool divideByWidth, bool needsQCD,
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
    if(needsQCD) h_qcd[var_num] = (TH1D*)DivideByBinWidth(h_qcd[var_num]);

    for(unsigned int i = 0; i < mcHistograms.size(); i++) mcHistograms[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_btagWeightUp.size(); i++) mcHistograms_btagWeightUp[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_btagWeightUp[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_btagWeightDown.size(); i++) mcHistograms_btagWeightDown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_btagWeightDown[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_scaleUp.size(); i++) mcHistograms_scaleUp[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_scaleUp[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_scaleDown.size(); i++) mcHistograms_scaleDown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_scaleDown[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_pdfUp.size(); i++) mcHistograms_pdfUp[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_pdfUp[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_pdfDown.size(); i++) mcHistograms_pdfDown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_pdfDown[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_topPtUp.size(); i++) mcHistograms_topPtUp[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_topPtUp[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_topPtDown.size(); i++) mcHistograms_topPtDown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_topPtDown[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_JECup.size(); i++) mcHistograms_JECup[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_JECup[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_JECdown.size(); i++) mcHistograms_JECdown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_JECdown[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_leptonSFup.size(); i++) mcHistograms_leptonSFup[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_leptonSFup[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_leptonSFdown.size(); i++) mcHistograms_leptonSFdown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_leptonSFdown[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_photonSFup.size(); i++) mcHistograms_photonSFup[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_photonSFup[i][var_num]);
    for(unsigned int i = 0; i < mcHistograms_photonSFdown.size(); i++) mcHistograms_photonSFdown[i][var_num] = (TH1D*)DivideByBinWidth(mcHistograms_photonSFdown[i][var_num]);

    for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) mcQCDHistograms[i][var_num] = (TH1D*)DivideByBinWidth(mcQCDHistograms[i][var_num]);

    h_siga[var_num] = (TH1D*)DivideByBinWidth(h_siga[var_num]);
    h_sigb[var_num] = (TH1D*)DivideByBinWidth(h_sigb[var_num]);
  }

  DrawPlot(var_num, variable, needsQCD,
	   xaxisTitle, yaxisTitle,
	   xmin, xmax,
	   ymin, ymax,
	   ratiomin, ratiomax,
	   drawSignal, drawLegend, drawPrelim,
	   out);

}

void PlotMaker::DrawPlot(int variableNumber, TString variable, bool needsQCD,
			 TString xaxisTitle, TString yaxisTitle,
			 Float_t xmin, Float_t xmax,
			 Float_t ymin, Float_t ymax,
			 Float_t ratiomin, Float_t ratiomax,
			 bool drawSignal, bool drawLegend, bool drawPrelim,
			 TFile*& out) {

  out->cd();
  h_gg[variableNumber]->Write();
  if(needsQCD) h_qcd[variableNumber]->Write();
  for(unsigned int i = 0; i < mcHistograms.size(); i++) mcHistograms[i][variableNumber]->Write();
  for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) mcQCDHistograms[i][variableNumber]->Write();

  TH1D *bkg, 
    *bkg_btagWeightUp, *bkg_btagWeightDown, 
    *bkg_scaleUp, *bkg_scaleDown, 
    *bkg_pdfUp, *bkg_pdfDown, 
    *bkg_topPtUp, *bkg_topPtDown, 
    *bkg_JECup, *bkg_JECdown,
    *bkg_leptonSFup, *bkg_leptonSFdown,
    *bkg_photonSFup, *bkg_photonSFdown;

  // Stack histograms; qcd is on top
  if(needsQCD) {
    bkg = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req);
    bkg_btagWeightUp = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_btagWeightUp");
    bkg_btagWeightDown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_btagWeightDown");
    bkg_scaleUp = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_scaleUp");
    bkg_scaleDown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_scaleDown");
    bkg_pdfUp = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_pdfUp");
    bkg_pdfDown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_pdfDown");
    bkg_topPtUp = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_topPtUp");
    bkg_topPtDown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_topPtDown");
    bkg_JECup = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_JECup");
    bkg_JECdown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_JECdown");
    bkg_leptonSFup = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_leptonSFup");
    bkg_leptonSFdown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_leptonSFdown");
    bkg_photonSFup = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_photonSFup");
    bkg_photonSFdown = (TH1D*)h_qcd[variableNumber]->Clone(variable+"_bkg_"+req+"_photonSFdown");

    for(unsigned int i = 0; i < mcHistograms.size(); i++) {
      bkg->Add(mcHistograms[i][variableNumber]);
      bkg_btagWeightUp->Add(mcHistograms_btagWeightUp[i][variableNumber]);
      bkg_btagWeightDown->Add(mcHistograms_btagWeightDown[i][variableNumber]);
      bkg_scaleUp->Add(mcHistograms_scaleUp[i][variableNumber]);
      bkg_scaleDown->Add(mcHistograms_scaleDown[i][variableNumber]);
      bkg_pdfUp->Add(mcHistograms_pdfUp[i][variableNumber]);
      bkg_pdfDown->Add(mcHistograms_pdfDown[i][variableNumber]);
      bkg_topPtUp->Add(mcHistograms_topPtUp[i][variableNumber]);
      bkg_topPtDown->Add(mcHistograms_topPtDown[i][variableNumber]);
      bkg_JECup->Add(mcHistograms_JECup[i][variableNumber]);
      bkg_JECdown->Add(mcHistograms_JECdown[i][variableNumber]);
      bkg_leptonSFup->Add(mcHistograms_leptonSFup[i][variableNumber]);
      bkg_leptonSFdown->Add(mcHistograms_leptonSFdown[i][variableNumber]);
      bkg_photonSFup->Add(mcHistograms_photonSFup[i][variableNumber]);
      bkg_photonSFdown->Add(mcHistograms_photonSFdown[i][variableNumber]);

      for(unsigned int j = i + 1; j < mcHistograms.size(); j++) {
	mcHistograms[i][variableNumber]->Add(mcHistograms[j][variableNumber]);
      }
    }
  }
  else {
    bkg = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req);
    bkg_btagWeightUp = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_btagWeightUp");
    bkg_btagWeightDown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_btagWeightDown");
    bkg_scaleUp = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_scaleUp");
    bkg_scaleDown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_scaleDown");
    bkg_pdfUp = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_pdfUp");
    bkg_pdfDown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_pdfDown");
    bkg_topPtUp = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_topPtUp");
    bkg_topPtDown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_topPtDown");
    bkg_JECup = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_JECup");
    bkg_JECdown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_JECdown");
    bkg_leptonSFup = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_leptonSFup");
    bkg_leptonSFdown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_leptonSFdown");
    bkg_photonSFup = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_photonSFup");
    bkg_photonSFdown = (TH1D*)mcHistograms[0][variableNumber]->Clone(variable+"_bkg_"+req+"_photonSFdown");

    for(unsigned int i = 1; i < mcHistograms.size(); i++) {
      bkg->Add(mcHistograms[i][variableNumber]);
      bkg_btagWeightUp->Add(mcHistograms_btagWeightUp[i][variableNumber]);
      bkg_btagWeightDown->Add(mcHistograms_btagWeightDown[i][variableNumber]);
      bkg_scaleUp->Add(mcHistograms_scaleUp[i][variableNumber]);
      bkg_scaleDown->Add(mcHistograms_scaleDown[i][variableNumber]);
      bkg_pdfUp->Add(mcHistograms_pdfUp[i][variableNumber]);
      bkg_pdfDown->Add(mcHistograms_pdfDown[i][variableNumber]);
      bkg_topPtUp->Add(mcHistograms_topPtUp[i][variableNumber]);
      bkg_topPtDown->Add(mcHistograms_topPtDown[i][variableNumber]);
      bkg_JECup->Add(mcHistograms_JECup[i][variableNumber]);
      bkg_JECdown->Add(mcHistograms_JECdown[i][variableNumber]);
      bkg_leptonSFup->Add(mcHistograms_leptonSFup[i][variableNumber]);
      bkg_leptonSFdown->Add(mcHistograms_leptonSFdown[i][variableNumber]);
      bkg_photonSFup->Add(mcHistograms_photonSFup[i][variableNumber]);
      bkg_photonSFdown->Add(mcHistograms_photonSFdown[i][variableNumber]);

      for(unsigned int j = i + 1; j < mcHistograms.size(); j++) {
	mcHistograms[i][variableNumber]->Add(mcHistograms[j][variableNumber]);
      }
    }
  }
   
  bkg->Write();
  bkg_btagWeightUp->Write();
  bkg_btagWeightDown->Write();
  bkg_scaleUp->Write();
  bkg_scaleDown->Write();
  bkg_pdfUp->Write();
  bkg_pdfDown->Write();
  bkg_topPtUp->Write();
  bkg_topPtDown->Write();
  bkg_JECup->Write();
  bkg_JECdown->Write();
  bkg_leptonSFup->Write();
  bkg_leptonSFdown->Write();
  bkg_photonSFup->Write();
  bkg_photonSFdown->Write();

  Double_t kolm = h_gg[variableNumber]->KolmogorovTest(bkg);
  TString kolmText = Form("KS test probability = %5.3g", kolm);
  TText * tt = new TText(0.92, 0.5, kolmText);
  tt->SetTextAngle(90.);
  tt->SetNDC(); tt->SetTextSize( 0.032 );

  TH1D * errors_stat = (TH1D*)bkg->Clone("errors_stat");

  TH1D * errors_sys = (TH1D*)bkg->Clone("errors_sys");
  for(int i = 0; i < errors_sys->GetNbinsX(); i++) {

    Double_t stat = bkg->GetBinError(i+1);

    Double_t btagUp = bkg_btagWeightUp->GetBinContent(i+1);
    Double_t btagDown = bkg_btagWeightDown->GetBinContent(i+1);
    Double_t btag_sys = fabs(btagUp - btagDown) / 2.;

    Double_t scaleUp = bkg_scaleUp->GetBinContent(i+1);
    Double_t scaleDown = bkg_scaleDown->GetBinContent(i+1);
    Double_t scale_sys = fabs(scaleUp - scaleDown) / 2.;
    
    Double_t pdfUp = bkg_pdfUp->GetBinContent(i+1);
    Double_t pdfDown = bkg_pdfDown->GetBinContent(i+1);
    Double_t pdf_sys = fabs(pdfUp - pdfDown) / 2.;

    Double_t topPtUp = bkg_topPtUp->GetBinContent(i+1);
    Double_t topPtDown = bkg_topPtDown->GetBinContent(i+1);
    Double_t topPt_sys = fabs(topPtUp - topPtDown) / 2.;

    Double_t JECup = bkg_JECup->GetBinContent(i+1);
    Double_t JECdown = bkg_JECdown->GetBinContent(i+1);
    Double_t JEC_sys = fabs(JECup - JECdown) / 2.;

    Double_t leptonSFup = bkg_leptonSFup->GetBinContent(i+1);
    Double_t leptonSFdown = bkg_leptonSFdown->GetBinContent(i+1);
    Double_t leptonSF_sys = fabs(leptonSFup - leptonSFdown) / 2.;

    Double_t photonSFup = bkg_photonSFup->GetBinContent(i+1);
    Double_t photonSFdown = bkg_photonSFdown->GetBinContent(i+1);
    Double_t photonSF_sys = fabs(photonSFup - photonSFdown) / 2.;

    Double_t totalError2 = stat*stat + 
      btag_sys*btag_sys + 
      scale_sys*scale_sys + 
      pdf_sys*pdf_sys + 
      topPt_sys*topPt_sys + 
      JEC_sys*JEC_sys + 
      leptonSF_sys*leptonSF_sys + 
      photonSF_sys*photonSF_sys;

    if(bkg->GetBinContent(i+1) == 0.) errors_sys->SetBinError(i+1, 0.);
    else errors_sys->SetBinError(i+1, sqrt(totalError2));
  }

  if(drawSignal) calculateROC(h_siga[variableNumber], h_sigb[variableNumber], bkg, req, variable);

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, NULL, "brNDC");
  leg->SetNColumns(2);
  leg->AddEntry(h_gg[variableNumber], "Data", "LP");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(errors_sys, "Stat. #oplus Syst. Errors", "F");
  leg->AddEntry((TObject*)0, "", "");
  if(needsQCD) leg->AddEntry(bkg, "QCD", "F");
  
  leg->AddEntry(mcHistograms[0][variableNumber], legendNames[0], "F");
  for(unsigned int i = 1; i < mcHistograms.size(); i++) {
    if(mcLayerNumbers[i] != mcLayerNumbers[i-1]) leg->AddEntry(mcHistograms[i][variableNumber], legendNames[i], "F");
  }
  if(needsQCD) leg->AddEntry((TObject*)0, "", "");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  TPaveText * reqText = new TPaveText(0.45, 0.47, 0.85, 0.57, "NDC");
  reqText->SetFillColor(0);
  reqText->SetFillStyle(0);
  reqText->SetLineColor(0);
  reqText->AddText(req+" Requirement");

  TPaveText * lumiHeader = new TPaveText(0.1, 0.901, 0.9, 0.94, "NDC");
  lumiHeader->SetFillColor(0);
  lumiHeader->SetFillStyle(0);
  lumiHeader->SetLineColor(0);
  lumiHeader->AddText("CMS Preliminary 2014     #sqrt{s} = 8 TeV     #intL = "+intLumi+" fb^{-1}");
  
  h_gg[variableNumber]->SetMarkerStyle(20); 
  h_gg[variableNumber]->SetMarkerSize(1.5);

  errors_stat->SetFillColor(kOrange+10);
  errors_stat->SetFillStyle(3154);
  errors_stat->SetMarkerSize(0);

  errors_sys->SetFillColor(kOrange+10);
  errors_sys->SetFillStyle(3154);
  errors_sys->SetMarkerSize(0);

  if(needsQCD) bkg->SetFillColor(kSpring-6);
  else bkg->SetFillColor(mcLayerColors[0]);
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
  if(needsQCD) mcHistograms[0][variableNumber]->Draw("same hist");
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    if(i != 0 && mcLayerNumbers[i] != mcLayerNumbers[i-1]) mcHistograms[i][variableNumber]->Draw("same hist");
  }
  //errors_stat->Draw("same e2");
  errors_sys->Draw("same e2");
  h_gg[variableNumber]->Draw("same e1");
  bkg->Draw("same axis");

  if(drawSignal) {
    h_siga[variableNumber]->SetLineColor(kMagenta);
    h_siga[variableNumber]->SetLineWidth(3);
    leg->AddEntry(h_siga[variableNumber], "GGM (460_175)", "L");
    h_siga[variableNumber]->Draw("same hist");
    
    h_sigb[variableNumber]->SetLineColor(kBlue);
    h_sigb[variableNumber]->SetLineWidth(3);
    leg->AddEntry(h_sigb[variableNumber], "GGM (560_325)", "L");
    h_sigb[variableNumber]->Draw("same hist");
  }

  lumiHeader->Draw("same");
  if(drawLegend) leg->Draw("same");
  if(drawPrelim && drawLegend) reqText->Draw("same");
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

  TH1D * ratio_stat = (TH1D*)bkg->Clone("ratio_stat");
  for(int i = 0; i < ratio_stat->GetNbinsX(); i++) {
    ratio_stat->SetBinContent(i+1, 1.);
    if(bkg->GetBinContent(i+1) == 0.) ratio_stat->SetBinError(i+1, 0.);
    else ratio_stat->SetBinError(i+1, bkg->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  TH1D * ratio_sys = (TH1D*)bkg->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {

    Double_t stat = bkg->GetBinError(i+1);

    Double_t btagUp = bkg_btagWeightUp->GetBinContent(i+1);
    Double_t btagDown = bkg_btagWeightDown->GetBinContent(i+1);
    Double_t btag_sys = fabs(btagUp - btagDown) / 2.;

    Double_t scaleUp = bkg_scaleUp->GetBinContent(i+1);
    Double_t scaleDown = bkg_scaleDown->GetBinContent(i+1);
    Double_t scale_sys = fabs(scaleUp - scaleDown) / 2.;
    
    Double_t pdfUp = bkg_pdfUp->GetBinContent(i+1);
    Double_t pdfDown = bkg_pdfDown->GetBinContent(i+1);
    Double_t pdf_sys = fabs(pdfUp - pdfDown) / 2.;

    Double_t topPtUp = bkg_topPtUp->GetBinContent(i+1);
    Double_t topPtDown = bkg_topPtDown->GetBinContent(i+1);
    Double_t topPt_sys = fabs(topPtUp - topPtDown) / 2.;

    Double_t JECup = bkg_JECup->GetBinContent(i+1);
    Double_t JECdown = bkg_JECdown->GetBinContent(i+1);
    Double_t JEC_sys = fabs(JECup - JECdown) / 2.;

    Double_t totalError2 = stat*stat + btag_sys*btag_sys + scale_sys*scale_sys + pdf_sys*pdf_sys + topPt_sys*topPt_sys + JEC_sys*JEC_sys;

    ratio_sys->SetBinContent(i+1, 1.);

    if(bkg->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, sqrt(totalError2) / bkg->GetBinContent(i+1));
  }

  if(xmax > xmin) ratio->GetXaxis()->SetRangeUser(xmin, xmax);

  ratio_stat->SetFillStyle(1001);
  ratio_stat->SetFillColor(kGray+1);
  ratio_stat->SetLineColor(kGray+1);
  ratio_stat->SetMarkerColor(kGray+1);

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

  TLegend * leg2 = new TLegend(0.78, 0.7, 0.88, 0.95, NULL, "brNDC");
  //leg2->SetNColumns(2);
  leg2->AddEntry(ratio_stat, "Stat.", "F");
  leg2->AddEntry(ratio_sys, "Stat. #oplus Syst.", "F");
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.032);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio_stat->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  leg2->Draw("same");

  TLine * oneLine = new TLine(xmin, 1, xmax, 1);
  oneLine->SetLineStyle(2);
  oneLine->Draw();  

  padhi->cd();
  padhi->SetLogy(true);
  can->SaveAs(variable+"_"+req+".pdf");

  delete can;

}

void PlotMaker::CreateTable() {

  // pfMET
  int variableNumber = 1;

  const int nBins = 5;
  Double_t rangeLow[nBins] = {0, 0, 50, 80, 100};
  Double_t rangeHigh[nBins] = {20, 50, -1, -1, -1};
  
  FILE * tableFile = fopen("errorTable_"+req+".temp", "w");
  FILE * datacardFile = fopen("datacard_"+req+".temp", "w");

  Double_t binLow[nBins], binHigh[nBins];
  for(int i = 0; i < nBins; i++) {
    binLow[i] = h_gg[variableNumber]->GetXaxis()->FindBin(rangeLow[i]);
    binHigh[i] = (rangeHigh[i] == -1) ? -1 : h_gg[variableNumber]->GetXaxis()->FindBin(rangeHigh[i]) - 1;
  }
  
  for(int i = 0; i < nBins; i++) {
    
    Double_t this_val, this_err, this_staterr2, this_syserr2_up, this_syserr2_down;
    Double_t this_btagUp, this_scaleUp, this_pdfUp, this_topPtUp, this_JECup, this_leptonSFup, this_photonSFup;
    Double_t this_btagDown, this_scaleDown, this_pdfDown, this_topPtDown, this_JECdown, this_leptonSFdown, this_photonSFdown;

    Double_t bkgval = 0;
    Double_t bkgstat2 = 0;
    Double_t bkgsys2_up = 0;
    Double_t bkgsys2_down = 0;

    Double_t bkg_btagUp2 = 0;
    Double_t bkg_scaleUp2 = 0;
    Double_t bkg_pdfUp2 = 0;
    Double_t bkg_topPtUp2 = 0;
    Double_t bkg_JECup2 = 0;
    Double_t bkg_leptonSFup2 = 0;
    Double_t bkg_photonSFup2 = 0;

    Double_t bkg_btagDown2 = 0;
    Double_t bkg_scaleDown2 = 0;
    Double_t bkg_pdfDown2 = 0;
    Double_t bkg_topPtDown2 = 0;
    Double_t bkg_JECdown2 = 0;
    Double_t bkg_leptonSFdown2 = 0;
    Double_t bkg_photonSFdown2 = 0;

    if(req.Contains("ele")) {
      this_val = h_qcd[variableNumber]->IntegralAndError(binLow[i], binHigh[i], this_err);
      bkgval += this_val;
      bkgstat2 += this_err*this_err;
      fprintf(tableFile, "qcdval%dx:%.1f\nqcdstat%dx:%.2f\n", i+1, this_val, i+1, this_err);
    }
    else fprintf(tableFile, "qcdval%dx:%.1f\nqcdstat%dx:%.2f\n", i+1, 0., i+1, 0.);
    
    for(unsigned int j = 0; j < mcHistograms.size(); j++) {

      if(j > 0 && tableNames[j] == tableNames[j-1]) continue;

      this_val = 0;
      this_staterr2 = 0;
      
      this_btagUp = 0;
      this_scaleUp = 0;
      this_pdfUp = 0;
      this_topPtUp = 0;
      this_JECup = 0;
      this_leptonSFup = 0;
      this_photonSFup = 0;

      this_btagDown = 0;
      this_scaleDown = 0;
      this_pdfDown = 0;
      this_topPtDown = 0;
      this_JECdown = 0;
      this_leptonSFdown = 0;
      this_photonSFdown = 0;

      for(unsigned int k = j; k < mcHistograms.size() && tableNames[k] == tableNames[j]; k++) {

	this_val += mcHistograms[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], this_err);
	this_staterr2 += this_err*this_err;

	Double_t temperr;
	this_btagUp += mcHistograms_btagWeightUp[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_scaleUp += mcHistograms_scaleUp[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_pdfUp += mcHistograms_pdfUp[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_topPtUp += mcHistograms_topPtUp[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_JECup += mcHistograms_JECup[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_leptonSFup += mcHistograms_leptonSFup[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_photonSFup += mcHistograms_photonSFup[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);

	this_btagDown += mcHistograms_btagWeightDown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_scaleDown += mcHistograms_scaleDown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_pdfDown += mcHistograms_pdfDown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_topPtDown += mcHistograms_topPtDown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_JECdown += mcHistograms_JECdown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_leptonSFdown += mcHistograms_leptonSFdown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
	this_photonSFdown += mcHistograms_photonSFdown[k][variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);

      }

      bkgval += this_val;
      bkgstat2 += this_staterr2;

      bkg_btagUp2 += (this_btagUp - this_val)*(this_btagUp - this_val);
      bkg_scaleUp2 += (this_scaleUp - this_val)*(this_scaleUp - this_val);
      bkg_pdfUp2 += (this_pdfUp - this_val)*(this_pdfUp - this_val);
      bkg_topPtUp2 += (this_topPtUp - this_val)*(this_topPtUp - this_val);
      bkg_JECup2 += (this_JECup - this_val)*(this_JECup - this_val);
      bkg_leptonSFup2 += (this_leptonSFup - this_val)*(this_leptonSFup - this_val);
      bkg_photonSFup2 += (this_photonSFup - this_val)*(this_photonSFup - this_val);

      bkg_btagDown2 += (this_btagDown - this_val)*(this_btagDown - this_val);
      bkg_scaleDown2 += (this_scaleDown - this_val)*(this_scaleDown - this_val);
      bkg_pdfDown2 += (this_pdfDown - this_val)*(this_pdfDown - this_val);
      bkg_topPtDown2 += (this_topPtDown - this_val)*(this_topPtDown - this_val);
      bkg_JECdown2 += (this_JECdown - this_val)*(this_JECdown - this_val);
      bkg_leptonSFdown2 += (this_leptonSFdown - this_val)*(this_leptonSFdown - this_val);
      bkg_photonSFdown2 += (this_photonSFdown - this_val)*(this_photonSFdown - this_val);

      this_syserr2_up = (this_btagUp - this_val)*(this_btagUp - this_val) +
	(this_scaleUp - this_val)*(this_scaleUp - this_val) +
	(this_pdfUp - this_val)*(this_pdfUp - this_val) +
	(this_topPtUp - this_val)*(this_topPtUp - this_val) +
	(this_JECup - this_val)*(this_JECup - this_val) +
	(this_leptonSFup - this_val)*(this_leptonSFup - this_val) +
	(this_photonSFup - this_val)*(this_photonSFup - this_val);

      this_syserr2_down = (this_btagDown - this_val)*(this_btagDown - this_val) +
	(this_scaleDown - this_val)*(this_scaleDown - this_val) +
	(this_pdfDown - this_val)*(this_pdfDown - this_val) +
	(this_topPtDown - this_val)*(this_topPtDown - this_val) +
	(this_JECdown - this_val)*(this_JECdown - this_val) +
	(this_leptonSFdown - this_val)*(this_leptonSFdown - this_val) +
	(this_photonSFdown - this_val)*(this_photonSFdown - this_val);

      bkgsys2_up += this_syserr2_up;
      bkgsys2_down += this_syserr2_down;

      TString fullTableLine = tableNames[j] + "val%%ux:%%.1f\n" + 
	tableNames[j] + "errorup%%ux:%%.2f\n" + 
	tableNames[j] + "errordown%%ux:%%.2f\n";
      char buffer[200];
      sprintf(buffer, fullTableLine.Data());

      fprintf(tableFile, buffer,
	      i+1, this_val, 
	      i+1, sqrt(this_syserr2_up + this_staterr2),
	      i+1, sqrt(this_syserr2_down + this_staterr2));

      if(rangeLow[i] == 100 && rangeHigh[i] == -1) {
	Float_t avg_error_stat = 1. + sqrt(this_staterr2) / this_val;
	Float_t avg_error_btag = 1. + (this_btagUp - this_btagDown) / 2. / this_val;
	Float_t avg_error_scale = 1. + (this_scaleUp - this_scaleDown) / 2. / this_val;
	Float_t avg_error_pdf = 1. + (this_pdfUp - this_pdfDown) / 2. / this_val;
	Float_t avg_error_topPt = 1. + (this_topPtUp - this_topPtDown) / 2. / this_val;
	Float_t avg_error_JEC = 1. + (this_JECup - this_JECdown) / 2. / this_val;
	Float_t avg_error_leptonSF = 1. + (this_leptonSFup - this_leptonSFdown) / 2. / this_val;
	Float_t avg_error_photonSF = 1. + (this_photonSFup - this_leptonSFdown) / 2. / this_val;

	TString fullDatacardLine = tableNames[j] + "val:%%.2f\n" + 
	  tableNames[j] + "stat:%%.2f\n" + 
	  tableNames[j] + "btag:%%.3f\n" +
	  tableNames[j] + "scale:%%.3f\n" +
	  tableNames[j] + "pdf:%%.3f\n" +
	  tableNames[j] + "topPt:%%.3f\n" +
	  tableNames[j] + "JEC:%%.3f\n" +
	  tableNames[j] + "leptonSF:%%.3f\n" +
	  tableNames[j] + "photonSF:%%.3f\n\n";
	sprintf(buffer, fullDatacardLine.Data());

	fprintf(datacardFile, buffer, this_val, 
		avg_error_stat,
		avg_error_btag,
		avg_error_scale,
		avg_error_pdf,
		avg_error_topPt,
		avg_error_JEC,
		avg_error_leptonSF,
		avg_error_photonSF);

      }
	

    }

    // total background
    fprintf(tableFile, "bkgval%dx:%.1f\nbkgerrorup%dx:%.2f\nbkgerrordown%dx:%.2f\n", i+1, bkgval, i+1, sqrt(bkgsys2_up + bkgstat2), i+1, sqrt(bkgsys2_down + bkgstat2));

    // Signal yields
    this_val = h_siga[variableNumber]->IntegralAndError(binLow[i], binHigh[i], this_err);
    this_staterr2 = this_err*this_err;
    
    Double_t temperr;
    this_btagUp = h_siga_btagWeightUp[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_topPtUp = h_siga_topPtUp[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_JECup = h_siga_JECup[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_leptonSFup = h_siga_leptonSFup[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_photonSFup = h_siga_photonSFup[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    
    this_btagDown = h_siga_btagWeightDown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_topPtDown = h_siga_topPtDown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_JECdown = h_siga_JECdown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_leptonSFdown = h_siga_leptonSFdown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_photonSFdown = h_siga_photonSFdown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);

    this_syserr2_up = (this_btagUp - this_val)*(this_btagUp - this_val) +
      .14499 * .14499 * this_val * this_val +
      (this_topPtUp - this_val)*(this_topPtUp - this_val) +
      (this_JECup - this_val)*(this_JECup - this_val) +
      (this_leptonSFup - this_val)*(this_leptonSFup - this_val) +
      (this_photonSFup - this_val)*(this_photonSFup - this_val);
    
    this_syserr2_down = (this_btagDown - this_val)*(this_btagDown - this_val) +
      .14499 * .14499 * this_val * this_val +
      (this_topPtDown - this_val)*(this_topPtDown - this_val) +
      (this_JECdown - this_val)*(this_JECdown - this_val) +
      (this_leptonSFdown - this_val)*(this_leptonSFdown - this_val) +
      (this_photonSFdown - this_val)*(this_photonSFdown - this_val);
    
    fprintf(tableFile, "sigaval%dx:%.1f\nsigaerrorup%dx:%.2f\nsigaerrordown%dx:%.2f\n", i+1, this_val, 
	    i+1, sqrt(this_staterr2 + this_syserr2_up), 
	    i+1, sqrt(this_staterr2 + this_syserr2_down));

    if(rangeLow[i] == 100 && rangeHigh[i] == -1) {
      Float_t avg_error_stat = 1. + sqrt(this_staterr2) / this_val;
      Float_t avg_error_btag = 1. + (this_btagUp - this_btagDown) / 2. / this_val;
      Float_t avg_error_scale = 1. + (this_scaleUp - this_scaleDown) / 2. / this_val;
      Float_t avg_error_pdf = 1. + (this_pdfUp - this_pdfDown) / 2. / this_val;
      Float_t avg_error_topPt = 1. + (this_topPtUp - this_topPtDown) / 2. / this_val;
      Float_t avg_error_JEC = 1. + (this_JECup - this_JECdown) / 2. / this_val;
      Float_t avg_error_leptonSF = 1. + (this_leptonSFup - this_leptonSFdown) / 2. / this_val;
      Float_t avg_error_photonSF = 1. + (this_photonSFup - this_leptonSFdown) / 2. / this_val;
      
      fprintf(datacardFile, "sigaval:%%.2f\nsigastat:%%.2f\nsigabtag:%%.3f\nsigascale:%%.3f\nsigapdf:%%.3f\nsigatopPt:%%.3f\nsigaJEC:%%.3f\nsigaleptonSF:%%.3f\nsigaphotonSF:%%.3f\n\n", 
	      this_val, 
	      avg_error_stat,
	      avg_error_btag,
	      avg_error_scale,
	      avg_error_pdf,
	      avg_error_topPt,
	      avg_error_JEC,
	      avg_error_leptonSF,
	      avg_error_photonSF);
      
    }

    this_val = h_sigb[variableNumber]->IntegralAndError(binLow[i], binHigh[i], this_err);
    this_staterr2 = this_err*this_err;
    
    this_btagUp = h_sigb_btagWeightUp[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_topPtUp = h_sigb_topPtUp[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_JECup = h_sigb_JECup[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_leptonSFup = h_sigb_leptonSFup[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_photonSFup = h_sigb_photonSFup[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    
    this_btagDown = h_sigb_btagWeightDown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_topPtDown = h_sigb_topPtDown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_JECdown = h_sigb_JECdown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_leptonSFdown = h_sigb_leptonSFdown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);
    this_photonSFdown = h_sigb_photonSFdown[variableNumber]->IntegralAndError(binLow[i], binHigh[i], temperr);

    this_syserr2_up = (this_btagUp - this_val)*(this_btagUp - this_val) +
      .16065 * .16065 * this_val * this_val +
      (this_topPtUp - this_val)*(this_topPtUp - this_val) +
      (this_JECup - this_val)*(this_JECup - this_val) +
      (this_leptonSFup - this_val)*(this_leptonSFup - this_val) +
      (this_photonSFup - this_val)*(this_photonSFup - this_val);
    
    this_syserr2_down = (this_btagDown - this_val)*(this_btagDown - this_val) +
      .16065 * .16065 * this_val * this_val +
      (this_topPtDown - this_val)*(this_topPtDown - this_val) +
      (this_JECdown - this_val)*(this_JECdown - this_val) +
      (this_leptonSFdown - this_val)*(this_leptonSFdown - this_val) +
      (this_photonSFdown - this_val)*(this_photonSFdown - this_val);
    
    fprintf(tableFile, "sigbval%dx:%.1f\nsigberrorup%dx:%.2f\nsigberrordown%dx:%.2f\n", i+1, this_val, 
	    i+1, sqrt(this_staterr2 + this_syserr2_up), 
	    i+1, sqrt(this_staterr2 + this_syserr2_down));

    if(rangeLow[i] == 100 && rangeHigh[i] == -1) {
      Float_t avg_error_stat = 1. + sqrt(this_staterr2) / this_val;
      Float_t avg_error_btag = 1. + (this_btagUp - this_btagDown) / 2. / this_val;
      Float_t avg_error_scale = 1. + (this_scaleUp - this_scaleDown) / 2. / this_val;
      Float_t avg_error_pdf = 1. + (this_pdfUp - this_pdfDown) / 2. / this_val;
      Float_t avg_error_topPt = 1. + (this_topPtUp - this_topPtDown) / 2. / this_val;
      Float_t avg_error_JEC = 1. + (this_JECup - this_JECdown) / 2. / this_val;
      Float_t avg_error_leptonSF = 1. + (this_leptonSFup - this_leptonSFdown) / 2. / this_val;
      Float_t avg_error_photonSF = 1. + (this_photonSFup - this_leptonSFdown) / 2. / this_val;
      
      fprintf(datacardFile, "sigbval:%%.2f\nsigbstat:%%.2f\nsigbbtag:%%.3f\nsigbscale:%%.3f\nsigbpdf:%%.3f\nsigbtopPt:%%.3f\nsigbJEC:%%.3f\nsigbleptonSF:%%.3f\nsigbphotonSF:%%.3f\n\n", 
	      this_val, 
	      avg_error_stat,
	      avg_error_btag,
	      avg_error_scale,
	      avg_error_pdf,
	      avg_error_topPt,
	      avg_error_JEC,
	      avg_error_leptonSF,
	      avg_error_photonSF);
      
    }

    // Data
    this_val = h_gg[variableNumber]->IntegralAndError(binLow[i], binHigh[i], this_err);
    fprintf(tableFile, "dataval%dx:%.0f\n", i+1, this_val);

    if(rangeLow[i] == 100 && rangeHigh[i] == -1) {
      fprintf(tableFile, "bkgstat5y:%.1f\n", 100. * sqrt(bkgstat2) / bkgval);
      fprintf(tableFile, "bkgsysup5y:%.1f\nbkgsysdown5y:%.1f\n", 100. * sqrt(bkgsys2_up) / bkgval, 100. * sqrt(bkgsys2_down) / bkgval);
      fprintf(tableFile, "bkgscaleup5y:%.1f\nbkgscaledown5y:%.1f\n", 100. * sqrt(bkg_scaleUp2) / bkgval, 100. * sqrt(bkg_scaleDown2) / bkgval);
      fprintf(tableFile, "bkgpdfup5y:%.1f\nbkgpdfdown5y:%.1f\n", 100. * sqrt(bkg_pdfUp2) / bkgval, 100. * sqrt(bkg_pdfDown2) / bkgval);
      fprintf(tableFile, "bkgjecup5y:%.1f\nbkgjecdown5y:%.1f\n", 100. * sqrt(bkg_JECup2) / bkgval, 100. * sqrt(bkg_JECdown2) / bkgval);
      fprintf(tableFile, "bkgtopptup5y:%.1f\nbkgtopptdown5y:%.1f\n", 100. * sqrt(bkg_topPtUp2) / bkgval, 100. * sqrt(bkg_topPtDown2) / bkgval);
      fprintf(tableFile, "bkgbtagup5y:%.1f\nbkgbtagdown5y:%.1f\n", 100. * sqrt(bkg_btagUp2) / bkgval, 100. * sqrt(bkg_btagDown2) / bkgval);
      fprintf(tableFile, "bkgleptonSFup5y:%.1f\nbkgleptonSFdown5y:%.1f\n", 100. * sqrt(bkg_leptonSFup2) / bkgval, 100. * sqrt(bkg_leptonSFdown2) / bkgval);
      fprintf(tableFile, "bkgphotonSFup5y:%.2f\nbkgphotonSFdown5y:%.2f\n", 100. * sqrt(bkg_photonSFup2) / bkgval, 100. * sqrt(bkg_photonSFdown2) / bkgval);
    }

  }

  fclose(tableFile);
  fclose(datacardFile);

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

void PlotMaker::GetLeptonSF(vector<Float_t> vars, int chan, Float_t& central, Float_t& up, Float_t& down) {

  Float_t pt, eta, error;

  if(chan < 2) {
    pt = min(vars[15], (float)199.);
    pt = max(pt, (float)15.);
    eta = min(fabs(vars[16]), (double)2.39);

    Float_t id_val = sf_electron->GetBinContent(sf_electron->FindBin(eta, pt));
    Float_t id_error = sf_electron->GetBinError(sf_electron->FindBin(eta, pt));

    Float_t trigger_val = sf_SingleElectronTrigger->GetBinContent(sf_SingleElectronTrigger->FindBin(eta, pt));
    Float_t trigger_error = sf_SingleElectronTrigger->GetBinError(sf_SingleElectronTrigger->FindBin(eta, pt));

    central = id_val * trigger_val;
    error = central * sqrt(id_error*id_error/(id_val*id_val) + trigger_error*trigger_error/(trigger_val*trigger_val));

    up = central + error;
    down = central - error;
  }

  else {
    pt = min(vars[17], (float)499.);
    pt = max(pt, (float)10.);
    eta = min(fabs(vars[18]), (double)2.09);

    central = sf_muon->GetBinContent(sf_muon->FindBin(pt, eta));
    error = sf_muon->GetBinError(sf_muon->FindBin(pt, eta));

    up = error;
    down = 2. * central - error;
  }

  return;  
}

void PlotMaker::GetPhotonSF(vector<Float_t> vars, Float_t& central, Float_t& up, Float_t& down) {

  if(vars[0] == 0) {
    central = 1.;
    up = 1.;
    down = 1.;
    return;
  }

  Float_t et, eta, error;

  if(vars[0] == 1 && vars.size() >= 21) {
    et = min(vars[19], (float)999.);
    et = min(et, (float)15.);
    eta = min(fabs(vars[20]), (double)2.49);

    Float_t id_val = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    central = id_val * veto_val;
    error = central * sqrt(id_error*id_error/(id_val*id_val) + veto_error*veto_error/(veto_val*veto_val));
  }

  else if(vars[0] >= 2 && vars.size() >= 27) {
    // lead photon
    et = min(vars[19], (float)999.);
    et = min(et, (float)15.);
    eta = min(fabs(vars[20]), (double)2.49);

    Float_t id_val_lead = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error_lead = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val_lead = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error_lead = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    // trail photon
    et = min(vars[24], (float)999.);
    et = min(et, (float)15.);
    eta = min(fabs(vars[26]), (double)2.49);

    Float_t id_val_trail = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error_trail = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val_trail = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error_trail = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    central = id_val_lead * veto_val_lead * id_val_trail * veto_val_trail;
    error = central * sqrt(id_error_lead*id_error_lead/(id_val_lead*id_val_lead) +
			   veto_error_lead*veto_error_lead/(veto_val_lead*veto_val_lead) +
			   id_error_trail*id_error_trail/(id_val_trail*id_val_trail) +
			   veto_error_trail*veto_error_trail/(veto_val_trail*veto_val_trail));
  }

  up = central + error;
  down = central - error;

  return;
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
    if(n > 0) h_acc->Fill(index1, index2, acceptance / (0.438/3.) / n);

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
