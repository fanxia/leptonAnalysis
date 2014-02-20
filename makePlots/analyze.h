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

TH1D * HistoFromTree(TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double metCut = -1.) {
  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();
  FillHistoFromTree(h, variable, metCut);
  return h;
}
  
TH1D * HistoFromTree(TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t* customBins, double metCut = -1.) {
  TH1D * h = new TH1D(name, title, nBins, customBins);
  h->Sumw2();
  FillHistoFromTree(h, tree, variable, metCut);
  return h;
}

TH1D * SignalHistoFromTree(Float_t scale, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double metCut = -1.) {
  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();
  FillSignalHistoFromTree(h, tree, variable, metCut, scale);
  return h;
}

TH1D * SignalHistoFromTree(Float_t scale, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t* customBins, double metCut = -1.) {
  TH1D * h = new TH1D(name, title, nBins, customBins);
  h->Sumw2();
  FillSignalHistoFromTree(h, tree, variable, metCut, scale);
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
  virtual ~PlotMaker() { 

    KSscores.clear();

    delete ggTree;
    delete ttHadronicTree;
    delete ttSemiLepTree;
    delete ttFullLepTree;
    delete wjetsTree;
    delete dyjetsTree;
    delete ttgjetsTree;
    delete ttggTree;
    delete sigaTree;
    delete sigbTree;
    
  }

  void SetTrees(TTree * gg,
		TTree * ttHadronic, TTree * ttSemiLep, TTree * ttFullLepTree,
		TTree * wjets, TTree * dy,
		TTree * ttgjets, TTree * ttgg,
		TTree * ttbar,
		TTree * sig_a, TTree * sig_b);

  void SetDisplayKStest(bool v) { displayKStest = v; }

  void CreatePlot(TString variable,
		  Int_t nBinsX, Float_t bin_lo, Float_t bin_hi,
		  TString xaxisTitle, TString yaxisTitle,
		  Float_t xmin, Float_t xmax,
		  Float_t ymin, Float_t ymax,
		  Float_t ratiomin, Float_t ratiomax,
		  bool drawSignal, bool drawLegend, bool drawPrelim,
		  TFile*& out, double metCut);

  void CreatePlot(TString variable,
		  Int_t nBinsX, Double_t* customBins,
		  TString xaxisTitle,
		  Float_t xmin, Float_t xmax,
		  Float_t ymin, Float_t ymax,
		  Float_t ratiomin, Float_t ratiomax,
		  bool drawSignal, bool drawLegend, bool drawPrelim,
		  TFile*& out, double metCut);

  void DrawPlot(TH1D * gg,
		TH1D * ttHadronic, TH1D * ttSemiLep, TH1D * ttFullLep,
		TH1D * wjets, TH1D * dyjets,
		TH1D * ttg, TH1D * ttgg,
		TH1D * sig_a, TH1D * sig_b,
		TString variable,
		TString xaxisTitle, TString yaxisTitle,
		Float_t xmin, Float_t xmax,
		Float_t ymin, Float_t ymax,
		Float_t ratiomin, Float_t ratiomax,
		bool drawSignal, bool drawLegend, bool drawPrelim,
		TFile*& out) {

  void CreateTable();

  void PlotKolmogorovValues();

 private:
  TTree * ggTree;
  
  TTree * ttHadronicTree;
  TTree * ttSemiLepTree;
  TTree * ttFullLepTree;

  TTree * wjetsTree;
  TTree * dyjetsTree;

  TTree * ttgjetsTree;
  TTree * ttggTree;

  TTree * sigaTree;
  TTree * sigbTree;

  Int_t intLumi_int;
  TString intLumi;
  TString req;

  bool displayKStest;
  bool blinded;

  vector<pair<TString, double> > KSscores;
  
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
}

void PlotMaker::SetTrees(TTree * gg,
			 TTree * ttHadronic, TTree * ttSemiLep, TTree * ttFullLepTree,
			 TTree * wjets, TTree * dy,
			 TTree * ttgjets, TTree * ttgg,
			 TTree * ttbar,
			 TTree * sig_a, TTree * sig_b) {

  ggTree = gg;
  
  ttHadronicTree = ttHadronic;
  ttSemiLepTree = ttSemiLep;
  ttFullLepTree = ttFullLep;

  wjetsTree = wjets;
  dyjetsTree = dy;

  ttgjetsTree = ttgjets;
  ttggTree = ttgg;

  sigaTree = sig_a;
  sigbTree = sig_b;

}

void PlotMaker::CreatePlot(TString variable,
			   Int_t nBinsX, Float_t bin_lo, Float_t bin_hi,
			   TString xaxisTitle, TString yaxisTitle,
			   Float_t xmin, Float_t xmax,
			   Float_t ymin, Float_t ymax,
			   Float_t ratiomin, Float_t ratiomax,
			   bool drawSignal, bool drawLegend, bool drawPrelim,
			   TFile*& out, double metCut) {

  TH1D * gg = HistoFromTree(variable, ggTree, variable+"_gg_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);

  TH1D * ttHadronic = SignalHistoFromTree(intLumi_int * 53.4 * 1.019 * 1.019 / 10537444., variable, ttHadronicTree, variable+"_ttHadronic_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * ttSemiLep = SignalHistoFromTree(intLumi_int * 53.2 * 1.019 * 1.019 / 25424818., variable, ttSemiLepTree, variable+"_ttSemiLep_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * ttFullLep = SignalHistoFromTree(intLumi_int * 13.43 * 1.019 * 1.019 / 12119013., variable, ttFullLepTree, variable+"_ttFullLep_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);

  TH1D * wjets = SignalHistoFromTree(intLumi_int *  1.019 * 1.019 * 37509. / 18393090., variable, wjetsTree, variable+"_wjets_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * dyjets = SignalHistoFromTree(intLumi_int *  1.019 * 1.019 * 3504. / 30459503., variable, dyjetsTree, variable+"_dyjets_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);

  TH1D * ttg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 14.0 / 1719954., variable, ttgjetsTree, variable+"_ttgjets_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * ttgg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 0.146 / 1719954., variable, ttggTree, variable+"_ttgg_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);

  TH1D * sig_a;
  TH1D * sig_b;
  if(drawSignal) {
    sig_a = SignalHistoFromTree(0.147492 * intLumi_int * 1.019 * 1.019 / 15000., variable, sigaTree, variable+"_a_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
    sig_b = SignalHistoFromTree(0.0399591 * intLumi_int * 1.019 * 1.019 / 15000., variable, sigbTree, variable+"_b_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  }

  DrawPlot(gg,
	   ttHadronic, ttSemiLep, ttFullLep,
	   wjets, dyjets,
	   ttg, ttgg,
	   sig_a, sig_b,
	   variable,
	   xaxisTitle, yaxisTitle,
	   xmin, xmax,
	   ymin, ymax,
	   ratiomin, ratiomax,
	   drawSignal, drawLegend, drawPrelim,
	   out);

}

void PlotMaker::CreatePlot(TString variable,
			   Int_t nBinsX, Double_t* customBins,
			   TString xaxisTitle,
			   Float_t xmin, Float_t xmax,
			   Float_t ymin, Float_t ymax,
			   Float_t ratiomin, Float_t ratiomax,
			   bool drawSignal, bool drawLegend, bool drawPrelim,
			   TFile*& out, double metCut) {

  TString yaxisTitle = "Number of Events / GeV";

  TH1D * gg = HistoFromTree(variable, ggTree, variable+"_gg_"+req, variable, nBinsX, customBins, metCut); 
  gg = (TH1D*)DivideByBinWidth(gg);

  TH1D * ttHadronic = SignalHistoFromTree(intLumi_int * 53.4 * 1.019 * 1.019 / 10537444., variable, ttHadronicTree, variable+"_ttHadronic_"+req, variable, nBinsX, customBins, metCut);
  ttHadronic = (TH1D*)DivideByBinWidth(ttHadronic);
  TH1D * ttSemiLep = SignalHistoFromTree(intLumi_int * 53.2 * 1.019 * 1.019 / 25424818., variable, ttSemiLepTree, variable+"_ttSemiLep_"+req, variable, nBinsX, customBins, metCut);
  ttSemiLep = (TH1D*)DivideByBinWidth(ttSemiLep);
  TH1D * ttFullLep = SignalHistoFromTree(intLumi_int * 13.43 * 1.019 * 1.019 / 12119013., variable, ttFullLepTree, variable+"_ttFullLep_"+req, variable, nBinsX, customBins, metCut);
  ttFullLep = (TH1D*)DivideByBinWidth(ttFullLep);

  TH1D * wjets = SignalHistoFromTree(intLumi_int *  1.019 * 1.019 * 37509. / 18393090., variable, wjetsTree, variable+"_wjets_"+req, variable, nBinsX, customBins, metCut);
  wjets = (TH1D*)DivideByBinWidth(wjets);
  TH1D * dyjets = SignalHistoFromTree(intLumi_int *  1.019 * 1.019 * 3504. / 30459503., variable, dyjetsTree, variable+"_dyjets_"+req, variable, nBinsX, customBins, metCut);
  dyjets = (TH1D*)DivideByBinWidth(dyjets);

  TH1D * ttg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 14.0 / 1719954., variable, ttgjetsTree, variable+"_ttgjets_"+req, variable, nBinsX, customBins, metCut);
  ttg = (TH1D*)DivideByBinWidth(ttg);
  TH1D * ttgg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 0.146 / 1719954., variable, ttggTree, variable+"_ttgg_"+req, variable, nBinsX, customBins, metCut);
  ttgg = (TH1D*)DivideByBinWidth(ttgg);
  
  TH1D * sig_a;
  TH1D * sig_b;
  if(drawSignal) {
    sig_a = SignalHistoFromTree(0.147492 * intLumi_int * 1.019 * 1.019 / 15000., true, variable, sigaTree, variable+"_a_"+req, variable, nBinsX, customBins, metCut);
    sig_a = (TH1D*)DivideByBinWidth(sig_a);
    sig_b = SignalHistoFromTree(0.0399591 * intLumi_int * 1.019 * 1.019 / 15000., true, variable, sigbTree, variable+"_b_"+req, variable, nBinsX, customBins, metCut);
    sig_b = (TH1D*)DivideByBinWidth(sig_b);
  }

  DrawPlot(gg,
	   ttHadronic, ttSemiLep, ttFullLep,
	   wjets, dyjets,
	   ttg, ttgg,
	   sig_a, sig_b,
	   variable,
	   xaxisTitle, yaxisTitle,
	   xmin, xmax,
	   ymin, ymax,
	   ratiomin, ratiomax,
	   drawSignal, drawLegend, drawPrelim,
	   out);

}

void PlotMaker::DrawPlot(TH1D * gg,
			 TH1D * ttHadronic, TH1D * ttSemiLep, TH1D * ttFullLep,
			 TH1D * wjets, TH1D * dyjets,
			 TH1D * ttg, TH1D * ttgg,
			 TH1D * sig_a, TH1D * sig_b,
			 TString variable,
			 TString xaxisTitle, TString yaxisTitle,
			 Float_t xmin, Float_t xmax,
			 Float_t ymin, Float_t ymax,
			 Float_t ratiomin, Float_t ratiomax,
			 bool drawSignal, bool drawLegend, bool drawPrelim,
			 TFile*& out) {

  if(blinded) for(int i = 0; i < gg->GetNbinsX(); i++) {
      gg->SetBinContent(i+1, 1.e-14);
      gg->SetBinError(i+1, 1.e-14);
    }

  TH1D * ttbar = (TH1D*)ttHadronic->Clone(variable+"_ttbar_"+req);
  ttbar->Add(ttSemiLep);
  ttbar->Add(ttFullLep);

  TH1D * bkg = (TH1D*)ttbar->Clone(variable+"_bkg_"+req);

  bkg->Add(wjets);
  bkg->Add(dyjets);
  bkg->Add(ttg);
  bkg->Add(ttgg);
  
  Double_t kolm = gg->KolmogorovTest(bkg);
  TString kolmText = Form("KS test probability = %5.3g", kolm);
  TText * tt = new TText(0.92, 0.5, kolmText);
  tt->SetTextAngle(90.);
  tt->SetNDC(); tt->SetTextSize( 0.032 );

  out->cd();
  gg->Write();
  ttbar->Write();
  wjets->Write();
  dyjets->Write();
  ttg->Write();
  ttgg->Write();
  bkg->Write();

  wjets->Add(dyjets);
  wjets->Add(ttg);
  wjets->Add(ttgg);
  
  dyjets->Add(ttg);
  dyjets->Add(ttgg);

  ttg->Add(ttgg);

  TH1D * errors = (TH1D*)bkg->Clone("errors");

  if(drawSignal) calculateROC(sig_a, sig_b, bkg, req, variable);

  TLegend * leg = new TLegend(0.50, 0.65, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(gg, "#gamma#gamma Candidate Sample", "LP");
  leg->AddEntry(errors, "Total Background Uncertainty", "F");
  leg->AddEntry(bkg, "t#bar{t} inclusive", "F");
  leg->AddEntry(wjets, "W + Jets", "F");
  leg->AddEntry(dyjets, "DY + Jets", "F");
  leg->AddEntry(ttg, "t#bar{t}+#gamma", "F");
  leg->AddEntry(ttgg, "t#bar{t}+#gamma#gamma", "F");
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

  gg->SetMarkerStyle(20); 
  gg->SetMarkerSize(1.5);

  errors->SetFillColor(kOrange+10);
  errors->SetFillStyle(3154);
  errors->SetMarkerSize(0);

  // new stack: ttbar, wjets, dyjets, ttg, ttgg
  bkg->SetFillColor(kGray);
  bkg->SetMarkerSize(0);
  bkg->SetLineColor(1);

  wjets->SetFillColor(kOrange-3);
  wjets->SetMarkerSize(0);
  wjets->SetLineColor(1);

  dyjets->SetFillColor(kYellow);
  dyjets->SetMarkerSize(0);
  dyjets->SetLineColor(1);

  ttg->SetFillColor(kViolet);
  ttg->SetMarkerSize(0);
  ttg->SetLineColor(1);

  ttgg->SetFillColor(8);
  ttgg->SetMarkerSize(0);
  ttgg->SetLineColor(1);

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

  // new stack: ttbar, wjets, dyjets, ttg, ttgg
  bkg->Draw("hist");
  wjets->Draw("same hist");
  dyjets->Draw("same hist");
  ttg->Draw("same hist");
  ttgg->Draw("same hist");
  errors->Draw("same e2");
  gg->Draw("same e1");
  bkg->Draw("same axis");

  if(drawSignal) {
    sig_a->SetLineColor(kMagenta);
    sig_a->SetLineWidth(3);
    leg->AddEntry(sig_a, "GGM #gamma#gamma (460_175)", "L");
    sig_a->Draw("same hist");
    
    sig_b->SetLineColor(kBlue);
    sig_b->SetLineWidth(3);
    leg->AddEntry(sig_b, "GGM #gamma#gamma (560_325)", "L");
    sig_b->Draw("same hist");
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

  TH1D * ratio = (TH1D*)gg->Clone("ratio");
  ratio->Reset();
  ratio->SetTitle("Data / Background");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(bkg->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, gg->GetBinContent(i+1) / bkg->GetBinContent(i+1));
    ratio->SetBinError(i+1, gg->GetBinError(i+1) / bkg->GetBinContent(i+1));
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

  const int nBins = 5;
  Double_t xbins[nBins+1] = {0, 20, 50, 80, 100, 1000};

  Double_t rangeLow[nBins] = {0, 0, 50, 80, 100};
  Double_t rangeHigh[nBins] = {20, 50, -1, -1, -1};

  TH1D * h_gg = HistoFromTree(true, "pfMET", ggTree, "pfMet2_gg_"+req, "pfMet2", nBins, xbins, -1.);
  TH1D * h_ewk = HistoFromTree(true, "pfMET", egTree, "pfMet2_eg_"+req, "pfMet2", nBins, xbins, -1.);
  
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

void prep_signal(TString req) {

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
  TH2D * h_contam = new TH2D("contam_"+req, "contam_"+req, 30, xbins, 32, ybins);

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
    
    TTree * ggTree = (TTree*)f->Get("gg_"+req+"_EvtTree"+code_t);
    TTree * gfTree = (TTree*)f->Get("gf_"+req+"_EvtTree"+code_t);
    TTree * ffTree = (TTree*)f->Get("ff_"+req+"_EvtTree"+code_t);

    TH1D * gg;
    TH1D * gf;
    TH1D * ff;

    if(ggTree->GetEntries() > 0) {
      gg = (TH1D*)SignalHistoFromTree(1.0, true, "pfMET", ggTree, "met_gg_"+req+code_t, "met_gg_"+req+code_t, 400, 0., 2000.);
      gf = (TH1D*)SignalHistoFromTree(1.0, true, "pfMET", gfTree, "met_gf_"+req+code_t, "met_gf_"+req+code_t, 400, 0., 2000.);
      ff = (TH1D*)SignalHistoFromTree(1.0, true, "pfMET", ffTree, "met_ff_"+req+code_t, "met_ff_"+req+code_t, 400, 0., 2000.);

      out->cd();
      gg->Write();
      gf->Write();
      ff->Write();
    }
    else {
      f->Close();
      continue;
    }

    //TH1D * h = (TH1D*)f->Get("ngen"+code_t);
    //double n = h->GetBinContent(1);
    double n = 15000.;

    double acceptance = gg->Integral();
    if(n > 0) h_acc->Fill(index1, index2, acceptance / n);

    double contamination = ff->Integral();
    if(n > 0) h_contam->Fill(index1, index2, contamination / acceptance);

    f->Close();
  }

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  h_acc->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_acc->GetXaxis()->SetRangeUser(0, 1600);
  h_acc->GetXaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetTitle("Bino} mass (GeV/c^{2})");
  h_acc->GetYaxis()->SetTitleOffset(1.3);
  h_acc->GetYaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetRangeUser(0, 1600);
  h_acc->GetZaxis()->SetLabelSize(0.02);
  h_acc->Draw("colz");
  can->SaveAs("acceptance_"+req+".pdf");
  
  h_contam->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_contam->GetXaxis()->SetRangeUser(0, 1600);
  h_contam->GetXaxis()->SetLabelSize(0.03);
  h_contam->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_contam->GetYaxis()->SetTitleOffset(1.3);
  h_contam->GetYaxis()->SetLabelSize(0.03);
  h_contam->GetYaxis()->SetRangeUser(0, 1600);
  h_contam->GetZaxis()->SetLabelSize(0.02);
  h_contam->Draw("colz");
  can->SaveAs("contamination_"+req+".pdf");

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
