#include "TH1.h"
#include "TString.h"

#include <vector>
#include <stdarg.h>

using namespace std;

//const TString gifOrPdf = ".pdf";

bool LargerHistogram(const TH1D* h1, const TH1D* h2) {
  return (h1->GetMaximum() > h2->GetMaximum());
}

void OverlayHistograms(int n, TString title, bool logy, bool scale, TString fileName, TH1D * hist, ...) {

  if(n < 1) return;

  va_list histos;
  va_start(histos, hist);

  vector<TH1D*> h;
  h.push_back( (TH1D*)hist->Clone() );
  
  for(int i = 1; i < n; i++) {
    h.push_back( (TH1D*)va_arg(histos, TH1D*)->Clone() );
  }
  va_end(histos);
 
  sort(h.begin(), h.end(), LargerHistogram);
  
  TLegend * leg = new TLegend(0.7, 0.7, 0.87, 0.87);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);

  TCanvas * canv = new TCanvas("canv", "Plot", 10, 10, 2000, 2000);

  Int_t furthestBin = -1;

  for(unsigned int i = 0; i < h.size(); i++) {

    for(int j = h[i]->GetNbinsX(); j > 0; j--) {
      if(h[i]->GetBinContent(j) > 0) {
        if(j > furthestBin) furthestBin = j;
	break;
      }
    }

    if(scale) h[i]->Scale(1./(h[i]->Integral()));

    if(i == 0) {
      //h[i]->SetTitle(title);
      h[i]->Draw();
    }
    else {
      h[i]->SetLineColor(i+1);
      h[i]->Draw("same");
    }

    leg->AddEntry(h[i], "", "PL");
  }

  leg->Draw("same");
  
  h[0]->SetTitle(title);
  h[0]->GetXaxis()->SetRangeUser(0, 1.3 * h[0]->GetXaxis()->GetBinCenter(furthestBin));

  canv->SetLogy(logy);  
  canv->SaveAs(fileName);

  delete canv;
  delete leg;
}

TH1D * DivideByBinWidth(TH1D * h) {

  for(Int_t i = 0; i < h->GetNbinsX(); i++) {
    Double_t val = h->GetBinContent(i+1);
    Double_t err = h->GetBinError(i+1);
    Double_t width = h->GetBinWidth(i+1);
    
    h->SetBinContent(i+1, val / width);
    h->SetBinError(i+1, err / width);
  }

  return h;
}

