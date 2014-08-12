using namespace std;
using namespace RooFit;

void makeFit(TString varname, double varmin, double varmax, TH1D * signalHist, TH1D * backgroundHist, TH1D * dataHist, double& value, double& error) {

  //RooFit variables
  RooRealVar var(varname, varname, varmin, varmax);

  // create PDFs
  RooDataHist signalDataHist("signalDataHist", "signal RooDataHist", RooArgList(var), signalHist);
  RooHistPdf signalPdf("signalPdf", varname+" of signal", RooArgSet(var), signalDataHist);

  RooDataHist backgroundDataHist("backgroundDataHist", "background RooDataHist", RooArgList(var), backgroundHist);
  RooHistPdf backgroundPdf("backgroundPdf", varname+" of background", RooArgSet(var), backgroundDataHist);

  // data
  RooDataHist dataDataHist("data "+varname, varname+"in Data", RooArgList(var), dataHist);

  // signal fraction parameter
  RooRealVar signalFractionVar("signal fraction", "signal fraction", 0.5, 0.0, 1.0);
  RooAddPdf sumPdf("totalPdf", "signal and background", signalPdf, backgroundPdf, signalFractionVar);

  // fit
  sumPdf.fitTo(dataDataHist, RooFit::SumW2Error(kFALSE), RooFit::PrintLevel(-1));
  
  value = signalFractionVar.getVal();
  error = signalFractionVar.getError();
  
  return;
}

TH1D * fillDataHistogram(TString variable, TTree * tree, int nPhotons_req,
			 TString title, Int_t nBins, Double_t xlo, Double_t xhi) {
			
  Float_t Nphotons, z_invmass, var;

  tree->SetBranchAddress("Nphotons", &Nphotons);
  tree->SetBranchAddress("z_invmass", &z_invmass);
  if(variable != "z_invmass") tree->SetBranchAddress(variable, &var);

  TH1D * h = new TH1D(title, title, nBins, xlo, xhi);
  h->Sumw2();

  for(Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if((int)Nphotons != nPhotons_req) continue;
    if(variable != "z_invmass" && (z_invmass <= 76 || z_invmass >= 106)) continue;

    if(variable == "z_invmass") var = z_invmass;

    h->Fill(var);
  }

  return h;
}

TH1D * fillMCHistogram(TString variable, TTree * tree, int nPhotons_req,
		       TString title, Int_t nBins, Double_t xlo, Double_t xhi,
		       double scale) {

  Float_t pileupWeight, pileupWeightErr,
    btagWeight, btagWeightErr,
    Nphotons, z_invmass, var;

  tree->SetBranchAddress("pileupWeight", &pileupWeight);
  tree->SetBranchAddress("pileupWeightErr", &pileupWeightErr);
  tree->SetBranchAddress("btagWeight", &btagWeight);
  tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  tree->SetBranchAddress("Nphotons", &Nphotons);
  tree->SetBranchAddress("z_invmass", &z_invmass);
  if(variable != "z_invmass") tree->SetBranchAddress(variable, &var);

  TH1D * h = new TH1D(title, title, nBins, xlo, xhi);
  h->Sumw2();

  for(Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(btagWeight != btagWeight) continue;

    if((int)Nphotons != nPhotons_req) continue;
    if(variable != "z_invmass" && (z_invmass <= 76 || z_invmass >= 106)) continue;

    if(variable == "z_invmass") var = z_invmass;

    Float_t olderror = h->GetBinError(h->FindBin(var));
    h->Fill(var, pileupWeight * btagWeight);

    // protection from weird 1200 weight errors...
    if(btagWeightErr > 20.) btagWeightErr = btagWeight;

    Float_t addError2 = pileupWeight*pileupWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*pileupWeightErr*pileupWeightErr;
    Float_t newerror = sqrt(olderror*olderror + addError2);

    h->SetBinError(h->FindBin(var), newerror);
  }

  h->Scale(scale);

  tree->ResetBranchAddresses();

  return h;
}

void plotZGG(bool runElectrons = false) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString channel = (runElectrons) ? "ele" : "muon";

  TString type = "noSigmaIetaIeta";

  TFile * fData;
  if(runElectrons) fData = new TFile("condor/SingleElectron_zg.root");
  else fData = new TFile("condor/SingleMuon_zg.root");
  
  TTree * dataTree = (TTree*)fData->Get(channel+"_"+type+"Tree");

  TFile * fDY1 = new TFile("condor/zg_"+channel+"/signal_contamination_dy1JetsToLL.root");
  TFile * fDY2 = new TFile("condor/zg_"+channel+"/signal_contamination_dy2JetsToLL.root");
  TFile * fDY3 = new TFile("condor/zg_"+channel+"/signal_contamination_dy3JetsToLL.root");
  TFile * fDY4 = new TFile("condor/zg_"+channel+"/signal_contamination_dy4JetsToLL.root");
  TFile * fZG = new TFile("condor/zg_"+channel+"/signal_contamination_ZGToLLG.root");

  TFile * fFullLep = new TFile("condor/zg_"+channel+"/signal_contamination_ttJetsFullLep.root");
  TFile * fSemiLep = new TFile("condor/zg_"+channel+"/signal_contamination_ttJetsSemiLep.root");
  TFile * fHadronic = new TFile("condor/zg_"+channel+"/signal_contamination_ttJetsHadronic.root");
  
  TFile * fWW = new TFile("condor/zg_"+channel+"/signal_contamination_WW.root");
  TFile * fWZ = new TFile("condor/zg_"+channel+"/signal_contamination_WZ.root");
  TFile * fZZ = new TFile("condor/zg_"+channel+"/signal_contamination_ZZ.root");

  TTree * tDY1 = (TTree*)fDY1->Get(channel+"_"+type+"Tree");
  TTree * tDY2 = (TTree*)fDY2->Get(channel+"_"+type+"Tree");
  TTree * tDY3 = (TTree*)fDY3->Get(channel+"_"+type+"Tree");
  TTree * tDY4 = (TTree*)fDY4->Get(channel+"_"+type+"Tree");
  TTree * tZG = (TTree*)fZG->Get(channel+"_"+type+"Tree");

  TTree * tFullLep = (TTree*)fFullLep->Get(channel+"_"+type+"Tree");
  TTree * tSemiLep = (TTree*)fSemiLep->Get(channel+"_"+type+"Tree");
  TTree * tHadronic = (TTree*)fHadronic->Get(channel+"_"+type+"Tree");

  TTree * tWW = (TTree*)fWW->Get(channel+"_"+type+"Tree");
  TTree * tWZ = (TTree*)fWZ->Get(channel+"_"+type+"Tree");
  TTree * tZZ = (TTree*)fZZ->Get(channel+"_"+type+"Tree");

  TH1D * hData = fillDataHistogram("leadPhoton_sIetaIeta", dataTree, 1,
				   "data_sIetaIeta", 100, 0., 0.05);

  TH1D * ngen_dy1 = (TH1D*)fDY1->Get("nEvents_dy1JetsToLL");
  double scale_dy1 = 666.7 * 1177.3 * 3 / 3503.71 * 19712. / ngen_dy1->Integral();
  TH1D * hDY1 = fillMCHistogram("leadPhoton_sIetaIeta", tDY1, 1,
				"dy1_sIetaIeta", 100, 0., 0.05,
				scale_dy1);

  TH1D * ngen_dy2 = (TH1D*)fDY2->Get("nEvents_dy2JetsToLL");
  double scale_dy2 = 215.1 * 1177.3 * 3 / 3503.71 * 19712. / ngen_dy2->Integral();
  TH1D * hDY2 = fillMCHistogram("leadPhoton_sIetaIeta", tDY2, 1,
				"dy2_sIetaIeta", 100, 0., 0.05,
				scale_dy2);

  TH1D * ngen_dy3 = (TH1D*)fDY3->Get("nEvents_dy3JetsToLL");
  double scale_dy3 = 66.07 * 1177.3 * 3 / 3503.71 * 19712. / ngen_dy3->Integral();
  TH1D * hDY3 = fillMCHistogram("leadPhoton_sIetaIeta", tDY3, 1,
				"dy3_sIetaIeta", 100, 0., 0.05,
				scale_dy3);

  TH1D * ngen_dy4 = (TH1D*)fDY4->Get("nEvents_dy4JetsToLL");
  double scale_dy4 = 27.38 * 1177.3 * 3 / 3503.71 * 19712. / ngen_dy4->Integral();
  TH1D * hDY4 = fillMCHistogram("leadPhoton_sIetaIeta", tDY4, 1,
				"dy4_sIetaIeta", 100, 0., 0.05,
				scale_dy4);

  TH1D * hDY = (TH1D*)hDY1->Clone("dy");
  hDY->Add(hDY2);
  hDY->Add(hDY3);
  hDY->Add(hDY4);

  TH1D * ngen_zg = (TH1D*)fZG->Get("nEvents_ZGToLLG");
  double scale_zg = 132.6 * 19712. / ngen_zg->Integral();
  TH1D * hZG = fillMCHistogram("leadPhoton_sIetaIeta", tZG, 1,
			       "zg_sIetaIeta", 100, 0., 0.05,
			       scale_zg);

  TH1D * ngen_fullLep = (TH1D*)fFullLep->Get("nEvents_ttJetsFullLep");
  double scale_fullLep = 245.8 * 0.105 * 19712. / ngen_fullLep->Integral();
  TH1D * hFullLep = fillMCHistogram("leadPhoton_sIetaIeta", tFullLep, 1,
				    "fullLep_sIetaIeta", 100, 0., 0.05,
				    scale_fullLep);
  
  TH1D * ngen_semiLep = (TH1D*)fSemiLep->Get("nEvents_ttJetsSemiLep");
  double scale_semiLep = 245.8 * 0.438 * 19712. / ngen_semiLep->Integral();
  TH1D * hSemiLep = fillMCHistogram("leadPhoton_sIetaIeta", tSemiLep, 1,
				    "semiLep_sIetaIeta", 100, 0., 0.05,
				    scale_semiLep);
  
  TH1D * ngen_hadronic = (TH1D*)fHadronic->Get("nEvents_ttJetsHadronic");
  double scale_hadronic = 245.8 * 0.457 * 19712. / ngen_hadronic->Integral();
  TH1D * hHadronic = fillMCHistogram("leadPhoton_sIetaIeta", tHadronic, 1,
				     "hadronic_sIetaIeta", 100, 0., 0.05,
				     scale_hadronic);
  
  TH1D * ngen_WW = (TH1D*)fWW->Get("nEvents_WW");
  double scale_WW = 57.1097 * 19712. / ngen_WW->Integral();
  TH1D * hWW = fillMCHistogram("leadPhoton_sIetaIeta", tWW, 1,
			       "WW_sIetaIeta", 100, 0., 0.05,
			       scale_WW);

  TH1D * ngen_WZ = (TH1D*)fWZ->Get("nEvents_WZ");
  double scale_WZ = 32.3161 * 19712. / ngen_WZ->Integral();
  TH1D * hWZ = fillMCHistogram("leadPhoton_sIetaIeta", tWZ, 1,
			       "WZ_sIetaIeta", 100, 0., 0.05,
			       scale_WZ);
  
  TH1D * ngen_ZZ = (TH1D*)fZZ->Get("nEvents_ZZ");
  double scale_ZZ = 8.25561 * 19712. / ngen_ZZ->Integral();
  TH1D * hZZ = fillMCHistogram("leadPhoton_sIetaIeta", tZZ, 1,
			       "ZZ_sIetaIeta", 100, 0., 0.05,
			       scale_ZZ);

  /*
    hData->Add(hFullLep, -1.);
    hData->Add(hSemiLep, -1.);
    hData->Add(hHadronic, -1.);
    hData->Add(hWW, -1.);
    hData->Add(hWZ, -1.);
    hData->Add(hZZ, -1.);
  */

  double fitVal, fitError;
  makeFit("sIetaIeta", 0.005, 0.02, hZG, hDY, hData, fitVal, fitError);

  cout << endl << "ZG Fit returned ZG fraction = " << fitVal << " +/- " << fitError << endl;

  double dataInt = hData->Integral();
  double zgInt = hZG->Integral();
  double dyInt = hDY->Integral();
    
  double zgSF = fitVal * dataInt / zgInt;
  double zgSFerror = fitError * dataInt / zgInt;

  double dySF = (1. - fitVal) * dataInt / dyInt;
  double dySFerror = fitError * dataInt / dyInt;

  cout << "-------------------------------------------------------------" << endl;
  cout << "zgSF = " << zgSF << " +/- " << zgSFerror << endl;
  cout << "dySF = " << dySF << " +/- " << dySFerror << endl;
  cout << "-------------------------------------------------------------" << endl << endl;

  hZG->Scale(zgSF);
  hDY->Scale(dySF);

  TCanvas * can = new TCanvas("fit_can", "Plot", 10, 10, 2000, 2000);

  hData->SetLineWidth(3);
  hData->GetXaxis()->SetTitle("#gamma #sigma_{i#etai#eta}");

  hDY->SetFillColor(kGray);
  hDY->Add(hZG);
  hDY->Add(hWW);
  hDY->Add(hWZ);
  hDY->Add(hZZ);
  hDY->Add(hSemiLep);
  hDY->Add(hFullLep);
  hDY->Add(hHadronic);

  hZG->SetFillColor(8);
  hZG->Add(hWW);
  hZG->Add(hWZ);
  hZG->Add(hZZ);
  hZG->Add(hSemiLep);
  hZG->Add(hFullLep);
  hZG->Add(hHadronic);

  hWW->SetFillColor(kRed);
  hWW->Add(hWZ);
  hWW->Add(hZZ);
  hWW->Add(hSemiLep);
  hWW->Add(hFullLep);
  hWW->Add(hHadronic);
  
  hSemiLep->SetFillColor(kBlue);
  hSemiLep->Add(hFullLep);
  hSemiLep->Add(hHadronic);

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(hData, "Data", "LP");
  leg->AddEntry(hDY, "Drell-Yan", "F");
  leg->AddEntry(hZG, "Z#gamma", "F");
  leg->AddEntry(hWW, "Di-boson", "F");
  leg->AddEntry(hSemiLep, "t#bar{t}", "F");

  TPaveText * reqText = new TPaveText(0.45, 0.47, 0.85, 0.57, "NDC");
  reqText->SetFillColor(0);
  reqText->SetFillStyle(0);
  reqText->SetLineColor(0);
  if(runElectrons) reqText->AddText("ee+#gamma");
  else reqText->AddText("#mu#mu+#gamma");
  reqText->Draw("same");

  hData->Draw("axis");
  hDY->Draw("hist same");
  hZG->Draw("hist same");
  hWW->Draw("hist same");
  hSemiLep->Draw("hist same");
  hData->Draw("e1 same");
  hData->Draw("axis same");
  leg->Draw("same");
  reqText->Draw("same");

  can->SaveAs("ZG_sIetaIeta.gif");
  can->SetLogy(true);
  can->SaveAs("ZG_sIetaIeta_log.gif");
  can->SetLogy(false);
  
  TH1D * h_zmass_Data = fillDataHistogram("z_invmass", dataTree, 1,
					  "data_z_invmass", 100, 0., 500.);
  
  TH1D * h_zmass_DY1 = fillMCHistogram("z_invmass", tDY1, 1,
				       "dy1_z_invmass", 100, 0., 500.,
				       scale_dy1);
  
  TH1D * h_zmass_DY2 = fillMCHistogram("z_invmass", tDY2, 1,
				       "dy2_z_invmass", 100, 0., 500.,
				       scale_dy2);
  
  TH1D * h_zmass_DY3 = fillMCHistogram("z_invmass", tDY3, 1,
				     "dy3_z_invmass", 100, 0., 500.,
				     scale_dy3);
  
  TH1D * h_zmass_DY4 = fillMCHistogram("z_invmass", tDY4, 1,
				     "dy4_z_invmass", 100, 0., 500.,
				     scale_dy4);
  
  TH1D * h_zmass_DY = (TH1D*)h_zmass_DY1->Clone("dy_z_invmass");
  h_zmass_DY->Add(h_zmass_DY2);
  h_zmass_DY->Add(h_zmass_DY3);
  h_zmass_DY->Add(h_zmass_DY4);

  h_zmass_DY->Scale(dySF);

  TH1D * h_zmass_ZG = fillMCHistogram("z_invmass", tZG, 1,
				    "zg_z_invmass", 100, 0., 500.,
				    scale_zg);

  h_zmass_ZG->Scale(zgSF);

  TH1D * h_zmass_FullLep = fillMCHistogram("z_invmass", tFullLep, 1,
					 "fullLep_z_invmass", 100, 0., 500.,
					 scale_fullLep);
  
  TH1D * h_zmass_SemiLep = fillMCHistogram("z_invmass", tSemiLep, 1,
					 "semiLep_z_invmass", 100, 0., 500.,
					 scale_semiLep);
  
  TH1D * h_zmass_Hadronic = fillMCHistogram("z_invmass", tHadronic, 1,
					  "hadronic_z_invmass", 100, 0., 500.,
					  scale_hadronic);
  
  TH1D * h_zmass_WW = fillMCHistogram("z_invmass", tWW, 1,
				    "WW_z_invmass", 100, 0., 500.,
				    scale_WW);
  
  TH1D * h_zmass_WZ = fillMCHistogram("z_invmass", tWZ, 1,
				    "WZ_z_invmass", 100, 0., 500.,
				    scale_WZ);
  
  TH1D * h_zmass_ZZ = fillMCHistogram("z_invmass", tZZ, 1,
				    "ZZ_z_invmass", 100, 0., 500.,
				    scale_ZZ);
  
  h_zmass_Data->SetLineWidth(3);
  if(runElectrons) h_zmass_Data->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  else h_zmass_Data->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");

  h_zmass_DY->SetFillColor(kGray);
  h_zmass_DY->Add(h_zmass_ZG);
  h_zmass_DY->Add(h_zmass_WW);
  h_zmass_DY->Add(h_zmass_WZ);
  h_zmass_DY->Add(h_zmass_ZZ);
  h_zmass_DY->Add(h_zmass_SemiLep);
  h_zmass_DY->Add(h_zmass_FullLep);
  h_zmass_DY->Add(h_zmass_Hadronic);

  h_zmass_ZG->SetFillColor(8);
  h_zmass_ZG->Add(h_zmass_WW);
  h_zmass_ZG->Add(h_zmass_WZ);
  h_zmass_ZG->Add(h_zmass_ZZ);
  h_zmass_ZG->Add(h_zmass_SemiLep);
  h_zmass_ZG->Add(h_zmass_FullLep);
  h_zmass_ZG->Add(h_zmass_Hadronic);

  h_zmass_WW->SetFillColor(kRed);
  h_zmass_WW->Add(h_zmass_WZ);
  h_zmass_WW->Add(h_zmass_ZZ);
  h_zmass_WW->Add(h_zmass_SemiLep);
  h_zmass_WW->Add(h_zmass_FullLep);
  h_zmass_WW->Add(h_zmass_Hadronic);
  
  h_zmass_SemiLep->SetFillColor(kBlue);
  h_zmass_SemiLep->Add(h_zmass_FullLep);
  h_zmass_SemiLep->Add(h_zmass_Hadronic);

  h_zmass_Data->Draw("axis");
  h_zmass_DY->Draw("hist same");
  h_zmass_ZG->Draw("hist same");
  h_zmass_WW->Draw("hist same");
  h_zmass_SemiLep->Draw("hist same");
  h_zmass_Data->Draw("e1 same");
  h_zmass_Data->Draw("axis same");
  leg->Draw("same");
  reqText->Draw("same");

  TLine * lowEdge = new TLine(76, h_zmass_Data->GetMinimum(), 76, h_zmass_Data->GetMaximum());
  lowEdge->SetLineColor(kRed);
  lowEdge->SetLineWidth(3);
  lowEdge->Draw("same");

  TLine * highEdge = new TLine(106, h_zmass_Data->GetMinimum(), 106, h_zmass_Data->GetMaximum());
  highEdge->SetLineColor(kRed);
  highEdge->SetLineWidth(3);
  highEdge->Draw("same");

  can->SaveAs("ZG_zmass.gif");
  can->SetLogy(true);
  can->SaveAs("ZG_zmass_log.gif");
  can->SetLogy(false);

  TH1D * h_met_Data = fillDataHistogram("pfMET", dataTree, 2,
					"data_pfMET", 12, 0., 300.);

  TH1D * h_met_DY1 = fillMCHistogram("pfMET", tDY1, 2,
				     "dy1_pfMET", 12, 0., 300.,
				     scale_dy1);
  
  TH1D * h_met_DY2 = fillMCHistogram("pfMET", tDY2, 2,
				     "dy2_pfMET", 12, 0., 300.,
				     scale_dy2);
  
  TH1D * h_met_DY3 = fillMCHistogram("pfMET", tDY3, 2,
				     "dy3_pfMET", 12, 0., 300.,
				     scale_dy3);
  
  TH1D * h_met_DY4 = fillMCHistogram("pfMET", tDY4, 2,
				     "dy4_pfMET", 12, 0., 300.,
				     scale_dy4);
  
  TH1D * h_met_DY = (TH1D*)h_met_DY1->Clone("dy_pfMET");
  h_met_DY->Add(h_met_DY2);
  h_met_DY->Add(h_met_DY3);
  h_met_DY->Add(h_met_DY4);

  h_met_DY->Scale(dySF);

  TH1D * h_met_ZG = fillMCHistogram("pfMET", tZG, 2,
				    "zg_pfMET", 12, 0., 300.,
				    scale_zg);

  h_met_ZG->Scale(zgSF);

  TH1D * h_met_FullLep = fillMCHistogram("pfMET", tFullLep, 2,
					 "fullLep_pfMET", 12, 0., 300.,
					 scale_fullLep);
  
  TH1D * h_met_SemiLep = fillMCHistogram("pfMET", tSemiLep, 2,
					 "semiLep_pfMET", 12, 0., 300.,
					 scale_semiLep);
  
  TH1D * h_met_Hadronic = fillMCHistogram("pfMET", tHadronic, 2,
					  "hadronic_pfMET", 12, 0., 300.,
					  scale_hadronic);
  
  TH1D * h_met_WW = fillMCHistogram("pfMET", tWW, 2,
				    "WW_pfMET", 12, 0., 300.,
				    scale_WW);
  
  TH1D * h_met_WZ = fillMCHistogram("pfMET", tWZ, 2,
				    "WZ_pfMET", 12, 0., 300.,
				    scale_WZ);
  
  TH1D * h_met_ZZ = fillMCHistogram("pfMET", tZZ, 2,
				    "ZZ_pfMET", 12, 0., 300.,
				    scale_ZZ);
  
  h_met_Data->SetLineWidth(3);
  h_met_Data->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");

  h_met_DY->SetFillColor(kGray);
  h_met_DY->Add(h_met_ZG);
  h_met_DY->Add(h_met_WW);
  h_met_DY->Add(h_met_WZ);
  h_met_DY->Add(h_met_ZZ);
  h_met_DY->Add(h_met_SemiLep);
  h_met_DY->Add(h_met_FullLep);
  h_met_DY->Add(h_met_Hadronic);

  h_met_ZG->SetFillColor(8);
  h_met_ZG->Add(h_met_WW);
  h_met_ZG->Add(h_met_WZ);
  h_met_ZG->Add(h_met_ZZ);
  h_met_ZG->Add(h_met_SemiLep);
  h_met_ZG->Add(h_met_FullLep);
  h_met_ZG->Add(h_met_Hadronic);

  h_met_WW->SetFillColor(kRed);
  h_met_WW->Add(h_met_WZ);
  h_met_WW->Add(h_met_ZZ);
  h_met_WW->Add(h_met_SemiLep);
  h_met_WW->Add(h_met_FullLep);
  h_met_WW->Add(h_met_Hadronic);
  
  h_met_SemiLep->SetFillColor(kBlue);
  h_met_SemiLep->Add(h_met_FullLep);
  h_met_SemiLep->Add(h_met_Hadronic);

  reqText->Clear();
  if(runElectrons) reqText->AddText("ee+#gamma#gamma");
  else reqText->AddText("#mu#mu+#gamma#gamma");

  h_met_Data->Draw("axis");
  h_met_DY->Draw("hist same");
  h_met_ZG->Draw("hist same");
  h_met_WW->Draw("hist same");
  h_met_SemiLep->Draw("hist same");
  h_met_Data->Draw("e1 same");
  h_met_Data->Draw("axis same");
  leg->Draw("same");
  reqText->Draw("same");

  can->SaveAs("ZGG_met.gif");
  can->SetLogy(true);
  can->SaveAs("ZGG_met_log.gif");
  can->SetLogy(false);

  TH1D * h_photonEt_Data = fillDataHistogram("leadPhotonEt", dataTree, 2,
					     "data_leadPhotonEt", 6, 0., 300.);

  TH1D * h_photonEt_DY1 = fillMCHistogram("leadPhotonEt", tDY1, 2,
					  "dy1_leadPhotonEt", 6, 0., 300.,
					  scale_dy1);
  
  TH1D * h_photonEt_DY2 = fillMCHistogram("leadPhotonEt", tDY2, 2,
					  "dy2_leadPhotonEt", 6, 0., 300.,
					  scale_dy2);
  
  TH1D * h_photonEt_DY3 = fillMCHistogram("leadPhotonEt", tDY3, 2,
					  "dy3_leadPhotonEt", 6, 0., 300.,
					  scale_dy3);
  
  TH1D * h_photonEt_DY4 = fillMCHistogram("leadPhotonEt", tDY4, 2,
					  "dy4_leadPhotonEt", 6, 0., 300.,
					  scale_dy4);
  
  TH1D * h_photonEt_DY = (TH1D*)h_photonEt_DY1->Clone("dy_leadPhotonEt");
  h_photonEt_DY->Add(h_photonEt_DY2);
  h_photonEt_DY->Add(h_photonEt_DY3);
  h_photonEt_DY->Add(h_photonEt_DY4);

  h_photonEt_DY->Scale(dySF);

  TH1D * h_photonEt_ZG = fillMCHistogram("leadPhotonEt", tZG, 2,
					 "zg_leadPhotonEt", 6, 0., 300.,
					 scale_zg);

  h_photonEt_ZG->Scale(zgSF);

  TH1D * h_photonEt_FullLep = fillMCHistogram("leadPhotonEt", tFullLep, 2,
					      "fullLep_leadPhotonEt", 6, 0., 300.,
					      scale_fullLep);
  
  TH1D * h_photonEt_SemiLep = fillMCHistogram("leadPhotonEt", tSemiLep, 2,
					      "semiLep_leadPhotonEt", 6, 0., 300.,
					      scale_semiLep);
  
  TH1D * h_photonEt_Hadronic = fillMCHistogram("leadPhotonEt", tHadronic, 2,
					       "hadronic_leadPhotonEt", 6, 0., 300.,
					       scale_hadronic);
  
  TH1D * h_photonEt_WW = fillMCHistogram("leadPhotonEt", tWW, 2,
					 "WW_leadPhotonEt", 6, 0., 300.,
					 scale_WW);
  
  TH1D * h_photonEt_WZ = fillMCHistogram("leadPhotonEt", tWZ, 2,
					 "WZ_leadPhotonEt", 6, 0., 300.,
					 scale_WZ);
  
  TH1D * h_photonEt_ZZ = fillMCHistogram("leadPhotonEt", tZZ, 2,
					 "ZZ_leadPhotonEt", 6, 0., 300.,
					 scale_ZZ);

  h_photonEt_Data->SetLineWidth(3);
  h_photonEt_Data->GetXaxis()->SetTitle("E_{T} of leading #gamma");

  h_photonEt_DY->SetFillColor(kGray);
  h_photonEt_DY->Add(h_photonEt_ZG);
  h_photonEt_DY->Add(h_photonEt_WW);
  h_photonEt_DY->Add(h_photonEt_WZ);
  h_photonEt_DY->Add(h_photonEt_ZZ);
  h_photonEt_DY->Add(h_photonEt_SemiLep);
  h_photonEt_DY->Add(h_photonEt_FullLep);
  h_photonEt_DY->Add(h_photonEt_Hadronic);

  h_photonEt_ZG->SetFillColor(8);
  h_photonEt_ZG->Add(h_photonEt_WW);
  h_photonEt_ZG->Add(h_photonEt_WZ);
  h_photonEt_ZG->Add(h_photonEt_ZZ);
  h_photonEt_ZG->Add(h_photonEt_SemiLep);
  h_photonEt_ZG->Add(h_photonEt_FullLep);
  h_photonEt_ZG->Add(h_photonEt_Hadronic);

  h_photonEt_WW->SetFillColor(kRed);
  h_photonEt_WW->Add(h_photonEt_WZ);
  h_photonEt_WW->Add(h_photonEt_ZZ);
  h_photonEt_WW->Add(h_photonEt_SemiLep);
  h_photonEt_WW->Add(h_photonEt_FullLep);
  h_photonEt_WW->Add(h_photonEt_Hadronic);
  
  h_photonEt_SemiLep->SetFillColor(kBlue);
  h_photonEt_SemiLep->Add(h_photonEt_FullLep);
  h_photonEt_SemiLep->Add(h_photonEt_Hadronic);

  h_photonEt_Data->Draw("axis");
  h_photonEt_DY->Draw("hist same");
  h_photonEt_ZG->Draw("hist same");
  h_photonEt_WW->Draw("hist same");
  h_photonEt_SemiLep->Draw("hist same");
  h_photonEt_Data->Draw("e1 same");
  h_photonEt_Data->Draw("axis same");
  leg->Draw("same");
  reqText->Draw("same");

  can->SaveAs("ZGG_photonEt.gif");
  can->SetLogy(true);
  can->SaveAs("ZGG_photonEt_log.gif");
  can->SetLogy(false);

  delete can;

  fDY1->Close();
  fDY2->Close();
  fDY3->Close();
  fDY4->Close();
  fZG->Close();

  fFullLep->Close();
  fSemiLep->Close();
  fHadronic->Close();
  
  fWW->Close();
  fWZ->Close();
  fZZ->Close();

  fData->Close();

}
