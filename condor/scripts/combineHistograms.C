void combineHistograms(TString datasetname) {

  TFile * puFile = new TFile("pileupReweighting_temp_"+datasetname+".root", "READ");

  TH1D * data = (TH1D*)puFile->Get("pu_data_nonorm_"+datasetname);
  TH1D * mc = (TH1D*)puFile->Get("pu_mc_nonorm_"+datasetname);

  data->Scale(1./data->Integral());
  mc->Scale(1./mc->Integral());

  TH1D * weights = (TH1D*)data->Clone("puWeights_"+datasetname);
  weights->Divide(mc);

  TFile * outFile = new TFile("pileupReweighting_"+datasetname+".root", "RECREATE");
  outFile->cd();

  weights->Write();
  data->Write("pu_data_"+datasetname);
  mc->Write("pu_mc_"+datasetname);
  outFile->Close();

  puFile->Close();


  TFile * btagFile = new TFile("btagEfficiency_temp_"+datasetname+".root", "READ");

  TH1F * bjets = (TH1F*)btagFile->Get("bjets_"+datasetname);
  TH1F * btags = (TH1F*)btagFile->Get("btags_"+datasetname);
  TH1F * bEff = (TH1F*)btags->Clone("bEff_"+datasetname);
  bEff->Divide(bjets);

  TH1F * cjets = (TH1F*)btagFile->Get("cjets_"+datasetname);
  TH1F * ctags = (TH1F*)btagFile->Get("ctags_"+datasetname);
  TH1F * cEff = (TH1F*)ctags->Clone("cEff_"+datasetname);
  cEff->Divide(cjets);

  TH1F * ljets = (TH1F*)btagFile->Get("ljets_"+datasetname);
  TH1F * ltags = (TH1F*)btagFile->Get("ltags_"+datasetname);
  TH1F * lEff = (TH1F*)ltags->Clone("lEff_"+datasetname);
  lEff->Divide(ljets);
  
  TFile * btagOutFile = new TFile("btagEfficiency_"+datasetname+".root", "RECREATE");
  btagOutFile->cd();

  bjets->Write();
  btags->Write();
  bEff->Write();

  cjets->Write();
  ctags->Write();
  cEff->Write();

  ljets->Write();
  ltags->Write();
  lEff->Write();

  btagOutFile->Close();

  btagFile->Close();

}

