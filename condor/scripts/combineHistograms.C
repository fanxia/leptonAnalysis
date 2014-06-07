void combineHistograms(TString datasetname) {

  TFile * puFile = new TFile("pileupReweighting_temp_"+datasetname+".root", "READ");

  TH1D * data = (TH1D*)puFile->Get("pu_data_nonorm_"+datasetname);
  TH1D * data_up = (TH1D*)puFile->Get("pu_data_up_nonorm_"+datasetname);
  TH1D * data_down = (TH1D*)puFile->Get("pu_data_down_nonorm_"+datasetname);

  TH1D * mc = (TH1D*)puFile->Get("pu_mc_nonorm_"+datasetname);

  data->Scale(1./data->Integral());
  data_up->Scale(1./data_up->Integral());
  data_down->Scale(1./data_down->Integral());
  mc->Scale(1./mc->Integral());

  TH1D * weights = (TH1D*)data->Clone("puWeights_"+datasetname);
  weights->Divide(mc);

  TH1D * weights_up = (TH1D*)data_up->Clone("puWeights_up_"+datasetname);
  weights_up->Divide(mc);

  TH1D * weights_down = (TH1D*)data_down->Clone("puWeights_down_"+datasetname);
  weights_down->Divide(mc);

  TFile * outFile = new TFile("pileupReweighting_"+datasetname+".root", "RECREATE");
  outFile->cd();

  weights->Write();
  weights_up->Write();
  weights_down->Write();
  data->Write("pu_data_"+datasetname);
  data_up->Write("pu_data_up_"+datasetname);
  data_down->Write("pu_data_down_"+datasetname);
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

