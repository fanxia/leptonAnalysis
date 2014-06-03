void createMuonSF() {

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs#22Jan2013_ReReco_of_2012_data_re

  TFile * fIso = new TFile("MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "READ");
  TFile * fID = new TFile("MuonEfficiencies_Run2012ReReco_53X.root", "READ");
  TFile * fTrigger = new TFile("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root", "READ");

  const int nBinsPt = 11;
  const int nBinsEta = 3;

  Double_t ptBins[nBinsPt + 1] = {0, 10, 20, 25, 30, 35, 40, 50, 60, 90, 140, 300};
  Double_t etaBins[nBinsEta + 1] = {0, 0.9, 1.2, 2.1};

  Double_t x, y;

  TFile * fOut = new TFile("muon_scaleFactors_8TeV_53x.root", "RECREATE");
  fOut->cd();

  TH2D * h_muTight_iso = new TH2D("muTight_SF_iso", "muTight_SF_iso", nBinsPt, ptBins, nBinsEta, etaBins); h_muTight_iso->Sumw2();
  TH2D * h_muTight_id = new TH2D("muTight_SF_id", "muTight_SF_id", nBinsPt, ptBins, nBinsEta, etaBins); h_muTight_id->Sumw2();
  TH2D * h_mu_trigger = new TH2D("mu_SF_trigger", "mu_SF_trigger", nBinsPt, ptBins, nBinsEta, etaBins); h_mu_trigger->Sumw2();
  TH2D * h_muTight_iso_id = new TH2D("muTight_SF_iso_id", "muTight_SF_iso_id", nBinsPt, ptBins, nBinsEta, etaBins); h_muTight_iso_id->Sumw2();
  TH2D * h_muTight_full = new TH2D("muTight_SF_full", "muTight_SF_full", nBinsPt, ptBins, nBinsEta, etaBins); h_muTight_full->Sumw2();

  TH2D * h_muLoose_iso = new TH2D("muLoose_SF_iso", "muLoose_SF_iso", nBinsPt, ptBins, nBinsEta, etaBins); h_muLoose_iso->Sumw2();
  TH2D * h_muLoose_id = new TH2D("muLoose_SF_id", "muLoose_SF_id", nBinsPt, ptBins, nBinsEta, etaBins); h_muLoose_id->Sumw2();
  TH2D * h_muLoose_full = new TH2D("muLoose_SF_full", "muLoose_SF_full", nBinsPt, ptBins, nBinsEta, etaBins); h_muLoose_full->Sumw2();
  
  TGraphAsymmErrors * gr_id_tight_0to09 = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Tight_pt_abseta<0.9");
  TGraphAsymmErrors * gr_id_tight_09to12 = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
  TGraphAsymmErrors * gr_id_tight_12to21 = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");

  for(int i = 0; i < gr_id_tight_0to09->GetN(); i++) {
    gr_id_tight_0to09->GetPoint(i, x, y);
    h_muTight_id->SetBinContent(h_muTight_id->FindBin(x, 0.4), y);
    h_muTight_id->SetBinError(h_muTight_id->FindBin(x, 0.4), gr_id_tight_0to09->GetErrorY(i));
  }

  for(int i = 0; i < gr_id_tight_09to12->GetN(); i++) {
    gr_id_tight_09to12->GetPoint(i, x, y);
    h_muTight_id->SetBinContent(h_muTight_id->FindBin(x, 1), y);
    h_muTight_id->SetBinError(h_muTight_id->FindBin(x, 1), gr_id_tight_09to12->GetErrorY(i));
  }

  for(int i = 0; i < gr_id_tight_12to21->GetN(); i++) {
    gr_id_tight_12to21->GetPoint(i, x, y);
    h_muTight_id->SetBinContent(h_muTight_id->FindBin(x, 2), y);
    h_muTight_id->SetBinError(h_muTight_id->FindBin(x, 2), gr_id_tight_12to21->GetErrorY(i));
  }

  TGraphAsymmErrors * gr_id_loose_0to09 = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Loose_pt_abseta<0.9");
  TGraphAsymmErrors * gr_id_loose_09to12 = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Loose_pt_abseta0.9-1.2");
  TGraphAsymmErrors * gr_id_loose_12to21 = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Loose_pt_abseta1.2-2.1");

  for(int i = 0; i < gr_id_loose_0to09->GetN(); i++) {
    gr_id_loose_0to09->GetPoint(i, x, y);
    h_muLoose_id->SetBinContent(h_muLoose_id->FindBin(x, 0.4), y);
    h_muLoose_id->SetBinError(h_muLoose_id->FindBin(x, 0.4), gr_id_loose_0to09->GetErrorY(i));
  }

  for(int i = 0; i < gr_id_loose_09to12->GetN(); i++) {
    gr_id_loose_09to12->GetPoint(i, x, y);
    h_muLoose_id->SetBinContent(h_muLoose_id->FindBin(x, 1), y);
    h_muLoose_id->SetBinError(h_muLoose_id->FindBin(x, 1), gr_id_loose_09to12->GetErrorY(i));
  }

  for(int i = 0; i < gr_id_loose_12to21->GetN(); i++) {
    gr_id_loose_12to21->GetPoint(i, x, y);
    h_muLoose_id->SetBinContent(h_muLoose_id->FindBin(x, 2), y);
    h_muLoose_id->SetBinError(h_muLoose_id->FindBin(x, 2), gr_id_loose_12to21->GetErrorY(i));
  }

  TGraphAsymmErrors * gr_iso_tight_0to09 = (TGraphAsymmErrors*)fIso->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9");
  TGraphAsymmErrors * gr_iso_tight_09to12 = (TGraphAsymmErrors*)fIso->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2");
  TGraphAsymmErrors * gr_iso_tight_12to21 = (TGraphAsymmErrors*)fIso->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1");

  for(int i = 0; i < gr_iso_tight_0to09->GetN(); i++) {
    gr_iso_tight_0to09->GetPoint(i, x, y);
    h_muTight_iso->SetBinContent(h_muTight_iso->FindBin(x, 0.4), y);
    h_muTight_iso->SetBinError(h_muTight_iso->FindBin(x, 0.4), gr_iso_tight_0to09->GetErrorY(i));
  }

  for(int i = 0; i < gr_iso_tight_09to12->GetN(); i++) {
    gr_iso_tight_09to12->GetPoint(i, x, y);
    h_muTight_iso->SetBinContent(h_muTight_iso->FindBin(x, 1), y);
    h_muTight_iso->SetBinError(h_muTight_iso->FindBin(x, 1), gr_iso_tight_09to12->GetErrorY(i));
  }

  for(int i = 0; i < gr_iso_tight_12to21->GetN(); i++) {
    gr_iso_tight_12to21->GetPoint(i, x, y);
    h_muTight_iso->SetBinContent(h_muTight_iso->FindBin(x, 2), y);
    h_muTight_iso->SetBinError(h_muTight_iso->FindBin(x, 2), gr_iso_tight_12to21->GetErrorY(i));
  }

  TGraphAsymmErrors * gr_iso_loose_0to09 = (TGraphAsymmErrors*)fIso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta<0.9");
  TGraphAsymmErrors * gr_iso_loose_09to12 = (TGraphAsymmErrors*)fIso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta0.9-1.2");
  TGraphAsymmErrors * gr_iso_loose_12to21 = (TGraphAsymmErrors*)fIso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Loose_pt_abseta1.2-2.1");

  for(int i = 0; i < gr_iso_loose_0to09->GetN(); i++) {
    gr_iso_loose_0to09->GetPoint(i, x, y);
    h_muLoose_iso->SetBinContent(h_muLoose_iso->FindBin(x, 0.4), y);
    h_muLoose_iso->SetBinError(h_muLoose_iso->FindBin(x, 0.4), gr_iso_loose_0to09->GetErrorY(i));
  }

  for(int i = 0; i < gr_iso_loose_09to12->GetN(); i++) {
    gr_iso_loose_09to12->GetPoint(i, x, y);
    h_muLoose_iso->SetBinContent(h_muLoose_iso->FindBin(x, 1), y);
    h_muLoose_iso->SetBinError(h_muLoose_iso->FindBin(x, 1), gr_iso_loose_09to12->GetErrorY(i));
  }

  for(int i = 0; i < gr_iso_loose_12to21->GetN(); i++) {
    gr_iso_loose_12to21->GetPoint(i, x, y);
    h_muLoose_iso->SetBinContent(h_muLoose_iso->FindBin(x, 2), y);
    h_muLoose_iso->SetBinError(h_muLoose_iso->FindBin(x, 2), gr_iso_loose_12to21->GetErrorY(i));
  }

  TGraphAsymmErrors * gr_trigger_tight_0to09 = (TGraphAsymmErrors*)fTrigger->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
  TGraphAsymmErrors * gr_trigger_tight_09to12 = (TGraphAsymmErrors*)fTrigger->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
  TGraphAsymmErrors * gr_trigger_tight_12to21 = (TGraphAsymmErrors*)fTrigger->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");
  
  for(int i = 0; i < gr_trigger_tight_0to09->GetN(); i++) {
    gr_trigger_tight_0to09->GetPoint(i, x, y);
    h_mu_trigger->SetBinContent(h_mu_trigger->FindBin(x, 0.4), y);
    h_mu_trigger->SetBinError(h_mu_trigger->FindBin(x, 0.4), gr_trigger_tight_0to09->GetErrorY(i));
  }

  for(int i = 0; i < gr_trigger_tight_09to12->GetN(); i++) {
    gr_trigger_tight_09to12->GetPoint(i, x, y);
    h_mu_trigger->SetBinContent(h_mu_trigger->FindBin(x, 1), y);
    h_mu_trigger->SetBinError(h_mu_trigger->FindBin(x, 1), gr_trigger_tight_09to12->GetErrorY(i));
  }

  for(int i = 0; i < gr_trigger_tight_12to21->GetN(); i++) {
    gr_trigger_tight_12to21->GetPoint(i, x, y);
    h_mu_trigger->SetBinContent(h_mu_trigger->FindBin(x, 2), y);
    h_mu_trigger->SetBinError(h_mu_trigger->FindBin(x, 2), gr_trigger_tight_12to21->GetErrorY(i));
  }

  Double_t ptBins[nBinsPt + 1] = {0, 10, 20, 25, 30, 35, 40, 50, 60, 90, 140, 300};
  Double_t etaBins[nBinsEta + 1] = {0, 0.9, 1.2, 2.1};

  for(int i = 0; i < nBinsPt; i++) {
    for(int j = 0; j < nBinsEta; j++) {

      Double_t val_iso = h_muTight_iso->GetBinContent(i+1, j+1);
      Double_t err_iso = h_muTight_iso->GetBinError(i+1, j+1);
      if(val_iso == 0) {
	val_iso = 1;
	err_iso = 0;
      }
      
      Double_t val_id = h_muTight_id->GetBinContent(i+1, j+1);
      Double_t err_id = h_muTight_id->GetBinError(i+1, j+1);
      if(val_id == 0) {
	val_id = 1;
	err_id = 0;
      }
      
      Double_t val_trig = h_mu_trigger->GetBinContent(i+1, j+1);
      Double_t err_trig = h_mu_trigger->GetBinError(i+1, j+1);
      if(val_trig == 0) {
	val_trig = 1;
	err_trig = 0;
      }
      
      Double_t err_iso_id = sqrt(err_iso*err_iso + err_id*err_id);
      Double_t err_full = sqrt(err_iso*err_iso + err_id*err_id + err_trig*err_trig);
      
      if(i > 3) {
	h_muTight_iso_id->SetBinContent(i+1, j+1, val_iso * val_id);
	h_muTight_iso_id->SetBinError(i+1, j+1, err_iso_id);
	
	h_muTight_full->SetBinContent(i+1, j+1, val_iso * val_id * val_trig);
	h_muTight_full->SetBinError(i+1, j+1, err_full);
      }
    }
  }

  for(int i = 0; i < nBinsPt; i++) {
    for(int j = 0; j < nBinsEta; j++) {

      Double_t val_iso = h_muLoose_iso->GetBinContent(i+1, j+1);
      Double_t err_iso = h_muLoose_iso->GetBinError(i+1, j+1);
      if(val_iso == 0) {
	val_iso = 1;
	err_iso = 0;
      }

      Double_t val_id = h_muLoose_id->GetBinContent(i+1, j+1);
      Double_t err_id = h_muLoose_id->GetBinError(i+1, j+1);
      if(val_id == 0) {
	val_id = 1;
	err_id = 0;
      }

      Double_t err_full = sqrt(err_iso*err_iso + err_id*err_id);

      if(i > 3) {
	h_muLoose_full->SetBinContent(i+1, j+1, val_iso * val_id);
	h_muLoose_full->SetBinError(i+1, j+1, err_full);
      }
    }
  }

  h_muTight_iso_id->GetZaxis()->SetRangeUser(0.9, 1.1);
  h_muTight_full->GetZaxis()->SetRangeUser(0.9, 1.1);
  h_muLoose_full->GetZaxis()->SetRangeUser(0.9, 1.1);

  fOut->Write();
  fOut->Close();

  fIso->Close();
  fID->Close();
  fTrigger->Close();

}
