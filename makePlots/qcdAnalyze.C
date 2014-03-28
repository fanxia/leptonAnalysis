void qcdAnalyze() {

  const Double_t xsec_ttHadronic = 245.8 * 0.457;
  const Double_t xsec_ttSemiLep = 245.8 * 0.438;
  const Double_t xsec_ttFullLep = 245.8 * 0.105;
  
  const Double_t xsec_T_s = 3.79;
  const Double_t xsec_T_t = 56.4;
  const Double_t xsec_T_tW = 11.1;
  
  const Double_t xsec_Tbar_s = 1.76;
  const Double_t xsec_Tbar_t = 30.7;
  const Double_t xsec_Tbar_tW = 11.1;
  
  const Double_t xsec_TTW = 0.232;
  const Double_t xsec_TTZ = 0.2057;
  
  const Double_t xsec_wjets = 36257.2;
  
  const Double_t xsec_dyjets = 3503.71;
  const Double_t xsec_dy1jets = 666.7;
  const Double_t xsec_dy2jets = 215.1;
  const Double_t xsec_dy3jets = 66.07;
  const Double_t xsec_dy4jets = 27.38;
  
  const Double_t xsec_ttgjets = 14.0;
  const Double_t xsec_ttgg = 0.146;
  
  const Double_t xsec_siga = 0.147492;
  const Double_t xsec_sigb = 0.0399591;
  
  TString sigName = "ele_bjj_signalTree";
  TString qcdName = "ele_bjj_eQCDTree";

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * in = new TFile("inputs/SingleElectron.root", "READ");
  TTree * ggTree = (TTree*)in->Get(sigName);
  TTree * qcdTree = (TTree*)in->Get(qcdName);

  TFile * fTTHadronic = new TFile("inputs/signal_contamination_ttJetsHadronic.root", "READ");
  TTree * ttHadronicTree = (TTree*)fTTHadronic->Get(qcdName);
  TH1D * h_nGen = (TH1D*)fTTHadronic->Get("nEvents_ttJetsHadronic");
  double nGen_ttHadronic = h_nGen->Integral();

  TFile * fTTSemiLep = new TFile("inputs/signal_contamination_ttJetsSemiLep.root", "READ");
  TTree * ttSemiLepTree = (TTree*)fTTSemiLep->Get(qcdName);
  h_nGen = (TH1D*)fTTSemiLep->Get("nEvents_ttJetsSemiLep");
  double nGen_ttSemiLep = h_nGen->Integral();

  TFile * fTTFullLep = new TFile("inputs/signal_contamination_ttJetsFullLep.root", "READ");
  TTree * ttFullLepTree = (TTree*)fTTFullLep->Get(qcdName);
  h_nGen = (TH1D*)fTTFullLep->Get("nEvents_ttJetsFullLep");
  double nGen_ttFullLep = h_nGen->Integral();

  TFile * fTBar_s = new TFile("inputs/signal_contamination_TBar_s.root", "READ");
  TTree * tbar_sTree = (TTree*)fTBar_s->Get(qcdName);
  h_nGen = (TH1D*)fTBar_s->Get("nEvents_TBar_s");
  double nGen_tbar_s = h_nGen->Integral();

  TFile * fTBar_t = new TFile("inputs/signal_contamination_TBar_t.root", "READ");
  TTree * tbar_tTree = (TTree*)fTBar_t->Get(qcdName);
  h_nGen = (TH1D*)fTBar_t->Get("nEvents_TBar_t");
  double nGen_tbar_t = h_nGen->Integral();

  TFile * fTBar_tW = new TFile("inputs/signal_contamination_TBar_tW.root", "READ");
  TTree * tbar_tWTree = (TTree*)fTBar_tW->Get(qcdName);
  h_nGen = (TH1D*)fTBar_tW->Get("nEvents_TBar_tW");
  double nGen_tbar_tW = h_nGen->Integral();

  TFile * fT_s = new TFile("inputs/signal_contamination_T_s.root", "READ");
  TTree * t_sTree = (TTree*)fT_s->Get(qcdName);
  h_nGen = (TH1D*)fT_s->Get("nEvents_T_s");
  double nGen_t_s = h_nGen->Integral();

  TFile * fT_t = new TFile("inputs/signal_contamination_T_t.root", "READ");
  TTree * t_tTree = (TTree*)fT_t->Get(qcdName);
  h_nGen = (TH1D*)fT_t->Get("nEvents_T_t");
  double nGen_t_t = h_nGen->Integral();

  TFile * fT_tW = new TFile("inputs/signal_contamination_T_tW.root", "READ");
  TTree * t_tWTree = (TTree*)fT_tW->Get(qcdName);
  h_nGen = (TH1D*)fT_tW->Get("nEvents_T_tW");
  double nGen_t_tW = h_nGen->Integral();

  TFile * fTTW = new TFile("inputs/signal_contamination_TTWJets.root", "READ");
  TTree * ttwjetsTree = (TTree*)fTTW->Get(qcdName);
  h_nGen = (TH1D*)fTTW->Get("nEvents_TTWJets");
  double nGen_ttwjets = h_nGen->Integral();

  TFile * fTTZ = new TFile("inputs/signal_contamination_TTZJets.root", "READ");
  TTree * ttzjetsTree = (TTree*)fTTZ->Get(qcdName);
  h_nGen = (TH1D*)fTTZ->Get("nEvents_TTZJets");
  double nGen_ttzjets = h_nGen->Integral();

  TFile * fWJets = new TFile("inputs/signal_contamination_WJetsToLNu.root", "READ");
  TTree * wjetsTree = (TTree*)fWJets->Get(qcdName);
  h_nGen = (TH1D*)fWJets->Get("nEvents_WJetsToLNu");
  double nGen_wjets = h_nGen->Integral();

  TFile * fDYJets = new TFile("inputs/signal_contamination_dyJetsToLL.root", "READ");
  TTree * dyjetsTree = (TTree*)fDYJets->Get(qcdName);
  h_nGen = (TH1D*)fDYJets->Get("nEvents_dyJetsToLL");
  double nGen_dyjets = h_nGen->Integral();

  TFile * fDY1Jets = new TFile("inputs/signal_contamination_dy1JetsToLL.root", "READ");
  TTree * dy1jetsTree = (TTree*)fDY1Jets->Get(qcdName);
  h_nGen = (TH1D*)fDY1Jets->Get("nEvents_dy1JetsToLL");
  double nGen_dy1jets = h_nGen->Integral();

  TFile * fDY2Jets = new TFile("inputs/signal_contamination_dy2JetsToLL.root", "READ");
  TTree * dy2jetsTree = (TTree*)fDY2Jets->Get(qcdName);
  h_nGen = (TH1D*)fDY2Jets->Get("nEvents_dy2JetsToLL");
  double nGen_dy2jets = h_nGen->Integral();

  TFile * fDY3Jets = new TFile("inputs/signal_contamination_dy3JetsToLL.root", "READ");
  TTree * dy3jetsTree = (TTree*)fDY3Jets->Get(qcdName);
  h_nGen = (TH1D*)fDY3Jets->Get("nEvents_dy3JetsToLL");
  double nGen_dy3jets = h_nGen->Integral();

  TFile * fDY4Jets = new TFile("inputs/signal_contamination_dy4JetsToLL.root", "READ");
  TTree * dy4jetsTree = (TTree*)fDY4Jets->Get(qcdName);
  h_nGen = (TH1D*)fDY4Jets->Get("nEvents_dy4JetsToLL");
  double nGen_dy4jets = h_nGen->Integral();

  TFile * fTTGJets = new TFile("inputs/signal_contamination_ttgjets.root", "READ");
  TTree * ttgjetsTree = (TTree*)fTTGJets->Get(qcdName);
  h_nGen = (TH1D*)fTTGJets->Get("nEvents_ttgjets");
  double nGen_ttgjets = h_nGen->Integral();

  TFile * fTTGG = new TFile("inputs/signal_contamination_ttGG.root", "READ");
  TTree * ttggTree = (TTree*)fTTGG->Get(qcdName);
  h_nGen = (TH1D*)fTTGG->Get("nEvents_ttGG");
  double nGen_ttgg = h_nGen->Integral();

  delete h_nGen;

  Float_t mva, iso, nphotons;
  
  qcdTree->SetBranchAddress("Nphotons", &nphotons);
  ttHadronicTree->SetBranchAddress("Nphotons", &nphotons);
  ttSemiLepTree->SetBranchAddress("Nphotons", &nphotons);
  ttFullLepTree->SetBranchAddress("Nphotons", &nphotons);
  tbar_sTree->SetBranchAddress("Nphotons", &nphotons);
  tbar_tTree->SetBranchAddress("Nphotons", &nphotons);
  tbar_tWTree->SetBranchAddress("Nphotons", &nphotons);
  t_sTree->SetBranchAddress("Nphotons", &nphotons);
  t_tTree->SetBranchAddress("Nphotons", &nphotons);
  t_tWTree->SetBranchAddress("Nphotons", &nphotons);
  ttwjetsTree->SetBranchAddress("Nphotons", &nphotons);
  ttzjetsTree->SetBranchAddress("Nphotons", &nphotons);
  wjetsTree->SetBranchAddress("Nphotons", &nphotons);
  dyjetsTree->SetBranchAddress("Nphotons", &nphotons);
  dy1jetsTree->SetBranchAddress("Nphotons", &nphotons);
  dy2jetsTree->SetBranchAddress("Nphotons", &nphotons);
  dy3jetsTree->SetBranchAddress("Nphotons", &nphotons);
  dy4jetsTree->SetBranchAddress("Nphotons", &nphotons);
  ttgjetsTree->SetBranchAddress("Nphotons", &nphotons);
  ttggTree->SetBranchAddress("Nphotons", &nphotons);

  qcdTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttHadronicTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttSemiLepTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttFullLepTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  tbar_sTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  tbar_tTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  tbar_tWTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  t_sTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  t_tTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  t_tWTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttwjetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttzjetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  wjetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  dyjetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  dy1jetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  dy2jetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  dy3jetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  dy4jetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttgjetsTree->SetBranchAddress("ele_mvaTrigV0", &mva);
  ttggTree->SetBranchAddress("ele_mvaTrigV0", &mva);

  qcdTree->SetBranchAddress("ele_relIso", &iso);
  ttHadronicTree->SetBranchAddress("ele_relIso", &iso);
  ttSemiLepTree->SetBranchAddress("ele_relIso", &iso);
  ttFullLepTree->SetBranchAddress("ele_relIso", &iso);
  tbar_sTree->SetBranchAddress("ele_relIso", &iso);
  tbar_tTree->SetBranchAddress("ele_relIso", &iso);
  tbar_tWTree->SetBranchAddress("ele_relIso", &iso);
  t_sTree->SetBranchAddress("ele_relIso", &iso);
  t_tTree->SetBranchAddress("ele_relIso", &iso);
  t_tWTree->SetBranchAddress("ele_relIso", &iso);
  ttwjetsTree->SetBranchAddress("ele_relIso", &iso);
  ttzjetsTree->SetBranchAddress("ele_relIso", &iso);
  wjetsTree->SetBranchAddress("ele_relIso", &iso);
  dyjetsTree->SetBranchAddress("ele_relIso", &iso);
  dy1jetsTree->SetBranchAddress("ele_relIso", &iso);
  dy2jetsTree->SetBranchAddress("ele_relIso", &iso);
  dy3jetsTree->SetBranchAddress("ele_relIso", &iso);
  dy4jetsTree->SetBranchAddress("ele_relIso", &iso);
  ttgjetsTree->SetBranchAddress("ele_relIso", &iso);
  ttggTree->SetBranchAddress("ele_relIso", &iso);

  TH1D * h_qcd_mva = new TH1D("qcd_mva", "qcd_mva", 40, -1, 1); h_qcd_mva->Sumw2();
  TH1D * h_ttbar_mva = new TH1D("ttbar_mva", "ttbar_mva", 40, -1, 1); h_ttbar_mva->Sumw2();
  TH1D * h_singleTop_mva = new TH1D("singleTop_mva", "singleTop_mva", 40, -1, 1); h_singleTop_mva->Sumw2();
  TH1D * h_ttwz_mva = new TH1D("ttwz_mva", "ttwz_mva", 40, -1, 1); h_ttwz_mva->Sumw2();
  TH1D * h_wjets_mva = new TH1D("wjets_mva", "wjets_mva", 40, -1, 1); h_wjets_mva->Sumw2();
  TH1D * h_dyjets_mva = new TH1D("dyjets_mva", "dyjets_mva", 40, -1, 1); h_dyjets_mva->Sumw2();
  TH1D * h_ttgjets_mva = new TH1D("ttgjets_mva", "ttgjets_mva", 40, -1, 1); h_ttgjets_mva->Sumw2();
  TH1D * h_ttgg_mva = new TH1D("ttgg_mva", "ttgg_mva", 40, -1, 1); h_ttgg_mva->Sumw2();

  TH1D * h_qcd_iso = new TH1D("qcd_iso", "qcd_iso", 100, 0, 10); h_qcd_iso->Sumw2();
  TH1D * h_ttbar_iso = new TH1D("ttbar_iso", "ttbar_iso", 100, 0, 10); h_ttbar_iso->Sumw2();
  TH1D * h_singleTop_iso = new TH1D("singleTop_iso", "singleTop_iso", 100, 0, 10); h_singleTop_iso->Sumw2();
  TH1D * h_ttwz_iso = new TH1D("ttwz_iso", "ttwz_iso", 100, 0, 10); h_ttwz_iso->Sumw2();
  TH1D * h_wjets_iso = new TH1D("wjets_iso", "wjets_iso", 100, 0, 10); h_wjets_iso->Sumw2();
  TH1D * h_dyjets_iso = new TH1D("dyjets_iso", "dyjets_iso", 100, 0, 10); h_dyjets_iso->Sumw2();
  TH1D * h_ttgjets_iso = new TH1D("ttgjets_iso", "ttgjets_iso", 100, 0, 10); h_ttgjets_iso->Sumw2();
  TH1D * h_ttgg_iso = new TH1D("ttgg_iso", "ttgg_iso", 100, 0, 10); h_ttgg_iso->Sumw2();

  for(int i = 0; i < qcdTree->GetEntries(); i++) {
    qcdTree->GetEntry(i);
    if(nphotons == 0) {
      h_qcd_mva->Fill(mva);
      h_qcd_iso->Fill(iso);
    }
  }

  Double_t weight = 19712. * xsec_ttHadronic / nGen_ttHadronic;
  for(int i = 0; i < ttHadronicTree->GetEntries(); i++) {
    ttHadronicTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttbar_mva->Fill(mva, weight);
      h_ttbar_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_ttSemiLep / nGen_ttSemiLep;
  for(int i = 0; i < ttSemiLepTree->GetEntries(); i++) {
    ttSemiLepTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttbar_mva->Fill(mva, weight);
      h_ttbar_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_ttFullLep / nGen_ttFullLep;
  for(int i = 0; i < ttFullLepTree->GetEntries(); i++) {
    ttFullLepTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttbar_mva->Fill(mva, weight);
      h_ttbar_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_Tbar_s / nGen_tbar_s;
  for(int i = 0; i < tbar_sTree->GetEntries(); i++) {
    tbar_sTree->GetEntry(i);
    if(nphotons == 0) {
      h_singleTop_mva->Fill(mva, weight);
      h_singleTop_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_Tbar_t / nGen_tbar_t;
  for(int i = 0; i < tbar_tTree->GetEntries(); i++) {
    tbar_tTree->GetEntry(i);
    if(nphotons == 0) {
      h_singleTop_mva->Fill(mva, weight);
      h_singleTop_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_Tbar_tW / nGen_tbar_tW;
  for(int i = 0; i < tbar_tWTree->GetEntries(); i++) {
    tbar_tWTree->GetEntry(i);
    if(nphotons == 0) {
      h_singleTop_mva->Fill(mva, weight);
      h_singleTop_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_T_s / nGen_t_s;
  for(int i = 0; i < t_sTree->GetEntries(); i++) {
    t_sTree->GetEntry(i);
    if(nphotons == 0) {
      h_singleTop_mva->Fill(mva, weight);
      h_singleTop_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_T_t / nGen_t_t;
  for(int i = 0; i < t_tTree->GetEntries(); i++) {
    t_tTree->GetEntry(i);
    if(nphotons == 0) {
      h_singleTop_mva->Fill(mva, weight);
      h_singleTop_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_T_tW / nGen_t_tW;
  for(int i = 0; i < t_tWTree->GetEntries(); i++) {
    t_tWTree->GetEntry(i);
    if(nphotons == 0) {
      h_singleTop_mva->Fill(mva, weight);
      h_singleTop_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_TTW / nGen_ttwjets;
  for(int i = 0; i < ttwjetsTree->GetEntries(); i++) {
    ttwjetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttwz_mva->Fill(mva, weight);
      h_ttwz_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_TTZ / nGen_ttzjets;
  for(int i = 0; i < ttzjetsTree->GetEntries(); i++) {
    ttzjetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttwz_mva->Fill(mva, weight);
      h_ttwz_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_dy1jets / nGen_dy1jets;
  for(int i = 0; i < dy1jetsTree->GetEntries(); i++) {
    dy1jetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_dyjets_mva->Fill(mva, weight);
      h_dyjets_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_dy2jets / nGen_dy2jets;
  for(int i = 0; i < dy2jetsTree->GetEntries(); i++) {
    dy2jetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_dyjets_mva->Fill(mva, weight);
      h_dyjets_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_dy3jets / nGen_dy3jets;
  for(int i = 0; i < dy3jetsTree->GetEntries(); i++) {
    dy3jetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_dyjets_mva->Fill(mva, weight);
      h_dyjets_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_dy4jets / nGen_dy4jets;
  for(int i = 0; i < dy4jetsTree->GetEntries(); i++) {
    dy4jetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_dyjets_mva->Fill(mva, weight);
      h_dyjets_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_ttgjets / nGen_ttgjets;
  for(int i = 0; i < ttgjetsTree->GetEntries(); i++) {
    ttgjetsTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttgjets_mva->Fill(mva, weight);
      h_ttgjets_iso->Fill(iso, weight);
    }
  }

  weight = 19712. * xsec_ttgg / nGen_ttgg;
  for(int i = 0; i < ttggTree->GetEntries(); i++) {
    ttggTree->GetEntry(i);
    if(nphotons == 0) {
      h_ttgg_mva->Fill(mva, weight);
      h_ttgg_iso->Fill(iso, weight);
    }
  }

  h_ttbar_mva->SetFillColor(kGray);
  h_ttbar_mva->SetMarkerSize(0);
  h_ttbar_mva->SetLineColor(1);
  h_ttbar_mva->Add(h_singleTop_mva);
  h_ttbar_mva->Add(h_ttwz_mva);
  h_ttbar_mva->Add(h_wjets_mva);
  h_ttbar_mva->Add(h_dyjets_mva);
  h_ttbar_mva->Add(h_ttgjets_mva);
  h_ttbar_mva->Add(h_ttgg_mva);

  h_wjets_mva->SetFillColor(kOrange-3);
  h_wjets_mva->SetMarkerSize(0);
  h_wjets_mva->SetLineColor(1);
  h_wjets_mva->Add(h_singleTop_mva);
  h_wjets_mva->Add(h_ttwz_mva);
  h_ttbar_mva->Add(h_dyjets_mva);
  h_ttbar_mva->Add(h_ttgjets_mva);
  h_ttbar_mva->Add(h_ttgg_mva);

  h_dyjets_mva->SetFillColor(kYellow);
  h_dyjets_mva->SetMarkerSize(0);
  h_dyjets_mva->SetLineColor(1);
  h_dyjets_mva->Add(h_singleTop_mva);
  h_dyjets_mva->Add(h_ttwz_mva);
  h_dyjets_mva->Add(h_ttgjets_mva);
  h_dyjets_mva->Add(h_ttgg_mva);

  h_singleTop_mva->SetFillColor(kRed);
  h_singleTop_mva->SetMarkerSize(0);
  h_singleTop_mva->SetLineColor(1);
  h_singleTop_mva->Add(h_ttwz_mva);
  h_singleTop_mva->Add(h_ttgjets_mva);
  h_singleTop_mva->Add(h_ttgg_mva);

  h_ttwz_mva->SetFillColor(kAzure-2);
  h_ttwz_mva->SetMarkerSize(0);
  h_ttwz_mva->SetLineColor(1);
  h_ttwz_mva->Add(h_ttgjets_mva);
  h_ttwz_mva->Add(h_ttgg_mva);

  h_ttgjets_mva->SetFillColor(8);
  h_ttgjets_mva->SetMarkerSize(0);
  h_ttgjets_mva->SetLineColor(1);
  h_ttgjets_mva->Add(h_ttgg_mva);

  h_ttgg_mva->SetFillColor(kCyan+3);
  h_ttgg_mva->SetMarkerSize(0);
  h_ttgg_mva->SetLineColor(1);

  h_ttbar_iso->SetFillColor(kGray);
  h_ttbar_iso->SetMarkerSize(0);
  h_ttbar_iso->SetLineColor(1);
  h_ttbar_iso->Add(h_singleTop_iso);
  h_ttbar_iso->Add(h_ttwz_iso);
  h_ttbar_iso->Add(h_wjets_iso);
  h_ttbar_iso->Add(h_dyjets_iso);
  h_ttbar_iso->Add(h_ttgjets_iso);
  h_ttbar_iso->Add(h_ttgg_iso);

  h_wjets_iso->SetFillColor(kOrange-3);
  h_wjets_iso->SetMarkerSize(0);
  h_wjets_iso->SetLineColor(1);
  h_wjets_iso->Add(h_singleTop_iso);
  h_wjets_iso->Add(h_ttwz_iso);
  h_ttbar_iso->Add(h_dyjets_iso);
  h_ttbar_iso->Add(h_ttgjets_iso);
  h_ttbar_iso->Add(h_ttgg_iso);

  h_dyjets_iso->SetFillColor(kYellow);
  h_dyjets_iso->SetMarkerSize(0);
  h_dyjets_iso->SetLineColor(1);
  h_dyjets_iso->Add(h_singleTop_iso);
  h_dyjets_iso->Add(h_ttwz_iso);
  h_dyjets_iso->Add(h_ttgjets_iso);
  h_dyjets_iso->Add(h_ttgg_iso);

  h_singleTop_iso->SetFillColor(kRed);
  h_singleTop_iso->SetMarkerSize(0);
  h_singleTop_iso->SetLineColor(1);
  h_singleTop_iso->Add(h_ttwz_iso);
  h_singleTop_iso->Add(h_ttgjets_iso);
  h_singleTop_iso->Add(h_ttgg_iso);

  h_ttwz_iso->SetFillColor(kAzure-2);
  h_ttwz_iso->SetMarkerSize(0);
  h_ttwz_iso->SetLineColor(1);
  h_ttwz_iso->Add(h_ttgjets_iso);
  h_ttwz_iso->Add(h_ttgg_iso);

  h_ttgjets_iso->SetFillColor(8);
  h_ttgjets_iso->SetMarkerSize(0);
  h_ttgjets_iso->SetLineColor(1);
  h_ttgjets_iso->Add(h_ttgg_iso);

  h_ttgg_iso->SetFillColor(kCyan+3);
  h_ttgg_iso->SetMarkerSize(0);
  h_ttgg_iso->SetLineColor(1);

  TCanvas * canv = new TCanvas("canv", "plot", 10, 10, 2000, 2000);
  canv->SetLogy(true);

  h_qcd_mva->Draw("e1");
  h_ttbar_mva->Draw("hist same");
  h_wjets_mva->Draw("hist same");
  h_dyjets_mva->Draw("hist same");
  h_singleTop_mva->Draw("hist same");
  h_ttwz_mva->Draw("hist same");
  h_ttgjets_mva->Draw("hist same");
  h_ttgg_mva->Draw("hist same");

  canv->SaveAs("qcd_mva.pdf");

  h_qcd_iso->Draw("e1");
  h_ttbar_iso->Draw("hist same");
  h_wjets_iso->Draw("hist same");
  h_dyjets_iso->Draw("hist same");
  h_singleTop_iso->Draw("hist same");
  h_ttwz_iso->Draw("hist same");
  h_ttgjets_iso->Draw("hist same");
  h_ttgg_iso->Draw("hist same");

  canv->SaveAs("qcd_iso.pdf");


}
  
  
