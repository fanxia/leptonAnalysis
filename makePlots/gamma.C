#include <vector>

using namespace std;
using namespace RooFit;

void makeFit(TString varname, double varmin, double varmax, TH1D * signalHist, TH1D * backgroundHist, TH1D * dataHist, TString plotName, double& value, double& error) {

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
  
  TCanvas * can = new TCanvas("fit_can", "Plot", 10, 10, 2000, 2000);
  can->SetLogy(true);

  TH1D * h_sig = (TH1D*)signalHist->Clone("h_sig");
  TH1D * h_bkg = (TH1D*)backgroundHist->Clone("h_bkg");

  TH1D * h_sig_unscaled = (TH1D*)signalHist->Clone("h_sig_unscaled");
  TH1D * h_bkg_unscaled = (TH1D*)backgroundHist->Clone("h_bkg_unscaled");

  dataHist->Draw("e1");

  h_sig_unscaled->SetLineColor(kRed);
  h_sig_unscaled->SetLineStyle(2);
  h_sig_unscaled->Draw("hist same");

  h_bkg_unscaled->SetLineColor(kBlue);
  h_bkg_unscaled->SetLineStyle(2);
  h_bkg_unscaled->Draw("hist same");

  h_sig->Scale(value * dataHist->Integral() / h_sig->Integral());
  h_bkg->Scale((1.-value) * dataHist->Integral() / h_bkg->Integral());

  TH1D * h_sum = (TH1D*)h_sig->Clone("h_sum");
  h_sum->Add(h_bkg);
  h_sum->SetLineWidth(3);

  h_sig->SetLineColor(kRed);
  h_sig->SetLineWidth(3);

  h_bkg->SetLineColor(kBlue);
  h_bkg->SetLineWidth(3);
  
  h_sum->Draw("hist same");
  h_sig->Draw("hist same");
  h_bkg->Draw("hist same");
  
  can->SaveAs(plotName);

  delete can;

  return;
}

double GetLeptonSF(TH2D * muSF, TH2D * eleSF, TH2D * eleTriggerSF,
		   Float_t lepton_pt, Float_t lepton_eta, TString channel) {

  Float_t pt, eta, val;

  if(channel.Contains("ele")) {
    pt = min(lepton_pt, (float)199.);
    pt = max(pt, (float)15.);
    eta = min(fabs(lepton_eta), (double)2.39);

    val = eleSF->GetBinContent(eleSF->FindBin(eta, pt));
    val *= eleTriggerSF->GetBinContent(eleTriggerSF->FindBin(eta, pt));
  }
  else if(channel.Contains("muon")) {
    pt = min(lepton_pt, (float)499.);
    pt = max(pt, (float)10.);
    eta = min(fabs(lepton_eta), (double)2.09);

    val = muSF->GetBinContent(muSF->FindBin(pt, eta));
  }
  else val = 0.;

  return val;
}

double GetPhotonSF(TH2D * idSF, TH2D * vetoSF,
		   Float_t photon_et, Float_t photon_eta) {

  Float_t et, eta, val;

  et = min(photon_et, (float)999.);
  et = max(et, (float)15.);
  eta = min(fabs(photon_eta), (double)1.44441);

  val = idSF->GetBinContent(idSF->FindBin(et, eta));
  val *= vetoSF->GetBinContent(vetoSF->FindBin(et, eta));

  return val;
}

TH1D * FillMCHistogram(TString scanName, TString channel, TString version,
		       TH1D * blankHist,
		       vector<TFile*>& files,
		       TH2D * muSF, TH2D * eleSF, TH2D * eleTriggerSF,
		       TH2D * photonIdSF, TH2D * photonVetoSF,
		       Double_t xsec,
		       bool removeTTA, bool reweightTop,
		       double metCut, double& nLowMET) {

  files.push_back(new TFile("inputs/signal_contamination_"+scanName+".root", "READ"));

  TTree * tree = (TTree*)files.back()->Get(channel+"_"+version+"Tree");

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr;
  Float_t overlaps_ttA;
  Float_t topPtReweighting;

  tree->SetBranchAddress("pileupWeight", &puWeight);
  tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree->SetBranchAddress("btagWeight", &btagWeight);
  tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  if(reweightTop) tree->SetBranchAddress("TopPtReweighting", &topPtReweighting);
  if(removeTTA) tree->SetBranchAddress("overlaps_ttA", &overlaps_ttA);

  Float_t var;
  
  if(version == "noSigmaIetaIeta") tree->SetBranchAddress("leadSigmaIetaIeta", &var);
  else if(version == "noChHadIsoTree") tree->SetBranchAddress("leadChargedHadronIso", &var);

  Float_t nphotons, met;
  tree->SetBranchAddress("Nphotons", &nphotons);
  tree->SetBranchAddress("pfMET", &met);
  
  Float_t lepton_pt, lepton_eta;
  if(channel.Contains("ele")) {
    tree->SetBranchAddress("ele_pt", &lepton_pt);
    tree->SetBranchAddress("ele_eta", &lepton_eta);
  }
  else if(channel.Contains("muon")) {
    tree->SetBranchAddress("muon_pt", &lepton_pt);
    tree->SetBranchAddress("muon_eta", &lepton_eta);
  }

  Float_t photon_et, photon_eta;
  tree->SetBranchAddress("leadPhotonEt", &photon_et);
  tree->SetBranchAddress("leadPhotonEta", &photon_eta);

  TH1D * h = (TH1D*)blankHist->Clone(scanName+"_"+version);

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(nphotons != 1) continue;
    if(removeTTA && overlaps_ttA > 0.001) continue;
    if(metCut > 0. && met >= metCut) continue;

    if(btagWeight != btagWeight) continue;
    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    if(reweightTop && topPtReweighting < 0) topPtReweighting = 1.;

    double leptonSF = GetLeptonSF(muSF, eleSF, eleTriggerSF,
				  lepton_pt, lepton_eta, channel);

    double photonSF = GetPhotonSF(photonIdSF, photonVetoSF,
				  photon_et, photon_eta);

    Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
    Float_t addError2_puOnly = btagWeight*btagWeight*puWeightErr*puWeightErr;

    double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
    if(reweightTop) totalWeight *= topPtReweighting;
    
    Float_t oldError = h->GetBinError(h->FindBin(var));
    Float_t newerror = sqrt(oldError*oldError + addError2);
    h->Fill(var, totalWeight);
    h->SetBinError(h->FindBin(var), newerror);

    if(met < 20.) nLowMET += totalWeight;
  }

  TH1D * h_nGen = files.back()->Get("nEvents_"+scanName);
  h->Scale(19712. * xsec / h_nGen->Integral());

  return h;

}

TH1D * FillDataHistogram(TFile * file, TString channel, TString version,
			 TH1D * blankHist,
			 double metCut,
			 double& nLowMET) {

  TTree * tree = (TTree*)file->Get(channel+"_"+version+"Tree");

  Float_t var;
  
  if(version == "noSigmaIetaIeta") tree->SetBranchAddress("leadSigmaIetaIeta", &var);
  else if(version == "noChHadIsoTree") tree->SetBranchAddress("leadChargedHadronIso", &var);

  Float_t nphotons, met;
  tree->SetBranchAddress("Nphotons", &nphotons);
  tree->SetBranchAddress("pfMET", &met);
  
  TH1D * h = (TH1D*)blankHist->Clone("data_"+channel+"_"+version);

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(nphotons != 1) continue;
    if(metCut > 0. && met >= metCut) continue;

    if(met < 20.) nLowMET += 1.;

    h->Fill(var);
  }

  return h;

}

void gamma() {

  double metCut = 50.;

  TString channel = "ele_bjj";
  TString version = "noSigmaIetaIeta";

  vector<TFile*> mcFiles;

  TFile * fLeptonSF = new TFile("../data/lepton_SF_8TeV_53x_baseline.root", "READ");
  TFile * fPhotonSF = new TFile("../data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root", "READ");

  TH2D * sf_electron = (TH2D*)fLeptonSF->Get("TightEleIdIsoSF");
  TH2D * sf_SingleElectronTrigger = (TH2D*)fLeptonSF->Get("TightEleTriggerSF");

  TH2D * sf_muon = (TH2D*)fLeptonSF->Get("mu_pt_eta_full_id_iso_hlt_8TeV");

  TH2D * sf_photon_id = (TH2D*)fPhotonSF->Get("PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
  TH2D * sf_photon_veto = (TH2D*)fPhotonSF->Get("PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");

  TH1D * dummy = new TH1D("dummy", "dummy", 40, 0., 0.04);

  TFile * fData_ele = new TFile("inputs/SingleElectron.root", "READ");
  TFile * fData_muon = new TFile("inputs/SingleMu.root", "READ");

  TH1D * data;
  TH1D * qcd;

  double nLowMET_data = 0.;
  double nLowMET_qcd = 0.;

  if(channel.Contains("ele")) {
    data = (TH1D*)FillDataHistogram(fData_ele, channel, version, dummy, metCut, nLowMET_data);
    qcd = (TH1D*)FillDataHistogram(fData_ele, "ele_jjj_veto", "eQCD", dummy, metCut, nLowMET_qcd);
  }
  else {
    data = (TH1D*)FillDataHistogram(fData_muon, channel, version, dummy, metCut, nLowMET_data);
    qcd = (TH1D*)FillDataHistogram(fData_muon, "muon_jjj_veto", "muQCD", dummy, metCut, nLowMET_qcd);
  }

  double nLowMET_mc = 0.;

  data->Add(FillMCHistogram("W1JetsToLNu", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    12234.4 * 3 * 6662. / 37509., false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("W2JetsToLNu", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    12234.4 * 3 * 2159. / 37509., false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("W3JetsToLNu", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    12234.4 * 3 * 640. / 37509., false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("W4JetsToLNu", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    12234.4 * 3 * 264. / 37509., false, false, metCut, nLowMET_mc),
	    -1.);

  data->Add(FillMCHistogram("dy1JetsToLL", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    666.7 * 1177.3 * 3 / 3503.71, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("dy2JetsToLL", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    215.1 * 1177.3 * 3 / 3503.71, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("dy3JetsToLL", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    66.07 * 1177.3 * 3 / 3503.71, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("dy4JetsToLL", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    27.38 * 1177.3 * 3 / 3503.71, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("TBar_s", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    1.76, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("TBar_t", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    30.7, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("TBar_tW", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    11.1, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("T_s", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    3.79, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("T_t", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    56.4, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("T_tW", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    11.1, false, false, metCut, nLowMET_mc),
	    -1.);

  data->Add(FillMCHistogram("WW", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    57.1097, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("WZ", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    32.3161, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("ZZ", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    8.25561, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("TTWJets", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    0.232, false, false, metCut, nLowMET_mc),
	    -1.);
  data->Add(FillMCHistogram("TTZJets", channel, version, dummy, mcFiles, 
			    sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
			    0.2057, false, false, metCut, nLowMET_mc),
	    -1.);


  qcd->Scale((nLowMET_data - nLowMET_mc) / nLowMET_qcd);
  data->Add(qcd, -1.);

  TH1D * ttgamma = (TH1D*)FillMCHistogram("ttA_2to5", channel, version, dummy, mcFiles, 
					  sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 
					  .9081 * 2, false, true, metCut, nLowMET_mc);

  TH1D * ttjets = (TH1D*)FillMCHistogram("ttJetsSemiLep", channel, version, dummy, mcFiles,
					  sf_muon, sf_electron, sf_SingleElectronTrigger,
					  sf_photon_id, sf_photon_veto,
					  245.8 * 0.438, true, true, metCut, nLowMET_mc);
  ttjets->Add(FillMCHistogram("ttJetsFullLep", channel, version, dummy, mcFiles,
			      sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto,
			      245.8 * 0.105, true, true, metCut, nLowMET_mc));
  ttjets->Add(FillMCHistogram("ttJetsHadronic", channel, version, dummy, mcFiles,
			      sf_muon, sf_electron, sf_SingleElectronTrigger, sf_photon_id, sf_photon_veto, 245.8 * 0.457, true, true, metCut, nLowMET_mc));


  double fitVal, fitError;

  makeFit("sigmaIetaIeta", 0.005, 0.025, ttgamma, ttjets, data, "sIetaIeta_fit_"+channel+".pdf", fitVal, fitError);

  cout << endl << "sIetaIeta Fit returned ttgamma fraction = " << fitVal << " +/- " << fitError << endl;

  double dataInt = data->Integral();
  double ttjetsInt = ttjets->Integral();
  double ttgammaInt = ttgamma->Integral();
    
  double ttgammaSF = fitVal * dataInt / ttgammaInt;
  double ttgammaSFerror = fitError * dataInt / ttgammaInt;

  double ttjetsSF = (1. - fitVal) * dataInt / ttjetsInt;
  double ttjetsSFerror = fitError * dataInt / ttjetsInt;

  cout << "-------------------------------------------------------------" << endl;
  cout << "ttjetsSF = " << ttjetsSF << " +/- " << ttjetsSFerror << endl;
  cout << "ttgammaSF = " << ttgammaSF << " +/- " << ttgammaSFerror << endl;
  cout << "-------------------------------------------------------------" << endl << endl;
  

}
				    

