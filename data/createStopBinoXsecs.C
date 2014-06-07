void createStopBinoXsecs() {

  TFile * fOut = new TFile("stop-bino_xsecs.root", "RECREATE");

  Double_t mst_bins[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mBino_bins[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst_bins[i] + mst_bins[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino_bins[i] + mBino_bins[i-1])/2.;
  ybins[32] = 2175;

  TH2D * h_real_xsec = new TH2D("real_xsec", "real_xsec", 30, xbins, 32, ybins);
  h_real_xsec->Sumw2();

  TH2D * h_real_errors = (TH2D*)h_real_xsec->Clone("real_errors");

  TH1D * h_xsec = new TH1D("xsec", "xsec", 380, 100, 2000);
  h_xsec->Sumw2();
  
  TH1D * h_errors = (TH1D*)h_xsec->Clone("errors");

  ifstream fin;
  fin.open("stop-bino_xsecs.dat");

  while(1) {
    int mStop;
    double xsec, xsecErr;

    fin >> mStop >> xsec >> xsecErr;
    
    if(!fin.good()) break;

    h_xsec->SetBinContent(h_xsec->FindBin(mStop), xsec);
    h_errors->SetBinContent(h_errors->FindBin(mStop), xsecErr);
  }
  
  fin.close();

  fin.open("realStopMass.dat");

  while(1) {
    int mStop, mBino;
    double realMass;

    fin >> mStop >> mBino >> realMass;

    if(!fin.good()) break;

    if(mBino > mStop) continue;

    int bin = h_xsec->FindBin(realMass);

    double xsec_left = h_xsec->GetBinContent(bin);
    double error_left = h_xsec->GetBinError(bin);

    double xsec_right = h_xsec->GetBinContent(bin + 1);
    double error_right = h_xsec->GetBinContent(bin + 1);

    double left = h_xsec->GetBinLowEdge(bin);
    double right = h_xsec->GetBinLowEdge(bin + 1);

    double alpha = (realMass - left) / (right - left);
    double beta = (right - realMass) / (right - left);

    double xsec = alpha * xsec_right + beta * xsec_left;
    double error = sqrt(alpha*alpha*error_right*error_right + beta*beta*error_left*error_left);

    h_real_xsec->SetBinContent(h_real_xsec->FindBin(mStop, mBino), xsec);
    h_real_errors->SetBinContent(h_real_errors->FindBin(mStop, mBino), error);
  }

  fin.close();

  fOut->Write();
  fOut->Close();

}
