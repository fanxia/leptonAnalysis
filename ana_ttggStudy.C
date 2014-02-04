void ana_ttggStudy(TString scan = "stop-bino", TString discriminant = "CSVM", bool isMC = true, bool isFastSim = true) {

  gROOT->Reset();
  gSystem->Load("libSusyEvent.so");
  
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");

  TChain chain("susyTree");

  char* tmp = getenv("CONDOR_SECTION");
  cout << "tmp = " << tmp << endl;
  int index = atoi ( tmp );
  cout << "index = " << index << endl;

  //int index1 = int(index)/17*100 + 400;
  //int index2 = int(index)%17*100 + 420;

  Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mGluino = 5050;
  Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  Double_t mst_ext[22] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1500, 2000, 5000};
  Double_t mGluino_ext = 5025;
  Double_t mBino_ext[5] = {25, 50, 75, 100, 125};

  // stop-bino scan
  int index1 = 0;
  int index2 = 0;

  if(scan == "stop-bino") {
    index1 = mst[int(index)/31];
    index2 = mBino[int(index)%31];
  }
  else if(scan == "stop-bino extended") {
    index1 = mst_ext[int(index)/5];
    index2 = mBino_ext[int(index)%5];
  }

  cout << "index1 = " << index1 << endl;
  cout << "index2 = " << index2 << endl;

  char input_file[500];

  if(scan == "stop-bino") sprintf(input_file, "dcap:///pnfs/cms/WAX/resilient/bfrancis/SusyNtuples/cms538v1/naturalBinoNLSP/tree_naturalBinoNLSPout_mst_%d_M3_5050_M1_%d.root", index1, index2);
  else if(scan == "stop-bino extended") sprintf(input_file, "/eos/uscms/store/user/lpcpjm/PrivateMC/FastSim/533p3_full/naturalBinoNLSP_try3/SusyNtuple/cms533v1_v1/tree_naturalBinoNLSPout_mst_%d_M3_5025_M1_%d.root", index1, index2);

  //sprintf(input_file, "/eos/uscms/store/user/lpcpjm/PrivateMC/FastSim/525p1v3/Spectra_gsq_B/SusyNtuple/cms533v1_v1/tree_%d_%d_375.root", index1, index2);
  //sprintf(input_file, "/eos/uscms/store/user/lpcpjm/PrivateMC/FastSim/525p1/Spectra_gsq_W/SusyNtuple/cms525v2_v1/tree_%d_%d_375.root", index1, index2);

  cout << "input_file = " << input_file << endl;

  chain.Add(input_file);

  chain.SetBranchStatus("*", 1);

  if(chain.LoadTree(0) != 0) {
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  SusyEventAnalyzer * sea = new SusyEventAnalyzer(chain);
  sea->SetUseDPhiCut(false);

  sea->SetScanName(scan);
  sea->SetPrintInterval(1e4);
  sea->SetPrintLevel(0);
  sea->SetUseTrigger(true);

  std::vector<TString> eg_names;
  eg_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50");
  std::vector<int> eg_types;
  eg_types.push_back(1);
  eg_types.push_back(2);
  eg_types.push_back(3);
  sea->AddHlt(eg_names, eg_types);

  std::vector<TString> f_names;
  f_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50");
  f_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85");
  f_names.push_back("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50");
  f_names.push_back("HLT_Photon36_R9Id85_Photon22_R9Id85");
  std::vector<int> f_types;
  f_types.push_back(4);
  f_types.push_back(5);
  f_types.push_back(6);
  sea->AddHlt(f_names, f_types);

  sea->SetProcessNEvents(-1);

  sea->SetUseJson(false);
  sea->SetIsMC(isMC);
  sea->SetIsFastSim(isFastSim);
  sea->SetDoBtagScaling(false);
  sea->SetBtagTechnicalStop("ABCD");

  sea->SetRejectFakeElectrons(false);

  sea->SetUseSyncFile(false);
  //sea->IncludeSyncFile("utility/syncList.txt");
  sea->SetCheckSingleEvent(false);
  sea->AddCheckSingleEvent(166890, 399, 436836287);

  sea->SetBtagger(discriminant);
  sea->SetDoPileupReweighting(false);

  TStopwatch ts;

  ts.Start();

  sea->ttggStudy();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}

