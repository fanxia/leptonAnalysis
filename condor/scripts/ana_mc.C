void ana_mc(TString scan = "DATASETNAME", TString discriminant = "CSVM", bool isMC = true, bool isFastSim = false) {

  gROOT->Reset();
  gSystem->Load("libSusyEvent.so");

  gROOT->LoadMacro("SusyEventAnalyzer.cc++");

  char* tmp = getenv("CONDOR_SECTION");
  int index = atoi ( tmp );
  cout << "index = " << index << endl;

  TChain chain("susyTree");

  

  chain.SetBranchStatus("*", 1);

  if(chain.LoadTree(0) != 0) {
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);
  sea->SetUseDPhiCut(false);

  // configuration parameters
  // any values given here will replace the default values
  sea->SetPrintInterval(1e4);             // print frequency
  sea->SetPrintLevel(0);                  // print level for event contents

  std::vector<TString> eg_names;
  eg_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  std::vector<int> eg_types;
  eg_types.push_back(1);
  eg_types.push_back(2);
  eg_types.push_back(3);
  sea->AddHlt(eg_names, eg_types);

  std::vector<TString> f_names;
  f_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  f_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  f_names.push_back("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  f_names.push_back("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
  std::vector<int> f_types;
  f_types.push_back(4);
  f_types.push_back(5);
  f_types.push_back(6);
  sea->AddHlt(f_names, f_types);

  sea->SetUseTrigger(true);

  sea->SetProcessNEvents(-1);      	  // number of events to be processed
  
  sea->SetIsMC(isMC);
  sea->SetIsFastSim(isFastSim);
  sea->SetScanName(scan);
  sea->SetDoPileupReweighting(true);
  sea->SetUseJson(false);
  sea->SetDoBtagScaling(false);
  sea->SetBtagTechnicalStop("ABCD");  

  sea->SetRejectFakeElectrons(true);

  sea->SetUseSyncFile(false);
  //sea->IncludeSyncFile("synchro/dmorse_ff.txt_not_brian_ff_nojet.txt");
  sea->SetCheckSingleEvent(false);
  sea->AddCheckSingleEvent(196203, 33, 27883630);

  sea->SetBtagger(discriminant);

  sea->AddValidTagger("TCHPT");
  sea->AddValidTagger("JPL");
  sea->AddValidTagger("JPM");
  sea->AddValidTagger("JPT");
  sea->AddValidTagger("CSVL");
  sea->AddValidTagger("CSVM");
  sea->AddValidTagger("CSVT");
  
  TStopwatch ts;

  ts.Start();

  TString stage = "STAGING";

  if(stage == "pileup") {
    std::cout << std::endl << "PileupWeights()" << std::endl;
    sea->PileupWeights("jan3_pileup.root");
  }
  else if(stage == "btag") {
    std::cout << std::endl << "CalculateBtagEfficiency()" << std::endl;
    sea->CalculateBtagEfficiency();
  }
  else if(stage == "acceptance") {
    std::cout << std::endl << "Acceptance()" << std::endl;
    sea->Acceptance();
  }
  else if(stage == "all") {
    std::cout << std::endl << "PileupWeights()" << std::endl;
    sea->PileupWeights("jan3_pileup.root");
    std::cout << std::endl << "CalculateBtagEfficiency()" << std::endl;
    sea->CalculateBtagEfficiency();
    std::cout << std::endl << "Acceptance()" << std::endl;
    sea->Acceptance();
  }
  
  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;
  
}
