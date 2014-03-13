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

  std::vector<TString> ele_trigger;
  ele_trigger.push_back("HLT_Ele27_WP80_v");
  std::vector<int> ele_type;
  ele_type.push_back(1);
  sea->AddHlt(ele_trigger, ele_type);

  std::vector<TString> mu_trigger;
  mu_trigger.push_back("HLT_IsoMu24_eta2p1_v");
  std::vector<int> mu_type;
  mu_type.push_back(2);
  sea->AddHlt(mu_trigger, mu_type);

  std::vector<TString> qcd_mu_trigger;
  qcd_mu_trigger.push_back("HLT_Mu24_eta2p1_v");
  std::vector<int> qcd_mu_type;
  qcd_mu_type.push_back(3);
  sea->AddHlt(qcd_mu_trigger, qcd_mu_type);

  sea->SetUseTrigger(true);

  sea->SetProcessNEvents(-1);      	  // number of events to be processed
  
  sea->SetIsMC(isMC);
  sea->SetIsFastSim(isFastSim);
  sea->SetScanName(scan);
  sea->SetDoPileupReweighting(true);
  sea->SetUseJson(false);
  sea->SetDoBtagScaling(false);
  sea->SetBtagTechnicalStop("ABCD");  

  sea->SetRejectFakeElectrons(false);

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
