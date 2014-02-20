void ana_condor(TString discriminant = "CSVM", bool isMC = false) {

  gROOT->Reset();
  gSystem->Load("libSusyEvent.so");

  gROOT->LoadMacro("SusyEventAnalyzer.cc++");

  TChain chain("susyTree");
  chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/singleElectronA_v1/susyEvents_100_1_9iw.root");
  
  chain.SetBranchStatus("*", 1);

  if(chain.LoadTree(0) != 0) {
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);

  sea->SetUseDPhiCut(false);

  // configuration parameters
  // any values given here will replace the default values
  sea->SetPrintInterval(1e5);             // print frequency
  sea->SetPrintLevel(0);                  // print level for event contents

  std::vector<TString> ele_trigger;
  ele_trigger.push_back("HLT_Ele27_WP80_v");
  //ele_trigger.push_back("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  std::vector<int> ele_type;
  ele_type.push_back(1);
  sea->AddHlt(ele_trigger, ele_type);

  std::vector<TString> mu_trigger;
  mu_trigger.push_back("HLT_IsoMu24_eta2p1_v");
  std::vector<int> mu_type;
  mu_type.push_back(2);
  sea->AddHlt(mu_trigger, mu_type);

  sea->SetUseTrigger(true);

  sea->SetProcessNEvents(-1);      	  // number of events to be processed
  
  sea->IncludeAJson("Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");

  sea->SetIsMC(isMC);
  
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
  
  sea->Data();
  
  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;
  
}
