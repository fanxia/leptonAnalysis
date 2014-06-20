void makePlots() {

  gSystem->AddIncludePath("-I$ROOFITSYS/include"); 
  gROOT->LoadMacro("analyze.C+");

  TStopwatch ts;
  ts.Start();

  TString input_ele = "ELE_FILE_TO_RUN";
  TString input_muon = "MUON_FILE_TO_RUN";
  bool addMC = true;
  int intLumi = 19712; // quote to 19.7

  double metCut = -1.;

  bool displayKStest = true;
  bool blinded = true;
  int nPhotons_req = NUM_PHOTONS_REQUIRED;

  const int nChannels = 4;
  TString channels[nChannels] = {"ele_jjj", "ele_bjj",
				 "muon_jjj", "muon_bjj"};
  int nBtagReq[nChannels] = {0, 1,
			     0, 1};

  for(int i = 0; i < nChannels; i++) {
    if(i != 1 && i != 3) continue;
    if(i < 2) {
      fitPhotons(input_ele, addMC, i, intLumi, metCut, nPhotons_req, nBtagReq[i], displayKStest, blinded, 0);
      analyze(input_ele, addMC, i, intLumi, metCut, nPhotons_req, nBtagReq[i], displayKStest, blinded, 0);
    }
    else {
      fitPhotons(input_muon, addMC, i, intLumi, metCut, nPhotons_req, nBtagReq[i], displayKStest, blinded, 0);
      analyze(input_muon, addMC, i, intLumi, metCut, nPhotons_req, nBtagReq[i], displayKStest, blinded, 0);
    }
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
