void makePlots() {

  gROOT->LoadMacro("analyze.C+");

  TStopwatch ts;
  ts.Start();

  TString input = "FILE_TO_RUN";
  bool addMC = true;
  //int intLumi = 19712; // quote to 19.7
  int intLumi = 876.225;

  double metCut = -1.;

  bool displayKStest = true;
  bool blinded = false;
  int nPhotons_req = 0;

  for(int i = 0; i < 8; i++) {
    analyze(input, addMC, i, intLumi, metCut, nPhotons_req, displayKStest, blinded);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
