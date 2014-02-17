void makePlots() {

  gROOT->LoadMacro("analyze_mc.C+");

  TStopwatch ts;
  ts.Start();

  TString input = "FILE_TO_RUN";
  bool addMC = true;
  int intLumi = 19712; // quote to 19.7

  double metCut = -1.;

  bool useTTbar = false;
  bool useTTMBD = false;
  bool displayKStest = true;
  bool blinded = true;

  for(int i = 0; i < 8; i++) {
    analyze(input, addMC, i, intLumi, metCut, useTTbar, useTTMBD, displayKStest, blinded);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
