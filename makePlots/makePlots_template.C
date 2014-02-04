void makePlots() {

  gROOT->LoadMacro("analyze.C+");

  TStopwatch ts;
  ts.Start();

  TString input = "FILE_TO_RUN";
  bool addMC = true;
  TString intLumi = "19.7";
  int intLumi_int = 19.712;
  
  bool useFF = true;
  bool useDifferenceSystematic = false;

  double metCut = -1.;

  bool displayKStest = true;

  for(int i = 0; i < 8; i++) {
    mvaTreeMaker(input, i);
    analyze(input, addMC, i, intLumi, intLumi_int, useFF, useDifferenceSystematic, metCut, displayKStest);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
