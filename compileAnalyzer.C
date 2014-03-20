void compileAnalyzer() {

  gROOT->Reset();
  gSystem->Load("libSusyEvent.so");
  gROOT->LoadMacro("SusyEventAnalyzer.cc++");

}
