void ana_genTreeViewerRA3() {

  gROOT->Reset();
  gSystem->Load("libSusyEvent.so");

  gROOT->LoadMacro("GenTreeViewerRA3.cc+");

  TChain chain("susyTree");
  chain.Add("/eos/uscms/store/user/bfrancis/SusyNtuples/cms538v1/LHE2EDM_WHIZARD_2to5_ttA/susyEvents_48_1_t0v.root");
  
  chain.SetBranchStatus("*", 1);

  if(chain.LoadTree(0) != 0) {
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  GenTreeViewerRA3 * gtv = new GenTreeViewerRA3(chain);

  gtv->SetProcessNEvents(10);
  gtv->SetShowMass(false);
  gtv->SetMinPt(2.);

  gtv->viewGenTreeRA3();

}
