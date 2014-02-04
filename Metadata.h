#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <vector>

#include "../src/SusyEvent.h"

using namespace std;

enum eventTypes {cNothing, cGG, cEG, cEE, cFF, cGF, cEF};

const int nCategories = 6;
TString categories[nCategories] = {"gg", "eg", "ff", "gf", "ee", "ef"};
const int nChannels = 8;

TString channels[nChannels] = {"nojet", "j", "b", "jj", "bj", "muJets", "eleJets", "hadronic"};
unsigned int nJetReq[nChannels] = {0, 1, 1, 2, 2, 2, 2, 4};
unsigned int nBtagReq[nChannels] = {0, 0, 1, 0, 1, 1, 1, 1};
int nEleReq[nChannels] = {-1, -1, -1, -1, -1, 0, 1, 0};
int nMuonReq[nChannels] = {-1, -1, -1, -1, -1, 1, 0, 0};

typedef std::vector<std::vector<TH1F*> > VTH1F;
typedef std::vector<std::vector<TH2F*> > VTH2F;

Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
Double_t mGluino = 5050;
Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

Double_t mst_ext[22] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1500, 2000, 5000};
Double_t mGluino_ext = 5025;
Double_t mBino_ext[5] = {25, 50, 75, 100, 125};

VTH1F BookTH1FVector(TString name, TString title, Int_t nBinsX, Double_t xlo, Double_t xhi, int nCategs, TString categs[], int nChans, TString chans[]) {
  
  VTH1F vv;
  for(int i = 0; i < nCategs; i++) {
    vector<TH1F*> v;
    for(int j = 0; j < nChans; j++) {
      v.push_back(new TH1F(name+"_"+categs[i]+"_"+chans[j], title, nBinsX, xlo, xhi));
      v[j]->Sumw2();
    }
    vv.push_back(v);
    v.clear();
  }
  
  return vv;
}

VTH2F BookTH2FVector(TString name, TString title, Int_t nBinsX, Double_t xlo, Double_t xhi, Int_t nBinsY, Double_t ylo, Double_t yhi, int nCategs, TString categs[], int nChans, TString chans[]) {

  VTH2F vv;
  for(int i = 0; i < nCategs; i++) {
    vector<TH2F*> v;
    for(int j = 0; j < nChans; j++) {
      v.push_back(new TH2F(name+"_"+categs[i]+"_"+chans[j], title, nBinsX, xlo, xhi, nBinsY, ylo, yhi));
    }
    vv.push_back(v);
    v.clear();
  }

  return vv;
}

TString FormatName(TString scan) {

  char * tmp = getenv("CONDOR_SECTION");
  int index = atoi(tmp);

  int index1, index2;

  bool stopScan = (scan == "stop-bino");
  bool stopScan_ext = (scan == "stop-bino extended");
  bool squarkGluinoScan = (scan == "squark-gluino");

  if(stopScan) {
    index1 = mst[int(index)/31];
    index2 = mBino[int(index)%31];
  }
  else if(stopScan_ext) {
    index1 = mst_ext[int(index)/5];
    index2 = mBino_ext[int(index)%5];
  }
  else if(squarkGluinoScan) {
    index1 = int(index)/17*100 + 400;
    index2 = int(index)%17*100 + 420;
  }
  else {
    index1 = 0;
    index2 = 0;
  }

  char output_code[100];
  if(stopScan || stopScan_ext) sprintf(output_code, "_mst_%d_m1_%d", index1, index2);
  else if(squarkGluinoScan) sprintf(output_code, "_mS%d_mG%d_mN375", index1, index2);
  
  TString output_code_t = output_code;
  if(!stopScan && !stopScan_ext && !squarkGluinoScan) output_code_t = "_"+scan;

  return output_code_t;
}

