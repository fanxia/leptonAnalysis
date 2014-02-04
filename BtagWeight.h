#ifndef BtagWeight_h
#define BtagWeight_h

#include <math.h>
#include <iostream>
#include <vector>

#include "BtagInfo.h"

// Adapted from https://twiki.cern.ch/twiki/pub/CMS/BTagWeight/BTagWeight2.cpp

using namespace std;

class BtagWeight {
 public:
 BtagWeight(int jmin) :
  minTags(jmin) {}

  bool filter(int t) { return (t >= minTags); };
  pair<float, float> weight(vector<BtagInfo> info, int tags, double scale, bool doScale);

 private:
  int minTags;
};

pair<float, float> BtagWeight::weight(vector<BtagInfo> info, int tags, double scale, bool doScale) {
  
  if(!filter(tags)) return make_pair(0., 0.);

  int njets = info.size();
  int comb = 1 << njets;
  float pMC = 0.;
  float pMC_err2 = 0.;
  float pData = 0.;
  float pData_err2 = 0.;

  for(int i = 0; i < comb; i++) {
    float mc = 1.;
    float mc_err2 = 0.;
    float data = 1.;
    float data_err2 = 0.;
    int ntagged = 0;
    
    for(int j = 0; j < njets; j++) {
      bool tagged = ((i >> j) & 0x1) == 1;

      float old_mc = mc;
      float old_data = data;
      float old_mc_err2 = mc_err2;
      float old_data_err2 = data_err2;

      if(tagged) {
	ntagged++;
	mc *= info[j].GetTaggingEfficiency();
	data *= info[j].GetTaggingEfficiency() * info[j].GetScaleFactor(scale, doScale);
      }
      else {
	mc *= (1. - info[j].GetTaggingEfficiency());
	data *= (1. - info[j].GetTaggingEfficiency() * info[j].GetScaleFactor(scale, doScale));
      }

      mc_err2 = old_mc*old_mc*info[j].GetTaggingEfficiencyError()*info[j].GetTaggingEfficiencyError();
      mc_err2 += info[j].GetTaggingEfficiency()*info[j].GetTaggingEfficiency()*old_mc_err2*old_mc_err2;
      
      data_err2 = old_data*old_data*info[j].GetTaggingEfficiency()*info[j].GetScaleFactor(scale, doScale)*info[j].GetTaggingEfficiency()*info[j].GetScaleFactor(scale, doScale);
      data_err2 += info[j].GetTaggingEfficiency()*info[j].GetScaleFactor(scale, doScale)*old_data_err2*old_data_err2;

    }

    if(filter(ntagged)) {
      pMC += mc;
      pMC_err2 += mc_err2;
      pData += data;
      pData_err2 += data_err2;
    }

  }

  float totalerr2 = pData*pData*pData_err2 + pMC*pMC*pMC_err2;

  if(pMC == 0) return make_pair(0, 0);
  return make_pair(pData/pMC, sqrt(totalerr2));
}

#endif
