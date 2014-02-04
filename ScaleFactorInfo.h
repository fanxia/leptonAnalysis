#ifndef ScaleFactorInfo_h
#define ScaleFactorInfo_h

#include <TString.h>
#include <TF1.h>
#include <TH1.h>

#include <iostream>
#include <vector>

using namespace std;

// Basic container for all scale factors, errors, and efficiencies
class ScaleFactorInfo {
 public:
  ScaleFactorInfo(TString tag);
  virtual ~ScaleFactorInfo() {};

  Float_t GetSFbottom(Float_t pt);
  Float_t GetSFcharm(Float_t pt) { return GetSFbottom(pt); }; // SFc = SFb with twice the quoted uncertainty

  // SFlightMin and SFlightMax correspond to SFlight +- (stat+syst error)
  Float_t GetSFlight(TString meanminmax, Float_t eta, TString DataPeriod, Float_t pt);
  Float_t GetSFlightMin(Float_t eta, TString DataPeriod, Float_t pt)  { return GetSFlight("min", eta, DataPeriod, pt); };
  Float_t GetSFlightMean(Float_t eta, TString DataPeriod, Float_t pt) { return GetSFlight("mean", eta, DataPeriod, pt); };
  Float_t GetSFlightMax(Float_t eta, TString DataPeriod, Float_t pt)  { return GetSFlight("max", eta, DataPeriod, pt); };
  
  // Tagging efficiencies. MC reweighting requires knowledge of all tagging efficiencies for that sample beforehand
  void SetTaggingEfficiencies(TH1F * eff_light, TH1F * eff_charm, TH1F * eff_bottom) {
    hasEfficiencies = true;
    lEff = (TH1F*)eff_light->Clone();
    cEff = (TH1F*)eff_charm->Clone();
    bEff = (TH1F*)eff_bottom->Clone();
  }

  Float_t GetEffLight(Float_t pt) { return (hasEfficiencies) ? lEff->GetBinContent(lEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffCharm(Float_t pt) { return (hasEfficiencies) ? cEff->GetBinContent(cEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffBottom(Float_t pt) { return (hasEfficiencies) ? bEff->GetBinContent(bEff->GetXaxis()->FindBin(pt)) : 0.; };

  Float_t GetEffErrorLight(Float_t pt) { return (hasEfficiencies) ? lEff->GetBinError(lEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffErrorCharm(Float_t pt) { return (hasEfficiencies) ? cEff->GetBinError(cEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffErrorBottom(Float_t pt) { return (hasEfficiencies) ? bEff->GetBinError(bEff->GetXaxis()->FindBin(pt)) : 0.; };

  // FastSim to FullSim correction factors
  vector<Float_t> GetCFbottom() { return CFb; };
  vector<Float_t> GetCFcharm() { return CFc; };
  vector<Float_t> GetCFlight(Float_t eta) {
    if(eta < 1.2) return CFudsg_12;
    return CFudsg_24;
  };

  vector<Float_t> GetCFbottomErrors() { return CFb_errors; };
  vector<Float_t> GetCFcharmErrors() { return CFc_errors; };
  vector<Float_t> GetCFlightErrors(Float_t eta) {
    if(eta < 1.2) return CFudsg_errors_12;
    return CFudsg_errors_24;
  };

  vector<Float_t> GetSFbottomErrors() { return SFb_error; };
  vector<Float_t> GetSFcharmErrors() { // SFc = SFb with twice the quoted uncertainty
    vector<Float_t> err = GetSFbottomErrors();
    for(unsigned int i = 0; i < err.size(); i++) err[i] *= 2.;
    return err;
  }
 
 private:
  TString tagger;

  bool hasEfficiencies;
  TH1F * lEff;
  TH1F * cEff;
  TH1F * bEff;

  vector<Float_t> CFb;
  vector<Float_t> CFc;
  vector<Float_t> CFudsg_12;
  vector<Float_t> CFudsg_24;

  vector<Float_t> CFb_errors;
  vector<Float_t> CFc_errors;
  vector<Float_t> CFudsg_errors_12;
  vector<Float_t> CFudsg_errors_24;

  vector<Float_t> SFb_error;

};

ScaleFactorInfo::ScaleFactorInfo(TString tag) {

  tagger = tag;
  hasEfficiencies = false;

  // DURP: FAST SIM CF COMING SOON FOR EPS13, STILL MORIOND13

  if(tagger == "CSVM") {
    CFb.push_back(0.982194); 
    CFb.push_back(0.980998); 
    CFb.push_back(0.992014); 
    CFb.push_back(0.994472); 
    CFb.push_back(0.996825); 
    CFb.push_back(0.999822); 
    CFb.push_back(1.00105); 
    CFb.push_back(1.00023); 
    CFb.push_back(0.991994); 
    CFb.push_back(0.979123); 
    CFb.push_back(0.947207); 
    CFb.push_back(0.928006); 
    CFb.push_back(0.874260); 
    CFb.push_back(0.839610);
  }
  else if(tagger == "CSVT") {
    CFb.push_back(0.968117); 
    CFb.push_back(0.967822); 
    CFb.push_back(0.978278); 
    CFb.push_back(0.981281); 
    CFb.push_back(0.987679); 
    CFb.push_back(0.986590); 
    CFb.push_back(0.990246); 
    CFb.push_back(0.984504); 
    CFb.push_back(0.967024); 
    CFb.push_back(0.940042); 
    CFb.push_back(0.873019); 
    CFb.push_back(0.850847); 
    CFb.push_back(0.769561); 
    CFb.push_back(0.650192);
  }
  else CFb.clear();

  if(tagger == "CSVM") {
    CFc.push_back(0.988545); 
    CFc.push_back(0.981714); 
    CFc.push_back(1.00946); 
    CFc.push_back(1.01591); 
    CFc.push_back(1.02810); 
    CFc.push_back(1.02195); 
    CFc.push_back(1.02590); 
    CFc.push_back(1.01936); 
    CFc.push_back(0.991228); 
    CFc.push_back(0.955343); 
    CFc.push_back(0.944433); 
    CFc.push_back(0.917282); 
    CFc.push_back(0.935018); 
    CFc.push_back(1.06375);
  }
  else if(tagger == "CSVT") {
    CFc.push_back(0.960959); 
    CFc.push_back(0.973876); 
    CFc.push_back(0.984323); 
    CFc.push_back(0.996344); 
    CFc.push_back(1.02418); 
    CFc.push_back(0.985580); 
    CFc.push_back(0.994745); 
    CFc.push_back(0.970489); 
    CFc.push_back(0.914155); 
    CFc.push_back(0.872072); 
    CFc.push_back(0.945289); 
    CFc.push_back(0.783816); 
    CFc.push_back(0.942773); 
    CFc.push_back(0.527354);
  }
  else CFc.clear();

  if(tagger == "CSVM") {
    CFudsg_12.push_back(1.21878); 
    CFudsg_12.push_back(1.28615); 
    CFudsg_12.push_back(1.37535); 
    CFudsg_12.push_back(1.38966); 
    CFudsg_12.push_back(1.40320); 
    CFudsg_12.push_back(1.49835); 
    CFudsg_12.push_back(1.44308); 
    CFudsg_12.push_back(1.58198); 
    CFudsg_12.push_back(1.55687); 
    CFudsg_12.push_back(1.65790); 
    CFudsg_12.push_back(1.90233); 
    CFudsg_12.push_back(1.92259); 
    CFudsg_12.push_back(2.66174); 
    CFudsg_12.push_back(3.08688);

    CFudsg_24.push_back(1.46970); 
    CFudsg_24.push_back(1.48732); 
    CFudsg_24.push_back(1.69024); 
    CFudsg_24.push_back(1.64494); 
    CFudsg_24.push_back(1.79297); 
    CFudsg_24.push_back(1.90760); 
    CFudsg_24.push_back(1.99867); 
    CFudsg_24.push_back(2.21659); 
    CFudsg_24.push_back(2.20103); 
    CFudsg_24.push_back(2.42645); 
    CFudsg_24.push_back(2.67594); 
    CFudsg_24.push_back(4.24735); 
    CFudsg_24.push_back(3.98979); 
    CFudsg_24.push_back(15.0457);
  }
  else if(tagger == "CSVT") {
    CFudsg_12.push_back(1.24890); 
    CFudsg_12.push_back(1.35145); 
    CFudsg_12.push_back(1.37205); 
    CFudsg_12.push_back(1.32472); 
    CFudsg_12.push_back(1.39976); 
    CFudsg_12.push_back(1.45884); 
    CFudsg_12.push_back(1.59912); 
    CFudsg_12.push_back(1.58971); 
    CFudsg_12.push_back(1.30841); 
    CFudsg_12.push_back(1.55936); 
    CFudsg_12.push_back(1.28346); 
    CFudsg_12.push_back(2.21265); 
    CFudsg_12.push_back(2.06927); 
    CFudsg_12.push_back(2.88109);

    CFudsg_24.push_back(1.67634); 
    CFudsg_24.push_back(1.70105); 
    CFudsg_24.push_back(1.75999); 
    CFudsg_24.push_back(1.78459); 
    CFudsg_24.push_back(2.19343); 
    CFudsg_24.push_back(2.73199); 
    CFudsg_24.push_back(3.49277); 
    CFudsg_24.push_back(2.58863); 
    CFudsg_24.push_back(2.48824); 
    CFudsg_24.push_back(4.01723); 
    CFudsg_24.push_back(3.86956); 
    CFudsg_24.push_back(0.000456049); 
    CFudsg_24.push_back(2.30988); 
    CFudsg_24.push_back(0.000855693);
  }
  else {
    CFudsg_12.clear();
    CFudsg_24.clear();
  }

  if(tagger == "CSVM") {
    CFb_errors.push_back(0.00253112); 
    CFb_errors.push_back(0.00296453); 
    CFb_errors.push_back(0.00113963); 
    CFb_errors.push_back(0.00128363); 
    CFb_errors.push_back(0.00232566); 
    CFb_errors.push_back(0.00232353); 
    CFb_errors.push_back(0.00219086); 
    CFb_errors.push_back(0.00156856); 
    CFb_errors.push_back(0.00322279); 
    CFb_errors.push_back(0.00400414); 
    CFb_errors.push_back(0.00737465); 
    CFb_errors.push_back(0.0105033); 
    CFb_errors.push_back(0.0171706); 
    CFb_errors.push_back(0.0344172);
  }
  else if(tagger == "CSVT") {
    CFb_errors.push_back(0.00223422); 
    CFb_errors.push_back(0.00367427); 
    CFb_errors.push_back(0.00145554); 
    CFb_errors.push_back(0.00337572); 
    CFb_errors.push_back(0.00344106); 
    CFb_errors.push_back(0.00591257); 
    CFb_errors.push_back(0.00218050); 
    CFb_errors.push_back(0.00472939); 
    CFb_errors.push_back(0.00353119); 
    CFb_errors.push_back(0.00739502); 
    CFb_errors.push_back(0.0193330); 
    CFb_errors.push_back(0.0158257); 
    CFb_errors.push_back(0.0306048); 
    CFb_errors.push_back(0.0603701);
  }
  else CFb_errors.clear();

  if(tagger == "CSVM") {
    CFc_errors.push_back(0.00746259); 
    CFc_errors.push_back(0.00661831); 
    CFc_errors.push_back(0.00968682); 
    CFc_errors.push_back(0.00751322); 
    CFc_errors.push_back(0.00675507); 
    CFc_errors.push_back(0.00562821); 
    CFc_errors.push_back(0.00862890); 
    CFc_errors.push_back(0.00768003); 
    CFc_errors.push_back(0.0188981); 
    CFc_errors.push_back(0.0261163); 
    CFc_errors.push_back(0.0450601); 
    CFc_errors.push_back(0.0448453); 
    CFc_errors.push_back(0.148805); 
    CFc_errors.push_back(0.177157);
  }
  else if(tagger == "CSVT") {
    CFc_errors.push_back(0.0155733); 
    CFc_errors.push_back(0.0121900); 
    CFc_errors.push_back(0.0131678); 
    CFc_errors.push_back(0.0113739); 
    CFc_errors.push_back(0.0213937); 
    CFc_errors.push_back(0.0123294); 
    CFc_errors.push_back(0.0153230); 
    CFc_errors.push_back(0.0156350); 
    CFc_errors.push_back(0.0409568); 
    CFc_errors.push_back(0.0654966); 
    CFc_errors.push_back(0.112785); 
    CFc_errors.push_back(0.187795); 
    CFc_errors.push_back(0.331301); 
    CFc_errors.push_back(0.162462);
  }
  else CFc_errors.clear();

  if(tagger == "CSVM") {
    CFudsg_errors_12.push_back(0.0182686); 
    CFudsg_errors_12.push_back(0.0373732); 
    CFudsg_errors_12.push_back(0.0461870); 
    CFudsg_errors_12.push_back(0.0288973); 
    CFudsg_errors_12.push_back(0.0333528); 
    CFudsg_errors_12.push_back(0.0513836); 
    CFudsg_errors_12.push_back(0.0420353); 
    CFudsg_errors_12.push_back(0.106627); 
    CFudsg_errors_12.push_back(0.0658359); 
    CFudsg_errors_12.push_back(0.117285); 
    CFudsg_errors_12.push_back(0.185533); 
    CFudsg_errors_12.push_back(0.214071); 
    CFudsg_errors_12.push_back(0.487274); 
    CFudsg_errors_12.push_back(0.871502);
    
    CFudsg_errors_24.push_back(0.104716); 
    CFudsg_errors_24.push_back(0.0392025); 
    CFudsg_errors_24.push_back(0.106315); 
    CFudsg_errors_24.push_back(0.115751); 
    CFudsg_errors_24.push_back(0.106807); 
    CFudsg_errors_24.push_back(0.0642086); 
    CFudsg_errors_24.push_back(0.138742); 
    CFudsg_errors_24.push_back(0.182345); 
    CFudsg_errors_24.push_back(0.169922); 
    CFudsg_errors_24.push_back(0.297889); 
    CFudsg_errors_24.push_back(0.320088); 
    CFudsg_errors_24.push_back(0.927736); 
    CFudsg_errors_24.push_back(1.24666); 
    CFudsg_errors_24.push_back(15.1860);
  }
  else if(tagger == "CSVT") {
    CFudsg_errors_12.push_back(0.0751438); 
    CFudsg_errors_12.push_back(0.0651619); 
    CFudsg_errors_12.push_back(0.0604241); 
    CFudsg_errors_12.push_back(0.0726285); 
    CFudsg_errors_12.push_back(0.0968158); 
    CFudsg_errors_12.push_back(0.0931768); 
    CFudsg_errors_12.push_back(0.163039); 
    CFudsg_errors_12.push_back(0.187749); 
    CFudsg_errors_12.push_back(0.198200); 
    CFudsg_errors_12.push_back(0.465354); 
    CFudsg_errors_12.push_back(0.339473); 
    CFudsg_errors_12.push_back(1.07079); 
    CFudsg_errors_12.push_back(1.07723); 
    CFudsg_errors_12.push_back(2.53188);
  
    CFudsg_errors_24.push_back(0.222165); 
    CFudsg_errors_24.push_back(0.161403); 
    CFudsg_errors_24.push_back(0.112342); 
    CFudsg_errors_24.push_back(0.275101); 
    CFudsg_errors_24.push_back(0.364229); 
    CFudsg_errors_24.push_back(0.330588); 
    CFudsg_errors_24.push_back(1.00953); 
    CFudsg_errors_24.push_back(0.404417); 
    CFudsg_errors_24.push_back(1.07731); 
    CFudsg_errors_24.push_back(2.65686); 
    CFudsg_errors_24.push_back(3.18286); 
    CFudsg_errors_24.push_back(5.25051e-05); 
    CFudsg_errors_24.push_back(2.38652); 
    CFudsg_errors_24.push_back(0.000438728);
  }
  else {
    CFudsg_errors_12.clear();
    CFudsg_errors_24.clear();
  }

  if(tagger == "CSVL") {
    Float_t SFb_error_EPS13[16] = {0.033299, 0.0146768, 0.013803, 0.0170145, 0.0166976, 0.0137879, 0.0149072, 0.0153068, 0.0133077, 0.0123737, 0.0157152, 0.0175161, 0.0209241, 0.0278605, 0.0346928, 0.0350099};
    for(int i = 0; i < 16; i++) SFb_error.push_back(SFb_error_EPS13[i]);
  }
  else if(tagger == "CSVM") {
    Float_t SFb_error_EPS13[16] = {0.0415707, 0.0204209, 0.0223227, 0.0206655, 0.0199325, 0.0174121, 0.0202332, 0.0182446, 0.0159777, 0.0218531, 0.0204688, 0.0265191, 0.0313175, 0.0415417, 0.0740446, 0.0596716};
    for(int i = 0; i < 16; i++) SFb_error.push_back(SFb_error_EPS13[i]);
  }
  else if(tagger == "CSVT") {
    Float_t SFb_error_EPS13[16] = {0.0515703, 0.0264008, 0.0272757, 0.0275565, 0.0248745, 0.0218456, 0.0253845, 0.0239588, 0.0271791, 0.0273912, 0.0379822, 0.0411624, 0.0786307, 0.0866832, 0.0942053, 0.102403};
    for(int i = 0; i < 16; i++) SFb_error.push_back(SFb_error_EPS13[i]);
  }
  else SFb_error.clear();

}

Float_t ScaleFactorInfo::GetSFbottom(Float_t pt) {

  Float_t x = pt;

  // EPS13 prescription with ttbar measurements
  if(tagger == "CSVL") return 0.997942*((1.+(0.00923753*x))/(1.+(0.0096119*x)));
  else if(tagger == "CSVM") return (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
  else if(tagger == "CSVT") return (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));

  return 0.;
}

Float_t ScaleFactorInfo::GetSFlight(TString meanminmax, Float_t eta, TString DataPeriod, Float_t pt) {
  
  Float_t eta_ = fabs(eta);
  Float_t x = pt;

  /*
  Double_t ptmax;
  if((tagger == "CSVL" && eta_ < 1.5) ||
     ((tagger == "CSVM" || tagger == "CSVT") && eta_ < 1.6)
     ) ptmax = 800.;
  else ptmax = 700.;
  */

  if(tagger == "CSVL") {
    if(eta_ < 0.5) {
      if( meanminmax == "mean" ) return ((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)));
    }
    else if(eta_ >= 0.5 && eta_ < 1.0) {
      if( meanminmax == "mean" ) return ((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)));
    }
    else if(eta_ >= 1.0 && eta_ < 1.5) {
      if( meanminmax == "mean" ) return ((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)));
    }
    else if(eta_ >= 1.5 && eta_ < 2.4) {
      if( meanminmax == "mean" ) return ((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)));
    }
  } // if CSVL
    
  else if(tagger == "CSVM") {
    if(eta_ < 0.8) {
      if( meanminmax == "mean" ) return ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
    }
    else if(eta_ >= 0.8 && eta_ < 1.6) {
      if( meanminmax == "mean" ) return ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
    }
    else if(eta_ >= 1.6 && eta_ < 2.4) {
      if( meanminmax == "mean" ) return ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
    }
  } // if CSVM

  else if(tagger == "CSVT") {
    if(eta_ < 2.4) {
      if( meanminmax == "mean" ) return ((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
      if( meanminmax == "min" ) return ((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
      if( meanminmax == "max" ) return ((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)));
    }
  } // if CSVT

  return 0.;

}

#endif
