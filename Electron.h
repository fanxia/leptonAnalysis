#include "math.h"
#include "TMath.h"

#include "../src/SusyEvent.h"

bool isLooseElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho, double d0, double dz) {

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());

  float ea;
  if(ele_eta < 1.0) ea = 0.13;
  else if(ele_eta < 1.479) ea = 0.14;
  else if(ele_eta < 2.0) ea = 0.07;
  else if(ele_eta < 2.2) ea = 0.09;
  else if(ele_eta < 2.3) ea = 0.11;
  else if(ele_eta < 2.4) ea = 0.11;
  else ea = 0.14;

  float ele_iso = max(0., (ele.photonIso + ele.neutralHadronIso - rho*ea));
  ele_iso += ele.chargedHadronIso;

  bool inCrack = ele_eta > 1.4442 && ele_eta < 1.566;

  bool passes = ele_eta < 2.5 && 
    !inCrack &&
    ele.momentum.Pt() > 10. &&
    ele.mvaTrig > 0.5 &&
    fabs(d0) < 0.04 &&
    ele_iso / ele.momentum.Pt() < 0.2 &&
    ele.passConversionVeto &&
    ele.nMissingHits <= 0;

  return passes;

}

bool isIsolatedElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho) {

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());

  float ea;
  if(ele_eta < 1.0) ea = 0.13;
  else if(ele_eta < 1.479) ea = 0.14;
  else if(ele_eta < 2.0) ea = 0.07;
  else if(ele_eta < 2.2) ea = 0.09;
  else if(ele_eta < 2.3) ea = 0.11;
  else if(ele_eta < 2.4) ea = 0.11;
  else ea = 0.14;

  float ele_iso = max(0., (ele.photonIso + ele.neutralHadronIso - rho*ea));
  ele_iso += ele.chargedHadronIso;
  
  return (ele_iso / ele.momentum.Pt() < 0.1);

}

bool isAntiIsolatedElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho) {

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());

  float ea;
  if(ele_eta < 1.0) ea = 0.13;
  else if(ele_eta < 1.479) ea = 0.14;
  else if(ele_eta < 2.0) ea = 0.07;
  else if(ele_eta < 2.2) ea = 0.09;
  else if(ele_eta < 2.3) ea = 0.11;
  else if(ele_eta < 2.4) ea = 0.11;
  else ea = 0.14;

  float ele_iso = max(0., (ele.photonIso + ele.neutralHadronIso - rho*ea));
  ele_iso += ele.chargedHadronIso;
  
  float relIso = ele_iso / ele.momentum.Pt();

  return (relIso >= 0.2 && relIso < 1.);

}

bool isTightElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho, double d0, double dz) {

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());

  float ea;
  if(ele_eta < 1.0) ea = 0.13;
  else if(ele_eta < 1.479) ea = 0.14;
  else if(ele_eta < 2.0) ea = 0.07;
  else if(ele_eta < 2.2) ea = 0.09;
  else if(ele_eta < 2.3) ea = 0.11;
  else if(ele_eta < 2.4) ea = 0.11;
  else ea = 0.14;

  float ele_iso = max(0., (ele.photonIso + ele.neutralHadronIso - rho*ea));
  ele_iso += ele.chargedHadronIso;

  bool inCrack = ele_eta > 1.4442 && ele_eta < 1.566;

  bool passes = ele_eta < 2.5 &&
    !inCrack &&
    ele.momentum.Pt() > 30. &&
    fabs(d0) < 0.02 &&
    fabs(dz) < 1.0 &&
    ele.passConversionVeto &&
    ele.nMissingHits <= 0;

  return passes;
}
