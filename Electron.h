#include "math.h"
#include "TMath.h"

#include "../src/SusyEvent.h"

bool isTightElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho, double d0, double dz) {

  if(ele.momentum.Pt() <= 30.) return false; // STUDY THIS (ttH(bb) for now)

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());
  float ea;
  if(ele_eta < 1.0) ea = 0.13;        // ± 0.001
  else if(ele_eta < 1.479) ea = 0.14; // ± 0.002
  else if(ele_eta < 2.0) ea = 0.07;   // ± 0.001
  else if(ele_eta < 2.2) ea = 0.09;   // ± 0.001
  else if(ele_eta < 2.3) ea = 0.11;   // ± 0.002
  else if(ele_eta < 2.4) ea = 0.11;   // ± 0.003
  else ea = 0.14;                     // ± 0.004

  float ele_iso = max(0., (ele.photonIso + ele.neutralHadronIso - rho*ea));
  ele_iso += ele.chargedHadronIso;

  bool tight_barrel = fabs(ele.deltaEtaSuperClusterTrackAtVtx) < 0.004 &&
    fabs(ele.deltaPhiSuperClusterTrackAtVtx) < 0.03 &&
    fabs(ele.sigmaIetaIeta) < 0.01 &&
    ele.hcalOverEcalBc < 0.12 &&
    fabs(d0) < 0.02 &&
    fabs(dz) < 0.1 &&
    fabs(1/(ele.ecalEnergy) - 1/(ele.ecalEnergy/ele.eSuperClusterOverP)) < 0.05 &&
    ele_iso / ele.momentum.Pt() < 0.10 &&
    ele.passConversionVeto &&
    ele.nMissingHits == 0;
     
  bool tight_endcap = fabs(ele.deltaEtaSuperClusterTrackAtVtx) < 0.005 &&
    ele.deltaPhiSuperClusterTrackAtVtx < 0.02 &&
    ele.sigmaIetaIeta < 0.03 &&
    ele.hcalOverEcalBc < 0.1 &&
    d0 < 0.02 &&
    dz < 0.1 &&
    fabs(1/(ele.ecalEnergy) - 1/(ele.ecalEnergy/ele.eSuperClusterOverP)) < 0.05 &&
    ((ele.momentum.Pt() > 20. && ele_iso / ele.momentum.Pt() < 0.10) || (ele.momentum.Pt() <= 20. && ele_iso / ele.momentum.Pt() < 0.07)) &&
    ele.passConversionVeto &&
    ele.nMissingHits == 0;
  
  // mvaNonTrigV0 > 0.9 if we don't use simple cuts-based

  if(ele_eta <= 1.479) return tight_barrel;
  else if(ele_eta > 1.479 && ele_eta < 2.5) return tight_endcap;
  
  return false;

}

bool isLooseElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho, double d0, double dz) {

  if(ele.momentum.Pt() <= 20.) return false; // STUDY THIS

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());

  if(ele_eta >= 2.5) return false;

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

  bool loose_barrel = fabs(ele.deltaEtaSuperClusterTrackAtVtx) < 0.007 &&
    fabs(ele.deltaPhiSuperClusterTrackAtVtx) < 0.15 &&
    fabs(ele.sigmaIetaIeta) < 0.01 &&
    ele.hcalOverEcalBc < 0.12 &&
    fabs(d0) < 0.02 &&
    fabs(dz) < 0.2 &&
    fabs(1/(ele.ecalEnergy) - 1/(ele.ecalEnergy/ele.eSuperClusterOverP)) < 0.05 &&
    ele_iso / ele.momentum.Pt() < 0.15 &&
    ele.passConversionVeto &&
    ele.nMissingHits <= 1;
      
  bool loose_endcap = fabs(ele.deltaEtaSuperClusterTrackAtVtx) < 0.009 &&
    ele.deltaPhiSuperClusterTrackAtVtx < 0.1 &&
    ele.sigmaIetaIeta < 0.03 &&
    ele.hcalOverEcalBc < 0.1 &&
    d0 < 0.02 &&
    dz < 0.2 &&
    fabs(1/(ele.ecalEnergy) - 1/(ele.ecalEnergy/ele.eSuperClusterOverP)) < 0.05 &&
    ((ele.momentum.Pt() > 20. && ele_iso / ele.momentum.Pt() < 0.15) || (ele.momentum.Pt() <= 20. && ele_iso / ele.momentum.Pt() < 0.1)) &&
    ele.passConversionVeto &&
    ele.nMissingHits <= 1;
  // mvaNonTrigV0 > 0.9 if we don't use simple cuts-based

  if(ele_eta <= 1.479) return loose_barrel;
  else if(ele_eta > 1.479 && ele_eta < 2.5) return loose_endcap;

  return false;

}

bool isJetVetoElectron(susy::Electron ele) {

  float iso = ele.chargedHadronIso + ele.neutralHadronIso + ele.photonIso;

  if(ele.passingMvaPreselection() &&
     fabs(ele.momentum.Eta()) < 2.6 &&
     ele.momentum.Pt() >= 15.0 &&
     iso / ele.momentum.Pt() < 0.2)
    return true;

  return false;
}

bool isVetoElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho, double d0, double dz) {

  if(ele.momentum.Pt() <= 20.) return false; // STUDY THIS

  if((int)ele.superClusterIndex >= (int)superClusters.size() || (int)ele.superClusterIndex < 0) return false;
  float ele_eta = fabs(superClusters[ele.superClusterIndex].position.Eta());

  if(ele_eta >= 2.5) return false;

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

  bool veto_barrel = fabs(ele.deltaEtaSuperClusterTrackAtVtx) < 0.007 &&
    fabs(ele.deltaPhiSuperClusterTrackAtVtx) < 0.8 &&
    fabs(ele.sigmaIetaIeta) < 0.01 &&
    ele.hcalOverEcalBc < 0.15 &&
    fabs(d0) < 0.04 &&
    fabs(dz) < 0.2 &&
    ele_iso / ele.momentum.Pt() < 0.15;
      
  bool veto_endcap = fabs(ele.deltaEtaSuperClusterTrackAtVtx) < 0.01 &&
    ele.deltaPhiSuperClusterTrackAtVtx < 0.7 &&
    ele.sigmaIetaIeta < 0.03 &&
    d0 < 0.04 &&
    dz < 0.2 &&
    ele_iso / ele.momentum.Pt() < 0.15;

  if(ele_eta <= 1.479) return veto_barrel;
  else if(ele_eta > 1.479 && ele_eta < 2.5) return veto_endcap;

  return false;
}

bool isMVAElectron(susy::Electron ele, vector<susy::SuperCluster> superClusters, double rho, double d0, double dz) {

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

  bool passes = ele_eta < 2.5 &&
    ele.momentum.Pt() > 20. &&
    ele.mvaNonTrig > 0.9 &&
    fabs(d0) < 0.02 &&
    fabs(dz) < 0.2 &&
    ele_iso / ele.momentum.Pt() < 0.15 &&
    ele.passConversionVeto &&
    ele.nMissingHits <= 1;

  return passes;
}
