using namespace std;

bool isGoodJet(susy::PFJet jet, TLorentzVector corrP4) {
  
  bool isGood = false;

  if((corrP4.Pt() > 30.0) && (fabs(corrP4.Eta()) < 2.6)) { // STUDY BOTH OF THESE

    if((jet.neutralHadronEnergy/jet.momentum.Energy() < 0.99) &&
       (jet.neutralEmEnergy/jet.momentum.Energy() < 0.99) &&
       ((unsigned int)jet.nConstituents > 1) &&
       jet.passPuJetIdTight(susy::kPUJetIdFull)) {

      if(fabs(corrP4.Eta()) < 2.4) {
	if((jet.chargedHadronEnergy > 0.0) &&
	   ((int)jet.chargedMultiplicity > 0) &&
	   (jet.chargedEmEnergy/jet.momentum.Energy() < 0.99))
	  isGood = true;
      }

      else isGood = true;
    }
  }

  return isGood;
}

float deltaR_jetLep(TLorentzVector& p1, TLorentzVector& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

bool JetOverlapsElectron(TLorentzVector corrP4, vector<susy::Electron*> isoEles, vector<susy::Electron*> looseEles) {

  bool same_ele = false;

  for(unsigned int i = 0; i < isoEles.size(); i++) {
    if(deltaR_jetLep(corrP4, isoEles[i]->momentum) < 0.5) {
      same_ele = true;
      break;
    }
  }

  for(unsigned int i = 0; i < looseEles.size(); i++) {
    if(deltaR_jetLep(corrP4, looseEles[i]->momentum) < 0.5) {
      same_ele = true;
      break;
    }
  }

  return same_ele;
}

bool JetOverlapsMuon(TLorentzVector corrP4, vector<susy::Muon*> isoMuons, vector<susy::Muon*> looseMuons) {

  bool same_muon = false;

  for(unsigned int i = 0; i < isoMuons.size(); i++) {
    if(deltaR_jetLep(corrP4, isoMuons[i]->momentum) < 0.5) {
      same_muon = true;
      break;
    }
  }

  for(unsigned int i = 0; i < looseMuons.size(); i++) {
    if(deltaR_jetLep(corrP4, looseMuons[i]->momentum) < 0.5) {
      same_muon = true;
      break;
    }
  }
  
  return same_muon;
}
