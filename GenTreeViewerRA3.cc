#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>

#include "TString.h"
#include "../src/SusyEvent.h"

using namespace std;

struct PNode {
  PNode() : mother(0), daughters(0), pdgId(0), status(0), mass(0.) {}
  ~PNode();
  PNode& operator=(PNode const& _rhs)
  {
    mother = _rhs.mother;
    daughters = _rhs.daughters;
    pdgId = _rhs.pdgId;
    status = _rhs.status;
    mass = _rhs.mass;
    pt = _rhs.pt;
    return *this;
  }
  PNode* mother;
  vector<PNode*> daughters;
  int pdgId;
  int status;
  double mass;
  double pt;
  string print(vector<bool>&, bool = false);
  bool hasMother() { return mother != 0; }
};

PNode::~PNode()
{
  for(unsigned iD(0); iD < daughters.size(); ++iD){
    delete daughters[iD];
  }
}

string
PNode::print(vector<bool>& _isLastChildAt, bool _showMass/* = false*/)
{
  using namespace std;

  stringstream ss;
  ss << "+" << setw(8) << pdgId;
  if(_showMass)
    ss << " [" << setw(6) << fixed << setprecision(2) << mass << "]";
  ss << " (" << setw(5) << fixed << setprecision(1) << pt << ") " << status << " ";
  for(unsigned iD(0); iD < daughters.size(); iD++){
    if(iD > 0){
      ss << endl;
      for(unsigned i(0); i < _isLastChildAt.size(); i++){
        ss << (_isLastChildAt[i] ? " " : "|") << "                   ";
        if(_showMass)
          ss << "         ";
      }
    }
    _isLastChildAt.push_back(iD == daughters.size() - 1);
    ss << daughters[iD]->print(_isLastChildAt, _showMass);
    _isLastChildAt.pop_back();
  }

  return ss.str();
}

PNode*
setDaughters(short _index, vector<susy::Particle>& _particles, float _minPt)
{
  susy::Particle& gen(_particles[_index]);
  if(gen.status == 1 && gen.momentum.Pt() < _minPt) return 0;

  PNode* node(new PNode);
  node->pdgId = gen.pdgId;
  node->status = gen.status;
  node->mass = gen.momentum.M();
  node->pt = gen.momentum.Pt();

  for(unsigned iP(0); iP < _particles.size(); ++iP){
    if(_particles[iP].motherIndex == _index){
      PNode* daughter(setDaughters(iP, _particles, _minPt));
      daughter->mother = node;
      node->daughters.push_back(daughter);
    }
  }

  return node;
}

class GenTreeViewerRA3 {
 public:
  GenTreeViewerRA3(TTree&);
  virtual ~GenTreeViewerRA3() { ; };

  virtual void viewGenTreeRA3();

  void SetProcessNEvents(int v) { processNEvents = v; }
  void SetShowMass(bool v) { showMass = v; }
  void SetMinPt(float v) { minPt = v; }

 protected:
  susy::Event event;
  TTree *fTree;
  int processNEvents;
  bool showMass;
  float minPt;

};

GenTreeViewerRA3::GenTreeViewerRA3(TTree& tree) :
  event(),
  fTree(&tree),
  processNEvents(-1),
  showMass(false),
  minPt(2.)
{
  event.setInput(tree);
}
  
void GenTreeViewerRA3::viewGenTreeRA3() {

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {
  
   vector<PNode*> rootNodes;

   vector<susy::Particle>& particles(event.genParticles);

   for(short iP(0); iP < short(particles.size()); ++iP){
     susy::Particle& gen(particles[iP]);

     if(gen.motherIndex == -1){
       PNode* node(setDaughters(iP, particles, minPt));
       rootNodes.push_back(node);
     }
   }

   cout << "--- CLEANED DECAY TREE ---" << endl << endl;
   for(unsigned iN(0); iN < rootNodes.size(); iN++){
     vector<bool> isLastChildAt(1, true);
     cout << rootNodes[iN]->print(isLastChildAt, showMass);
     cout << endl;

     delete rootNodes[iN];
   }

   cout << endl;
  }

}
