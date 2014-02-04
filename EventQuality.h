#include "math.h"
#include "TMath.h"

#include "../src/SusyEvent.h"

bool good_Vtx(susy::Vertex vtx) {
  if(!vtx.isFake() &&
     vtx.ndof > 4 &&
     fabs((vtx.position).z()) < 24.0 &&
     fabs((vtx.position).Perp()) < 2.0
     )
    return true;
  else return false;
}

int GetNumberPV(susy::Event& event) {

  int nPV = 0;

  for(vector<susy::Vertex>::iterator iPV = event.vertices.begin(); iPV != event.vertices.end(); iPV++) {
    if(good_Vtx(*iPV)) nPV++;
  }

  return nPV;
}
