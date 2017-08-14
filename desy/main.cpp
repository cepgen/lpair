#include <iostream>

#include "external/utils.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "lpair.h"
#include "../commons/TreeEvent.h"

using namespace std;

int main(int argc, char* argv[]) {
  int one = 1;

  Timer tmr;

  //const int maxevts = 2.5e6;
  //const int maxevts = 2.5e6;
  const int maxevts = 1.e5;
  //const int maxevts = 5;

  TTree *t;

  zduini_();

cout << "zduini passed" << endl;

  if (vegpar_.iend<2) return 0;

  Lpair::TreeEvent ev;
  t = new TTree("h4444", "A TTree containing information from the events produced from LPAIR");
  ev.create( t );

  ev.xsect = vgres_.s1;
  ev.errxsect = vgres_.s2;
  for (int i=0; i<vegpar_.ngen; i++) {
    tmr.reset();

    zduevt_(&one);

    ev.gen_time = tmr.elapsed();    

    if (i%10000==0 and i>0) 
      cout << "[" << 100.*i/maxevts << "%] Generating event #" << i << " / " << maxevts << endl;
    ev.np = 0;
    for (int j=0; j<lujets_.n; j++) {
      ev.PID[ev.np] = lujets_.k[1][j];
      ev.parentid[ev.np] = lujets_.k[2][j];
      /*ev.daughterid1[ev.np] = lujets_.k[3][j];
        ev.daughterid2[ev.np] = lujets_.k[4][j];*/
      ev.status[ev.np] = lujets_.k[0][j];
      TLorentzVector part( lujets_.p[0][j], lujets_.p[1][j], lujets_.p[2][j], lujets_.p[3][j] );
      cout << "part mass=" << part.M() << endl;
      ev.pt[ev.np] = part.Pt();
      ev.eta[ev.np] = part.Eta();
      ev.phi[ev.np] = part.Phi();
      ev.E[ev.np] = part.E();
      ev.M[ev.np] = part.M();
      ev.charge[ev.np] = luchge_(lujets_.k[1][j])/3.;
      switch (j+1) {
        case 1: ev.role[ev.np] = 1; break;
        case 2: ev.role[ev.np] = 2; break;
        case 3: ev.role[ev.np] = 41; break;
        case 4: ev.role[ev.np] = 42; break;
        case 5: ev.role[ev.np] = 3; break;
        case 6: ev.role[ev.np] = 6; break;
        case 7: ev.role[ev.np] = 0; break;
        case 8: ev.role[ev.np] = 7; break;
        //case 9: ev.role[ev.np] = 5; break;
        default: ev.role[ev.np] = -1; break; // proton remnants
      }

      ev.np++;
    }
    t->Fill();
  }

  TString filename = "events.root";
  if (argc>2) { filename = TString(argv[2]); }
  t->SaveAs(filename);

  delete t;
  
  return 0;
}
