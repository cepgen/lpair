#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "lpair.h"

#include "../commons/TreeInfo.h"
#include "../commons/Timer.h"

using namespace std;

int main( int argc, char* argv[] )
{
  int one = 1;

  Timer tmr;

  const char* filename = ( argc > 1 ) ? argv[1] : "events.root";

  TFile f( filename, "recreate" );

  zduini_();

  if ( vegpar_.iend < 2 )
    return 0;

  lpair::TreeRun run;
  run.create();
  run.xsect = vgres_.s1;
  run.errxsect = vgres_.s2;
  run.sqrt_s = beam_.inpe+beam_.inpp;
  run.fill();

  lpair::TreeEvent ev;
  ev.create( new TTree( "h4444", "A TTree containing information from the events produced from LPAIR" ) );

  for ( int i = 0; i < vegpar_.ngen; ++i ) {
    tmr.reset();

    zduevt_(&one);

    ev.gen_time = tmr.elapsed();

    if ( i % 10000 == 0 && i > 0 )
      cout << "[" << 100.*i/vegpar_.ngen << "%] Generating event #" << i << " / " << vegpar_.ngen << endl;
    ev.np = 0;
    ev.momentum.reserve( lujets_.n );
    for ( int j = 0; j < lujets_.n; ++j ) {
      ev.pdg_id[ev.np] = lujets_.k[1][j];
      ev.parent1[ev.np] = lujets_.k[2][j];
      ev.parent2[ev.np] = lujets_.k[3][j];
      ev.status[ev.np] = lujets_.k[0][j];
      ev.momentum[ev.np].SetPxPyPzE( lujets_.p[0][j], lujets_.p[1][j], lujets_.p[2][j], lujets_.p[3][j] );
      ev.pt[ev.np] = hypot( lujets_.p[0][j], lujets_.p[1][j] );
      const double p = hypot( ev.pt[ev.np], lujets_.p[2][j] );
      ev.eta[ev.np] = ( p != 0 ) ? atanh( lujets_.p[2][j]/p ) : 0;
      ev.phi[ev.np] = ( lujets_.p[0][j] * lujets_.p[1][j] != 0. ) ? atan2( lujets_.p[1][j], lujets_.p[0][j] ) : 0.;
      ev.E[ev.np] = lujets_.p[3][j];
      ev.m[ev.np] = lujets_.p[4][j];
      ev.charge[ev.np] = luchge_( lujets_.k[1][j] )/3.;
      switch ( j+1 ) {
        case 1: ev.role[ev.np] = 1; break;
        case 2: ev.role[ev.np] = 2; break;
        case 3: ev.role[ev.np] = 41; break;
        case 4: ev.role[ev.np] = 42; break;
        case 5: ev.role[ev.np] = 3; break;
        case 6: ev.role[ev.np] = 6; break;
        case 7: ev.role[ev.np] = 0; break;
        case 8: ev.role[ev.np] = 7; break;
        case 9: ev.role[ev.np] = 5; break;
        default: ev.role[ev.np] = -1; break; // proton remnants
      }

      ev.np++;
    }
    ev.fill();
  }

  f.Write();
  f.Close();

  return 0;
}
