#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "lpair.h"
#include "../commons/TreeInfo.h"
#include "../commons/Timer.h"

using namespace std;

int main( int argc, char* argv[] )
{
  const char* filename = ( argc > 1 ) ? argv[1] : "events.root";
  // Number of events to generate
  const unsigned long long nevents = ( argc > 2 ) ? atoll( argv[2] ) : 1e5;

  Timer tmr;

  TFile f( filename, "recreate" );

  fileini_();

  // Beam parameters
  /*datapar_.ipar[4] = 9;
  datapar_.lpar[2] = 3500.;

  // Outgoing leptons kinematics
  datapar_.lpar[4] = 2.5; // eta cut
  //datapar_.lpar[4] = 3.1313; // eta cut
  datapar_.lpar[5] = 0.; // energy cut
  datapar_.lpar[6] = 5.; // pt cut*/

  integrate_();

  lpair::TreeRun run;
  run.create();
  run.xsect = result_.s1;
  run.errxsect = result_.s2;
  run.sqrt_s = 2.*datapar_.lpar[2];
  run.fill();

  cout << "Pt > " << datapar_.lpar[6] << " GeV :" << endl
	  << "  xsec  = " << result_.s1 << endl
	  << "  error = " << result_.s2 << endl;

  lpair::TreeEvent ev;
  ev.create( new TTree( "h4444", "A TTree containing information from the events produced from LPAIR (CDF)" ) );

  unsigned long long i = 0;
  int evt = 1;
  do {
    tmr.reset();

    generate_( evt );

    if ( !event_.accepted )
      continue;

    if ( i % 10000 == 0 && i > 0 )
      cout << "[" <<( 100.*i/nevents ) << "%] Generating event #" << i << " / " << nevents << endl;

    ev.gen_time = tmr.elapsed();

    fragmentation_();
    ev.tot_time = tmr.elapsed();

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
        case 6: break; // dilepton system
        case 7: ev.role[ev.np] = 5; break;
        case 8: ev.role[ev.np] = 6; break;
        case 9: ev.role[ev.np] = 7; break;
        /*case 10: role[ev.np] = 3; break;
        case 11: role[ev.np] = 3; break;
        case 12: role[ev.np] = 5; break;
        case 13: role[ev.np] = 5; break;*/
        default: ev.role[ev.np] = -1; break; // FIXME need to figure out the origin of each remnant
      }

      ev.np++;
    }
    ev.fill();
  } while ( i<nevents );

  f.Write();
  f.Close();

  return 0;
}
