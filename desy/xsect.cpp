#include <iostream>
#include <fstream>
#include <cmath>
#include "lpair.h"

using namespace std;

double eta_to_theta(double eta) {
  return 2*atan(exp(eta))/acos(-1.)*180.;
}

int main() {

  gmuini_();
  gmucha_();

  // Beam parameters
  beam_.inpe = beam_.inpp = 3500.;
  beam_.pmod = 11;
  beam_.pmod = 2;
  beam_.emod = 2;

  // Vegas parameters
  vegpar_.iend = 1;
  vegpar_.nprn = 0;

  // Outgoing leptons kinematics
  cuts_.ecut = 0.;
  cuts_.thmin = eta_to_theta(-5.);
  cuts_.thmax = eta_to_theta( 5.);
  cuts_.q2mx = 1.e5;
  //cuts_.mxmx = 1000.;
  cuts_.ptcutmin = 3.;
  cout << "Theta in range [" << cuts_.thmin << ", " << cuts_.thmax << "]" << endl;

  cout << "modcut = " << cuts_.modcut << endl;

  ofstream cs("xsect_scan.dat");

  for (int i=0; i<100; i++) {
    cuts_.ptcutmin = 0.+i*0.5;
    gmubeg_();

    cs << cuts_.ptcutmin << "\t" << vgres_.s1 << "\t" << vgres_.s2 << endl;
    std::cout << "Pt > " << cuts_.ptcutmin << " GeV : xsec = " << vgres_.s1 << ", error = " << vgres_.s2 << std::endl;
  }

  return 0;
}
