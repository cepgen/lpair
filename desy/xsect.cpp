#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

extern "C" {
  //void zduini_();
  //void zduevt_(int* iwant);
  void gmuini_();
  void gmucha_();
  void gmubeg_();

  extern struct {
    double inpe, inpp;
    int intge, intgp, gpdf, spdf, pmod, emod, ipair, nquark;
  } beam_;

  extern struct {
    int ndim,ncvg,itmx,nprn,igraph,npoin,nprin,ntreat,ibeg,iend;
  } vegpar_;

  extern struct {
    float s1,s2,s3,s4;
  } vgres_;

  extern struct {
    double cotth1,cotth2,ecut,ptcutmin,ptcutmax,mxmin2,mxmax2;
    float thmax,thmin;
    double qp2min,qp2max;
    int modcut;
    float mxmn,mxmx,q2mn,q2mx;
  } cuts_;
}

int main() {

  gmuini_();
  gmucha_();

  // Beam parameters
  beam_.inpe = beam_.inpp = 3500.;
  beam_.pmod = 11;
  //beam_.pmod = 2;
  beam_.emod = 2;

  // Vegas parameters
  vegpar_.iend = 1;
  vegpar_.nprn = 0;

  // Outgoing leptons kinematics
  cuts_.ecut = 0.;
  cuts_.thmin = 2*atan(exp(-5.))/acos(-1.)*180.;
  cuts_.thmax = 2*atan(exp( 5.))/acos(-1.)*180.;
  cuts_.q2mx = 1.e5;
  //cuts_.mxmx = 1000.;
  cuts_.ptcutmin = 3.;
  cout << "Theta in range [" << cuts_.thmin << ", " << cuts_.thmax << "]" << endl;

  ofstream cs("xsect_scan_lpair_singleinelastic.dat");
  //ofstream cs("xsect_scan_lpair_elastic.dat");

  for (int i=0; i<50; i++) {
    cuts_.ptcutmin = 0.+i*0.5;
    //cuts_.q2mx = 2500.+i*2500;
    //cuts_.mxmx = 100.+i*100;
    gmubeg_();

    //cs << cuts_.q2mx << "\t" << vgres_.s1 << "\t" << vgres_.s2 << endl;
    //std::cout << "Q2<" << cuts_.q2mx << " GeV**2: xsec = " << vgres_.s1 << ", error = " << vgres_.s2 << std::endl;  
    //cs << cuts_.mxmx << "\t" << vgres_.s1 << "\t" << vgres_.s2 << endl;
    //std::cout << "MX<" << cuts_.mxmx << " GeV**2: xsec = " << vgres_.s1 << ", error = " << vgres_.s2 << std::endl;  
    cs << cuts_.ptcutmin << "\t" << vgres_.s1 << "\t" << vgres_.s2 << endl;
    std::cout << "pt>" << cuts_.ptcutmin << " GeV: xsec = " << vgres_.s1 << ", error = " << vgres_.s2 << std::endl;  
  }
  
  return 0;
}
