#include <iostream>
#include <fstream>

using namespace std;

#define nd 8

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
    double cotth1,cotth2,ecut,ptcut,mxmin2,mxmax2,thmax,thmin,qp2min,qp2max;
    int modcut;
    double mxmn,mxmx,q2mn,q2mx;
  } cuts_;
}

int main() {

  gmuini_();
  gmucha_();

  // Beam parameters
  beam_.inpe = beam_.inpp = 3500.;
  //beam_.pmod = 11;
  beam_.pmod = 2;
  beam_.emod = 2;

  // Vegas parameters
  vegpar_.nprn = 0;

  // Outgoing leptons kinematics
  cuts_.ecut = 0.;

  //ofstream cs("tmp/xsec_lpair_singleinelastic.dat");
  ofstream cs("tmp/xsec_lpair_elastic.dat");

  for (int i=0; i<100; i++) {
    cuts_.ptcut = 0.+i*0.1;
    gmubeg_();

    cs << cuts_.ptcut << "\t" << vgres_.s1 << "\t" << vgres_.s2 << endl;
    std::cout << "Pt > " << cuts_.ptcut << " GeV : xsec = " << vgres_.s1 << ", error = " << vgres_.s2 << std::endl;  
  }
  
  return 0;
}
