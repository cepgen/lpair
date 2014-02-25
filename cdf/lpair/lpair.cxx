#include <iostream>
#include <fstream>

using namespace std;

extern "C" {
  //void zduini_();
  //void zduevt_(int* iwant);
  void fileini_();
  void integrate_();
  void generate_(int& nevents);
  void pawfil1_();
  void doublefragmentation_();

  extern struct {
    int ipar[20];
    double lpar[20];
  } datapar_;

  extern struct {
    double s1,s2,s3,s4;
  } result_;

  extern struct {
    bool accepted;
  } event_;

  extern struct {
    double e, e1, e2, e3, e4, e5;
    double p, p3, p4, p5;
    double ct3, st3, ct4, st4, ct5, st5, cp3, sp3, cp5, sp5;
  } variab_;

  extern struct {
    double al3, al4, be4, be5, de3, de5;
    double pp3, pp4, pp5;
  } variac_;

  extern struct {
    double px3, py3, pz3;
    double px4, py4, pz4;
    double px5, py5, pz5;
    double px6, py6, pz6;
    double px7, py7, pz7;
  } mygenz_;

  /*
      common/inpu/me,mu,ebeam,const,sq
      common/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5
     1                                         ,st5,cp3,sp3,cp5,sp5
      common/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
      common/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
      common/lplot/xl(10),v1(2),v2(2),av(10)
      common/extra/s1,s2,t1,t2
      common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
      common/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
      common/civita/epsi,g5,g6,a5,a6,bb
      common/dotps/q1dq,q1dq2,w6
      common/tell/nn
      common/cuts/angcut,encut,etacut
      common/mygenz/px3,py3,pz3,px4,py4,pz4,px5,py5,pz5,
     +              px6,py6,pz6,px7,py7,pz7
   */
}

int main() {
  // Number of events to generate
  const int nevent = 2;
  int ev = 1;
  int i;

  fileini_();

  // Beam parameters
  datapar_.ipar[4] = 7;
  datapar_.lpar[2] = 3500.;

  // Outgoing leptons kinematics
  //datapar_.lpar[4] = 999.; // eta cut
  datapar_.lpar[5] = 0.; // energy cut
  datapar_.lpar[6] = 5.; // pt cut

  integrate_();

  std::cout << "Pt > " << datapar_.lpar[6] << " GeV : xsec = " << result_.s1 << ", error = " << result_.s2 << std::endl;  

  i = 0;
  do {
    generate_(ev);
    if (event_.accepted) i++;
    std::cout << "Event " << i << " generated !" << std::endl;
    pawfil1_();
    doublefragmentation_();
    /*// first outgoing proton
    double px3, py3, pz3;
    px3 = variac_.pp3*variab_.cp3;
    py3 = variac_.pp3*variab_.sp3;
    pz3 = variab_.p3*variab_.ct3;

    std::cout << px3 << "\t" << mygenz_.px3 << std::endl;
    std::cout << py3 << "\t" << mygenz_.py3 << std::endl;
    std::cout << pz3 << "\t" << mygenz_.pz3 << std::endl;*/

  } while (i<nevent);

  return 0;
}
