#include <iostream>

#include "external/utils.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;
#define nd 10

extern "C" {
  void zduini_();
  void zduevt_(int* iwant);
  int luchge_(int&);
  extern struct {
    double inpe, inpp;
    int intge, intgp, gpdf, spdf, pmod, emod, ipair, nquark;
  } beam_;
  extern struct {
    int n, k[5][4000];
    float p[5][4000], v[5][4000];
  } lujets_;
  extern struct {
    double t1min,t1max;
    double t2min,t2max;
    double d3;
  } photons_;
  extern struct {
    double s1,s2,t1,t2;
  } extra_;
  extern struct {
    float s1,s2,s3,s4;
  } vgres_;
  extern struct {
    int ndim,ncvg,itmx,nprn,igraph,npoin,nprin,ntreat,ibeg,iend,ngen;
  } vegpar_;
  extern struct {
    double w,valtreat, x[nd], z[nd];
  } treatb_;
  extern struct {
    double u1,u2,v1,v2;
    double t11,t12,t21,t22;
  } peric_;
  extern struct {
    double tmx;
  } mykin_;
}

int main(int argc, char* argv[]) {
  int one = 1;

  Timer tmr;

  //const int maxevts = 2.5e6;
  //const int maxevts = 2.5e6;
  const int maxevts = 1.e5;
  //const int maxevts = 5;
  const int maxpart = 1000;

  int npart;
  double xsect, errxsect;
  double px[maxpart], py[maxpart], pz[maxpart], E[maxpart], M[maxpart], charge[maxpart];
  int PID[maxpart], parentid[maxpart], daughterid1[maxpart], daughterid2[maxpart], status[maxpart], role[maxpart];
  float gen_time;

  TTree *t;

  zduini_();

  if (vegpar_.iend<2) return 0;

  t = new TTree("h4444", "A TTree containing information from the events produced from LPAIR");

  t->Branch("ip", &npart, "npart/I");
  t->Branch("xsect", &xsect, "xsect/D");
  t->Branch("errxsect", &errxsect, "errxsect/D");
  t->Branch("total_time", &gen_time, "total_time/F");
  t->Branch("px", px, "px[npart]/D");
  t->Branch("py", py, "py[npart]/D");
  t->Branch("pz", pz, "pz[npart]/D");
  t->Branch("charge", charge, "charge[npart]/D");
  t->Branch("icode", PID, "PID[npart]/I");
  t->Branch("parent", parentid, "parent[npart]/I");
  t->Branch("daughter1", daughterid1, "daughter1[npart]/I");
  t->Branch("daughter2", daughterid2, "daughter2[npart]/I");
  t->Branch("status", status, "status[npart]/I");
  t->Branch("role", role, "role[npart]/I");
  t->Branch("E", E, "E[npart]/D");
  t->Branch("m", M, "M[npart]/D");

  xsect = vgres_.s1;
  errxsect = vgres_.s2;
  for (int i=0; i<vegpar_.ngen; i++) {
    tmr.reset();

    zduevt_(&one);

    gen_time = tmr.elapsed();    

    if (i%10000==0 and i>0) 
      cout << "[" << 100.*i/maxevts << "%] Generating event #" << i << " / " << maxevts << endl;
    npart = 0;
    for (int j=0; j<lujets_.n; j++) {
      PID[npart] = lujets_.k[1][j];
      parentid[npart] = lujets_.k[2][j];
      daughterid1[npart] = lujets_.k[3][j];
      daughterid2[npart] = lujets_.k[4][j];
      status[npart] = lujets_.k[0][j];
      px[npart] = lujets_.p[0][j];
      py[npart] = lujets_.p[1][j];
      pz[npart] = lujets_.p[2][j];
      E[npart] = lujets_.p[3][j];
      M[npart] = lujets_.p[4][j];
      charge[npart] = luchge_(lujets_.k[1][j])/3.;
      switch (j+1) {
        case 1: role[npart] = 1; break;
        case 2: role[npart] = 2; break;
        case 3: role[npart] = 41; break;
        case 4: role[npart] = 42; break;
        case 5: role[npart] = 3; break;
        case 6: role[npart] = 6; break;
        case 7: role[npart] = 0; break;
        case 8: role[npart] = 7; break;
        //case 9: role[npart] = 5; break;
        default: role[npart] = -1; break; // proton remnants
      }

      npart++;
    }
    t->Fill();
  }

  TString filename = "events.root";
  if (argc>2) { filename = TString(argv[2]); }
  t->SaveAs(filename);

  delete t;
  
  return 0;
}
