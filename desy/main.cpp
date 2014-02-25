#include <iostream>

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
  void printlhe_(int* mode);
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
    int ndim,ncvg,itmx,nprn,igraph,npoin,nprin,ntreat,ibeg,iend;
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

Int_t main() {
  Int_t one = 1;
  Int_t two = 2;
  Int_t three = 3;

  Int_t firstdaughterid, lastdaughterid;

  const Int_t maxpart = 1000;

  Int_t npart;
  Double_t xsect, errxsect;
  Double_t eta[maxpart], phi[maxpart], rapidity[maxpart], px[maxpart], py[maxpart], pz[maxpart], pt[maxpart], E[maxpart], M[maxpart], charge[maxpart];
  Int_t PID[maxpart], isstable[maxpart], parentid[maxpart];
  Double_t mx;
  Double_t t1, t1min, t1max, t2, t2min, t2max;
  Double_t valtreat, wtreat;
  Double_t s1, s2;
  Double_t d3;
  Double_t u1, u2, v1, v2;
  Double_t t11, t12, t21, t22;

  TTree *t;
  TLorentzVector *mom;

  zduini_();

  Int_t ndim = vegpar_.ndim;
  Double_t xtreat[ndim], ztreat[ndim];

  t = new TTree("h4444", "A TTree containing information from the events produced from LPAIR");
  mom = new TLorentzVector();

  t->Branch("ip", &npart, "npart/I");
  t->Branch("xsect", &xsect, "xsect/D");
  t->Branch("errxsect", &errxsect, "errxsect/D");
  t->Branch("Eta", eta, "eta[npart]/D");
  t->Branch("phi", phi, "phi[npart]/D");
  t->Branch("rapidity", rapidity, "rapidity[npart]/D");
  t->Branch("px", px, "px[npart]/D");
  t->Branch("py", py, "py[npart]/D");
  t->Branch("pz", pz, "pz[npart]/D");
  t->Branch("pt", pt, "pt[npart]/D");
  t->Branch("charge", charge, "charge[npart]/D");
  t->Branch("icode", PID, "PID[npart]/I");
  t->Branch("parent", parentid, "parent[npart]/I");
  t->Branch("stable", isstable, "stable[npart]/I");
  t->Branch("E", E, "E[npart]/D");
  t->Branch("m", M, "M[npart]/D");
  t->Branch("MX", &mx, "MX/D");
  t->Branch("t1", &t1, "t1/D");
  t->Branch("t1min", &t1min, "t1min/D");
  t->Branch("t1max", &t1max, "t1max/D");
  t->Branch("t2", &t2, "t2/D");
  t->Branch("t2min", &t2min, "t2min/D");
  t->Branch("t2max", &t2max, "t2max/D");
  t->Branch("s1", &s1, "s1/D");
  t->Branch("s2", &s2, "s2/D");
  t->Branch("d3", &d3, "d3/D");
  t->Branch("ndim", &ndim, "ndim/I");
  t->Branch("wtreat", &wtreat, "wtreat/D");
  t->Branch("xtreat", xtreat, "xtreat[ndim]/D");
  t->Branch("ztreat", ztreat, "ztreat[ndim]/D");
  t->Branch("valtreat", &valtreat, "valtreat/D");
  t->Branch("u1", &u1, "u1/D");
  t->Branch("u2", &u2, "u2/D");
  t->Branch("v1", &v1, "u1/D");
  t->Branch("v2", &v2, "u2/D");
  t->Branch("t11", &t11, "t11/D");
  t->Branch("t12", &t12, "t12/D");
  t->Branch("t21", &t21, "t21/D");
  t->Branch("t22", &t22, "t22/D");

  //printlhe_(&one);

  xsect = vgres_.s1;
  errxsect = vgres_.s2;
  for (Int_t i=0; i<1e5; i++) {
    zduevt_(&one);
    if (i%10000==0 and i>0) {
      cout << "Generating event #" << i << endl;
    }
    npart = 0;
    mx = mykin_.tmx;
    t1 = extra_.t1;
    t1min = photons_.t1min;
    t1max = photons_.t1max;
    t2 = extra_.t2;
    t2min = photons_.t2min;
    t2max = photons_.t2max;
    s1 = extra_.s1;
    s2 = extra_.s2;
    d3 = photons_.d3;
    wtreat = treatb_.w;
    u1 = peric_.u1;
    u2 = peric_.u2;
    v1 = peric_.v1;
    v2 = peric_.v2;
    t11 = peric_.t11;
    t12 = peric_.t12;
    t21 = peric_.t21;
    t22 = peric_.t22;
    for (Int_t d=0; d<ndim; d++) {
      xtreat[d] = treatb_.x[d];
      ztreat[d] = treatb_.z[d];
    }
    //npart = pyjets_.n;
    for (Int_t p=0; p<lujets_.n; p++) {
      firstdaughterid = lujets_.k[3][p];
      lastdaughterid = lujets_.k[4][p];
      if (firstdaughterid!=0 or lastdaughterid!=0) continue;
      PID[npart] = lujets_.k[1][p];
      parentid[npart] = lujets_.k[2][p];
      //cout << p << "\t" << PID[npart] << "\t" << firstdaughterid << "\t" << lastdaughterid << endl;
      //isstable[npart] = lujets_.v[4][p]==0.;
      isstable[npart] = lujets_.k[0][p]==1;
      px[npart] = lujets_.p[0][p];
      py[npart] = lujets_.p[1][p];
      pz[npart] = lujets_.p[2][p];
      E[npart] = lujets_.p[3][p];
      M[npart] = lujets_.p[4][p];
      pt[npart] = sqrt(pow(px[npart],2)+pow(py[npart],2));
      charge[npart] = luchge_(lujets_.k[1][p])/3.;
      
      //cout << "--> PDG = " << PID[npart] << ", " << p << " = " << lujets_.p[0][p] << "\t" << lujets_.p[1][p] << "\t" << lujets_.p[2][p] << "\t" << lujets_.p[3][p] << "\t" << lujets_.p[4][p] << endl;
      mom->SetXYZM(px[npart], py[npart], pz[npart], M[npart]);
      if (pt[npart]==0.) eta[npart] = 9999.;
      else eta[npart] = mom->PseudoRapidity();
      rapidity[npart] = mom->Rapidity();
      phi[npart] = mom->Phi();
      npart++;
      //cout << "    " << lujets_.v[0][p] << "\t" << lujets_.v[1][p] << "\t" << lujets_.v[2][p] << "\t" << lujets_.v[3][p] << "\t" << lujets_.v[4][p] << endl;
      //cout << "    " << lujets_.k[0][p] << "\t" << lujets_.k[1][p] << "\t" << lujets_.k[2][p] << "\t" << lujets_.k[3][p] << "\t" << lujets_.k[4][p] << endl;
    }
    t->Fill();
    //printlhe_(&two);
  }
  //printlhe_(&three);

  t->SaveAs("events.root");

  delete mom;
  delete t;
  
  return 0;
}
