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

Int_t main() {
  Int_t one = 1;

  Timer tmr;

  const Int_t maxpart = 1000;

  Int_t npart;
  Double_t xsect, errxsect;
  Double_t eta[maxpart], phi[maxpart], rapidity[maxpart], px[maxpart], py[maxpart], pz[maxpart], pt[maxpart], p[maxpart], E[maxpart], M[maxpart], charge[maxpart];
  Int_t PID[maxpart], isstable[maxpart], parentid[maxpart], daughterid1[maxpart], daughterid2[maxpart], status[maxpart];
  Double_t mx, q2[2];
  Double_t t1, t1min, t1max, t2, t2min, t2max;
  Double_t valtreat, wtreat;
  Double_t s1, s2;
  Double_t d3;
  Double_t u1, u2, v1, v2;
  Double_t t11, t12, t21, t22;
  Double_t q2calc;
  Float_t gen_time;

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
  t->Branch("total_time", &gen_time, "total_time/F");
  t->Branch("Eta", eta, "eta[npart]/D");
  t->Branch("phi", phi, "phi[npart]/D");
  t->Branch("rapidity", rapidity, "rapidity[npart]/D");
  t->Branch("px", px, "px[npart]/D");
  t->Branch("py", py, "py[npart]/D");
  t->Branch("pz", pz, "pz[npart]/D");
  t->Branch("p", p, "p[npart]/D");
  t->Branch("pt", pt, "pt[npart]/D");
  t->Branch("charge", charge, "charge[npart]/D");
  t->Branch("icode", PID, "PID[npart]/I");
  t->Branch("parent", parentid, "parent[npart]/I");
  t->Branch("daughter1", daughterid1, "daughter1[npart]/I");
  t->Branch("daughter2", daughterid2, "daughter2[npart]/I");
  t->Branch("stable", isstable, "stable[npart]/I");
  t->Branch("status", status, "status[npart]/I");
  t->Branch("E", E, "E[npart]/D");
  t->Branch("m", M, "M[npart]/D");
  t->Branch("MX", &mx, "MX/D");
  t->Branch("t1", &t1, "t1/D");
  t->Branch("q2", q2, "q2[2]/D");

  xsect = vgres_.s1;
  errxsect = vgres_.s2;
  for (Int_t i=0; i<vegpar_.ngen; i++) {
    tmr.reset();

    zduevt_(&one);

    gen_time = tmr.elapsed();    

    if (i%10000==0 and i>0) 
      cout << "Generating event #" << i << endl;
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
    for (Int_t j=0; j<lujets_.n; j++) {
      PID[npart] = lujets_.k[1][j];
      parentid[npart] = lujets_.k[2][j];
      daughterid1[npart] = lujets_.k[3][j];
      daughterid2[npart] = lujets_.k[4][j];
      status[npart] = lujets_.k[0][j];
      isstable[npart] = lujets_.k[0][j]==1;
      px[npart] = lujets_.p[0][j];
      py[npart] = lujets_.p[1][j];
      pz[npart] = lujets_.p[2][j];
      E[npart] = lujets_.p[3][j];
      M[npart] = lujets_.p[4][j];
      pt[npart] = sqrt(pow(px[npart],2)+pow(py[npart],2));
      p[npart] = sqrt(pow(pt[npart],2)+pow(pz[npart],2));
      charge[npart] = luchge_(lujets_.k[1][j])/3.;

      if (PID[npart]==22 and status[npart]==21) {
	q2calc = pow(E[npart], 2)-pow(p[npart], 2);
	if (parentid[npart]==1) q2[0] = q2calc; // first proton's photon
	else if (parentid[npart]==2) q2[1] = q2calc; // second proton's photon
	//cout << "Photon found : parent = " << parentid[npart] << ", q2 = " << q2calc << endl;
      }
      //cout << "--> " << PID[npart] << " has status " << status[npart] << " and parent " << parentid[npart] << " (M=" << M[npart] << ")" << endl;
      
      //cout << "--> PDG = " << PID[npart] << ", " << j << " = " << lujets_.p[0][j] << "\t" << lujets_.p[1][j] << "\t" << lujets_.p[2][j] << "\t" << lujets_.p[3][j] << "\t" << lujets_.p[4][j] << endl;
      mom->SetXYZM(px[npart], py[npart], pz[npart], M[npart]);
      if (pt[npart]==0.) eta[npart] = 9999.;
      else eta[npart] = mom->PseudoRapidity();
      rapidity[npart] = mom->Rapidity();
      phi[npart] = mom->Phi();
      npart++;
      //cout << "    " << lujets_.v[0][j] << "\t" << lujets_.v[1][j] << "\t" << lujets_.v[2][j] << "\t" << lujets_.v[3][j] << "\t" << lujets_.v[4][j] << endl;
      //cout << "    " << lujets_.k[0][j] << "\t" << lujets_.k[1][j] << "\t" << lujets_.k[2][j] << "\t" << lujets_.k[3][j] << "\t" << lujets_.k[4][j] << endl;
    }
    t->Fill();
  }

  t->SaveAs("events.root");

  delete mom;
  delete t;
  
  return 0;
}
