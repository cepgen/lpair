#include <iostream>

#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

extern "C" {
  void zduini_();
  void zduevt_(int* iwant);
  extern struct {
    double inpe, inpp;
    int intge, intgp, gpdf, spdf, pmod, emod, ipair, nquark;
  } beam_;
  void printlhe_(int* mode);
  extern struct {
    int n, npad, k[5][4000];
    double p[5][4000], v[5][4000];
  } pyjets_;
}

Int_t main() {
  Int_t one = 1;
  /*Int_t two = 2;
  Int_t three = 3;*/

  Int_t parentid, firstdaughterid, lastdaughterid;

  const Int_t maxpart = 100;

  Int_t npart;
  Double_t eta[maxpart], phi[maxpart], rapidity[maxpart], px[maxpart], py[maxpart], pz[maxpart], pt[maxpart], E[maxpart], M[maxpart];
  Int_t PID[maxpart], isstable[maxpart];

  TTree *t;
  TLorentzVector *mom;

  t = new TTree("h4444", "A TTree containing information from the events produced from LPAIR");
  mom = new TLorentzVector();

  t->Branch("npart", &npart, "npart/I");
  t->Branch("eta", eta, "eta[npart]/D");
  t->Branch("phi", phi, "phi[npart]/D");
  t->Branch("rapidity", rapidity, "rapidity[npart]/D");
  t->Branch("px", px, "px[npart]/D");
  t->Branch("py", py, "py[npart]/D");
  t->Branch("pz", pz, "pz[npart]/D");
  t->Branch("pt", pt, "pt[npart]/D");
  t->Branch("PID", PID, "PID[npart]/I");
  t->Branch("stable", isstable, "isstable[npart]/I");
  t->Branch("E", E, "E[npart]/D");
  t->Branch("M", M, "M[npart]/D");

  zduini_();
  //printlhe_(&one);
  for (Int_t i=0; i<1e2; i++) {
    zduevt_(&one);
    npart = pyjets_.n;
    for (Int_t p=0; p<pyjets_.n; p++) {
      PID[p] = pyjets_.k[1][p];
      parentid = pyjets_.k[2][p];
      firstdaughterid = pyjets_.k[3][p];
      lastdaughterid = pyjets_.k[4][p];
      isstable[p] = pyjets_.v[4][p]==0.;
      px[p] = pyjets_.p[0][p];
      py[p] = pyjets_.p[1][p];
      pz[p] = pyjets_.p[2][p];
      E[p] = pyjets_.p[3][p];
      M[p] = pyjets_.p[4][p];
      pt[p] = sqrt(pow(px[p],2)+pow(py[p],2));
      
      cout << "--> " << p << " = " << pyjets_.p[0][p] << "\t" << pyjets_.p[1][p] << "\t" << pyjets_.p[2][p] << "\t" << pyjets_.p[3][p] << "\t" << pyjets_.p[4][p] << endl;
      mom->SetXYZM(px[p], py[p], pz[p], M[p]);
      eta[p] = mom->PseudoRapidity();
      rapidity[p] = mom->Rapidity();
      phi[p] = mom->Phi();
      //cout << "    " << pyjets_.v[0][p] << "\t" << pyjets_.v[1][p] << "\t" << pyjets_.v[2][p] << "\t" << pyjets_.v[3][p] << "\t" << pyjets_.v[4][p] << endl;
      //cout << "    " << pyjets_.k[0][p] << "\t" << pyjets_.k[1][p] << "\t" << pyjets_.k[2][p] << "\t" << pyjets_.k[3][p] << "\t" << pyjets_.k[4][p] << endl;
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
