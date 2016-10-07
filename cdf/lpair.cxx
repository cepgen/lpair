#include <iostream>
#include <fstream>

#include "utils.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

extern "C" {
  //void zduini_();
  //void zduevt_(int* iwant);
  void fileini_();
  void integrate_();
  void generate_(int& nevents);
  void fragmentation_();
  int luchge_(int&);

  extern struct {
    int ipar[20];
    double lpar[20];
  } datapar_;

  extern struct {
    double s1,s2,s3,s4;
  } result_;

  extern struct {
    bool accepted;
    int ndim;
    double x[10];
    int leppdg;
  } event_;

  extern struct {
    int n, k[5][4000];
    float p[5][4000],v[5][4000];
  } lujets_;
  extern struct {
    double mx1, mx2;
    double wx1, wx2;
    double w1, w6, w3, w8;
  } remnts_;
}

int main(int argc, char* argv[]) {
  // Number of events to generate
  const int nevent = 1e5;
  //const int nevent = 1e4;
  int ev = 1;
  int i;

  Timer tmr;

  fileini_();

  string filename = (argc>2) ? argv[2] : "events.root";

  // Beam parameters
  /*datapar_.ipar[4] = 9;
  datapar_.lpar[2] = 3500.;

  // Outgoing leptons kinematics
  datapar_.lpar[4] = 2.5; // eta cut
  //datapar_.lpar[4] = 3.1313; // eta cut
  datapar_.lpar[5] = 0.; // energy cut
  datapar_.lpar[6] = 5.; // pt cut*/

  integrate_();

  std::cout << "Pt > " << datapar_.lpar[6] << " GeV :" << std::endl
	    << "  xsec  = " << result_.s1 << std::endl
	    << "  error = " << result_.s2 << std::endl;  

  const int maxpart = 1000;

  double xsect, errxsect;
  int npart, ndim;
  int role[maxpart];
  double pt[maxpart], eta[maxpart], phi[maxpart], rapidity[maxpart];
  double E[maxpart], M[maxpart], charge[maxpart];
  int PID[maxpart], isstable[maxpart], status[maxpart], parentid[maxpart];
  double mx1, mx2;
  TLorentzVector *mom;
  TTree *t;
  float time_gen, time_tot;

  t = new TTree("h4444", "A TTree containing information from the events produced from LPAIR (CDF)");
  mom = new TLorentzVector();

  t->Branch("xsect", &xsect, "xsect/D");
  t->Branch("errxsect", &errxsect, "errxsect/D");
  t->Branch("npart", &npart, "npart/I");
  t->Branch("rapidity", rapidity, "rapidity[npart]/D");
  t->Branch("pt", pt, "pt[npart]/D");
  t->Branch("eta", eta, "eta[npart]/D");
  t->Branch("phi", phi, "phi[npart]/D");
  t->Branch("charge", charge, "charge[npart]/D");
  t->Branch("icode", PID, "icode[npart]/I");
  t->Branch("parent", parentid, "parent[npart]/I");
  t->Branch("stable", isstable, "stable[npart]/I");
  t->Branch("status", status, "status[npart]/I");
  t->Branch("role", role, "role[npart]/I");
  t->Branch("E", E, "E[npart]/D");
  t->Branch("m", M, "m[npart]/D");
  t->Branch("MX1", &mx1, "MX1/D");
  t->Branch("MX2", &mx2, "MX2/D");
  t->Branch("ndim", &ndim, "ndim/I");
  t->Branch("generation_time", &time_gen, "generation_time/F");
  t->Branch("total_time", &time_tot, "total_time/F");

  xsect = result_.s1;
  errxsect = result_.s2;
  ndim = event_.ndim;

  i = 0;
  do {
    tmr.reset();
    generate_(ev);
    if (event_.accepted) i++;
    if (i%5000==0) std::cout << "Event " << i << " generated!" << std::endl;
    time_gen = tmr.elapsed();
    fragmentation_();
    time_tot = tmr.elapsed();
    npart = 0;
    mx1 = remnts_.mx1;
    mx2 = remnts_.mx2;
    for (int p=0; p<lujets_.n; p++) {
      if (p>9) break;
      /*firstdaughterid[npart] = lujets_.k[3][p];
      lastdaughterid[npart] = lujets_.k[4][p];
      if (firstdaughterid!=0 or lastdaughterid!=0) continue;*/
      parentid[npart] = lujets_.k[2][p];

      PID[npart] = lujets_.k[1][p];
      isstable[npart] = lujets_.k[0][p]==1;
      status[npart] = lujets_.k[0][p];
      charge[npart] = luchge_(lujets_.k[1][p])/3.;
      switch (p+1) {
        case 1: { role[npart] = 1; break; }
        case 2: { role[npart] = 2; break; }
        case 3: { role[npart] = 41; break; }
        case 4: { role[npart] = 42; break; }
        case 5: { role[npart] = 3; break; }
        case 6: { role[npart] = 4; break; } // dilepton system
        case 7: { role[npart] = 5; break; }
        case 8: { role[npart] = 6; break; }
        case 9: { role[npart] = 7; break; }
        /*case 10: role[npart] = 3; break;
        case 11: role[npart] = 3; break;
        case 12: role[npart] = 5; break;
        case 13: role[npart] = 5; break;*/
        //default: role[npart] = -1; break; // FIXME need to figure out the origin of each remnant
        default: continue;
      }

      E[npart] = lujets_.p[3][p];
      M[npart] = lujets_.p[4][p];
      mom->SetXYZM(lujets_.p[0][p], lujets_.p[1][p], lujets_.p[2][p], lujets_.p[4][p]);
      pt[npart] = mom->Pt();

      if (pt[npart]!=0.) { eta[npart] = mom->PseudoRapidity(); }
      else eta[npart] = (mom->Pz()/fabs(mom->Pz()))*9999.;

      rapidity[npart] = mom->Rapidity();
      phi[npart] = mom->Phi();

      npart++;
    }
    t->Fill();
  } while (i<nevent);

  t->SaveAs(filename.c_str());

  delete t;
  delete mom;

  return 0;
}
