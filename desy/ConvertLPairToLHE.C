#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
using namespace std;

void ConvertLPairToLHE()
{
  const Double_t energy = 4000;
  const Int_t N = 200; // max number of particles in per event
  const Int_t max_events = 1e2;

  TFile *f1 = TFile::Open("events.root");
  TTree *t1 = (TTree*) f1->Get("h4444");
  Double_t xsec, errxsec;
  Double_t px[N],py[N],pz[N],en[N],m[N];
  Int_t partid[N], parent[N], daughter1[N], daughter2[N], ip, status[N];
  Double_t iz[N];

  t1->SetBranchAddress("xsect",&xsec);
  t1->SetBranchAddress("errxsect",&errxsec);
  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("E",en);
  t1->SetBranchAddress("icode",partid);
  t1->SetBranchAddress("m",m);
  t1->SetBranchAddress("status",status);
  t1->SetBranchAddress("parent",parent);
  t1->SetBranchAddress("daughter1",daughter1);
  t1->SetBranchAddress("daughter2",daughter2);
  //t1->SetBranchAddress("iz",iz);
  t1->SetBranchAddress("ip",&ip);

  ofstream output("gammagammamumu.lpair_inelel_pt15_8tev.lhe");
  //ofstream output("gammagammatautau.lpair_inelel_pt15_7tev.lhe");
  //ofstream output("gammagammatautau.lpair_inelel_tautau_pt25_8tev.lhe");
  //ofstream output("gammagammatautau.lpair_elel_tautau_pt40_7tev.lhe");

  Int_t nevts = t1->GetEntries();
  if(nevts<1) { std::cout << "no event in the file\n"; return;}

  int first_event = 0;

  output << "<LesHouchesEvents version=\"1.0\">"  << endl; 
  output << "<header>" << endl; 
  output << "This file was created from the output of the LPAIR generator" << endl; 
  output << "</header>" << endl; 

  t1->GetEntry(0);
  cout << "xsect = " << xsec << " +/- " << errxsec << endl;

  output << "<init>" << endl;
  output << "2212 2212 " << energy << " " << energy << " 0 0 10042 10042 2 1" << endl;
  output << xsec << " " << errxsec << " 0.26731120000E-03 0" << endl;
  output << "</init>" << endl;

  for(Int_t i = first_event;i < first_event+max_events;i++) {
    t1->GetEntry(i);
    if (i%10000==0)
      cout << i << ", Npart = " << ip << endl;
    
    output << "<event>" << endl;
    output << ip-3 << " 0 0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;
    //	cout << "there are " << ip << " particles in this event\n";
    
    for(int j=0; j<ip; j++) {
      //	Stupid trick to produce inelastic events in both directions!
      
      iz[j] = 0.;
      parent[j] -= 3;
      parent[j] = TMath::Max(parent[j], 0);
      
      if(i%2 == 0)
	pz[j] = -pz[j];
      
      if (partid[j]==2212 && status[j]==21) {
	/*if (fabs(pz[j])==energy) status[j] = -9; // incoming proton
	else*/
	continue;
      }
      else if (partid[j]==22 && status[j]==21) {
	if (j%2==0) iz[j] = 1.;
	else iz[j] = -1.; // helicity
	parent[j] = 0;
	status[j] = -1; // incoming photon
      }
      else if (partid[j]==92) status[j] = 3; // string
      else if ((partid[j]==1 || partid[j]==2 || partid[j]==2101 || partid[j]==2103 || partid[j]==2203) && (status[j]>=11 && status[j]<=13)) status[j] = 3; // quarks content
      else if (status[j]==11) status[j] = 2; // intermediate resonance
      
      //output << partid[j] << " 1 1 2 0 0 " << px[j] << " " << py[j] << " " << pz[j] << " " << en[j] << " " << m[j] << " 0. " << iz[j] << endl;        
      output << fixed << setw(6) << partid[j] << setw(4)
             << status[j] << setw(4)
             << parent[j] << "   0   0   0" << setw(14)
             << px[j] << setw(14)
             << py[j] << setw(14)
             << pz[j] << setw(14)
             << en[j] << setw(10)
             << m[j] << " 0. " << setw(9)
             << iz[j] << endl;
      //        output << "P "<< i*N+j << " " << partid[j] << " " << px[j] << " " << py[j] << " " << pz[j] << " " << en[j] << " 1 0 0 0 0" << endl;      
    }
 
    
    output << "</event>" << endl;
  }
  output << "</LesHouchesEvents>" << endl;
  output.close();
  
  cout << "Converted " << max_events << " events" << endl;
}
