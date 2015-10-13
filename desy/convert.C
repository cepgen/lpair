#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

void convert()
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
  stringstream ss;
  Int_t np;

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
    if (i%10000==0) cout << i << ", Npart = " << ip << endl;
    

    /*
      Event content :
      (0) : incoming proton 1 (elastic, or dissociating)
      (1) : incoming proton 2 (elastic only)
      (2) : photon from proton 1
      (3) : photon from proton 2
      (4) : outgoing proton (or remnant) 1
      (5) : muon 1
      (6) : dimuon system
      (7) : muon 2
      (8) : outgoing proton 2
      (9) : quark (from remnants)
      (10) : diquark (from remnants)
     */
    cout << "there are " << ip << " particles (" << ip-11 << " from remnants) in this event\n";

    ss.str("");

    for(int j=0; j<ip; j++) {
      // Stupid trick to produce inelastic events in both directions!
      if(i%2 == 0) pz[j] = -pz[j];
    }
    
    ss << partid[2] << " " << status[2] << " -1 0 0 0 0 " << px[2] << " " << py[2] << " " << pz[2] << " " << en[2] << " " << m[2] << " 0 0" << endl; np++;
    ss << partid[3] << " " << status[3] << " -1 0 0 0 0 " << px[3] << " " << py[3] << " " << pz[3] << " " << en[3] << " " << m[3] << " 0 0" << endl; np++;
    if (status[4]==1) { // elastic case
      ss << partid[4] << " " << status[4] << " 0 0 0 0 " << px[4] << " " << py[4] << " " << pz[4] << " " << en[4] << " " << m[4] << " 0 0" << endl; np++;
    }
    ss << partid[8] << " " << status[8] << " 0 0 0 0 " << px[8] << " " << py[8] << " " << pz[8] << " " << en[8] << " " << m[8] << " 0 0" << endl; np++;
    // Outgoing muons
    ss << partid[5] << " " << status[5] << " 1 0 0 0 " << px[5] << " " << py[5] << " " << pz[5] << " " << en[5] << " " << m[5] << " 0 0" << endl; np++;
    ss << partid[7] << " " << status[7] << " 2 0 0 0 " << px[7] << " " << py[7] << " " << pz[7] << " " << en[7] << " " << m[7] << " 0 0" << endl; np++;

    if (status[4]==21) { // single-dissociative case
      for (Int_t j=11; j<ip; j++) {
	if (status[j]!=1) continue;
	ss << partid[j] << " " << status[j] << " 1 2 0 0 " << px[j] << " " << py[j] << " " << pz[j] << " " << en[j] << " " << m[j] << " 0 0" << endl; np++;
      }
    }

    output << "<event>" << endl;
    output << np << " 0 0.298346E-04  0.91188E+02  0.7546772E-02  0.13E+00" << endl;
    output << ss.str();
    output << "</event>" << endl;
    
    /*
      
      if (partid[j]==2212 && status[j]==21) {
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
      output << partid[j] << " " 
	     << status[j] << " " 
	     << parent[j] << " 0 0 0 " 
	     << px[j] << " " 
	     << py[j] << " " 
	     << pz[j] << " " 
	     << en[j] << " " 
	     << m[j] << " 0. " 
	     << iz[j] << endl; 	
      //	output << "P "<< i*N+j << " " << partid[j] << " " << px[j] << " " << py[j] << " " << pz[j] << " " << en[j] << " 1 0 0 0 0" << endl;      
    }*/
    
  }
  output << "</LesHouchesEvents>" << endl;
  output.close();
  
  cout << "Converted " << max_events << " events" << endl;
}
