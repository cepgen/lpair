#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;

void ConvertLPairToLHE()
{
  //  TFile *f1 = TFile::Open("lpair-mumu-7tev-pt6pt5.root");
  //TFile *f1 = TFile::Open("lpair-el-tautau-7tev-pt19.root");  
  //TFile *f1 = TFile::Open("lpair-inelel-mumu-pt15-0908.root");
  //TFile *f1 = TFile::Open("lpair-inelel-tautau-pt15-0908.root");
  TFile *f1 = TFile::Open("lpair.root");
  TTree *t1 = (TTree*) f1->Get("h4444");
  const int N = 100; // max number of particles in per event
  Float_t px[N],py[N],pz[N],en[N],m[N];
  Int_t partid[N], ip;
  Float_t iz[N];

  t1->SetBranchAddress("px",px);
  t1->SetBranchAddress("py",py);
  t1->SetBranchAddress("pz",pz);
  t1->SetBranchAddress("E",en);
  t1->SetBranchAddress("icode",partid);
  t1->SetBranchAddress("m",m);
  t1->SetBranchAddress("iz",iz);
  t1->SetBranchAddress("ip",&ip);

  //ofstream output("gammagammamumu.lpair_inelel_pt15_7tev.lhe");
  //ofstream output("gammagammatautau.lpair_inelel_pt15_7tev.lhe");
  //ofstream output("gammagammatautau.lpair_inelel_tautau_pt40_7tev.lhe");
  ofstream output("gammagammatautau.lpair_elel_tautau_pt40_7tev.lhe");

  Int_t nevts = t1->GetEntries();
  if(nevts<1) { std::cout << "no event in the file\n"; return;}

  int first_event = 0;

  output << "<LesHouchesEvents version=\"1.0\">"  << endl; 
  output << "<header>" << endl; 
  output << "This file was created from the output of the LPAIR generator" << endl; 
  output << "</header>" << endl; 

  output << "<init>" << endl;
  output << "2212  2212  0.35000000000E+04  0.35000000000E+04 0 0 10042 10042 2  1" << endl;
  output << "0.10508723460E+01  0.96530000000E-02  0.26731120000E-03   0" << endl;
  output << "</init>" << endl;

  for(Int_t i = first_event;i < first_event+100000;i++)
    {
      t1->GetEntry(i);
      cout << i << ", Npart = " << ip << endl;

      // JH - for filtering on e/mu/had tau decays
      int founde = 0;
      int foundmu = 0;
      int foundpik = 0;
      for(int k = 0;k < ip;k++)
	{
	  if(fabs(partid[k]) == 13)
	    foundmu++;
	  if(fabs(partid[k]) == 11)
	    founde++;
	  if(fabs(partid[k]) == 211 || fabs(partid[k]) == 321)
	    foundpik++;
	}
      
      /* For selecting e+mu tautau events */
      /*
	if(founde == 0 || foundmu == 0)
	continue;
      */

      /* For selecting lepton + 1-prong tautau events - probably only works for ElEl */ 
      /*
      if((founde > 0) && (foundpik != 1))
      	continue;
      if((foundmu > 0) && (foundpik != 1))
      	continue;
      if((foundmu == 0) && (founde == 0))
      	continue;
      if((foundmu > 0) && (founde > 0))
	continue;
      */

      output << "<event>" << endl;
      output << ip+2 << "   0  0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;
      //	cout << "there are " << ip << " particles in this event\n";

      // JH - note here we add in two fake photons as the beam particles. The energies don't matter - this is only for 
      // the LHE event record.       
      output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+02  0.10000000000E+02  0.00000000000E+00 0.  1." << endl;
      output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+00  0.10000000000E+02  0.00000000000E+00 0. -1." << endl;   
  
      for(int j=0; j<ip; j++) {
	//	Stupid trick to produce inelastic events in both directions!
	if(i%2 == 0)
	  pz[j] = -1.0*pz[j];
	
	output << partid[j] << " 1 1 2 0 0 " << px[j] << " " << py[j] << " " << pz[j] << " " << en[j] << " " << m[j] << " 0. " << iz[j] << endl; 	
	//	output << "P "<< i*N+j << " " << partid[j] << " " << px[j] << " " << py[j] << " " << pz[j] << " " << en[j] << " 1 0 0 0 0" << endl;      
      }
      

      output << "</event>" << endl;
    }
  output << "</LesHouchesEvents>" << endl;
  output.close();

  cout << "Converted " << i << " events" << endl;
}
