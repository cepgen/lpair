#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "includes/OpenTree.C"
#include "includes/TDRStyle.C"
#include "includes/GlobalFunctions.C"

#define PI 3.14159265359
#define MUON_MASS 0.1057

#define NUM_CUTS 13
#define CUT_INVM_MIN 20.
//#define CUT_INVM_MIN 160. // FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUT_PT_TRAILING 20.
#define CUT_PT_LEADING 20.

#define CUT_ETA_MAX 2.4

#define CUT_VTX_CHISQ_MIN 1E-3
#define CUT_VTX_ZSEPARATION_MIN 0.2
#define CUT_VTX_Z_MAX 24.0

#define CUT_DPT_MAX 1.0
#define CUT_DPHI_MIN 0.9
#define CUT_OPENANGLE_MAX 0.96
#define CUT_ACOP_MAX 0.1

#define CUT_VTX_EXCLU_TRACKCOUNTING_MAX 0.1

#define CUT_ZDC_EM_MAX 9999918
#define CUT_ZDC_HAD_MAX 9999120

#define INV_MASS_SPLIT 20
#define INV_MASS_ZPEAK_LOW 70.
#define INV_MASS_ZPEAK_UP 106.

#define ZPEAK 0
#define HIGHM 1
#define FULLM 2

using namespace std;

void Processing() {

  TObjArray *hList = new TObjArray(0);
  TString mass_range[] = {"Zpeak", "highm", "full"};
  const Int_t masses = sizeof(mass_range)/sizeof(TString);

  stringstream ss;
  Double_t cut_pTPair_min, cut_pTPair_max, cut_acop, cut_dpt;
  TString cut_nExtTrk;
  TString output_file;
  Bool_t drawAfterwards, antisel;
  Int_t pair_p, pair_m, trailing, leading, num_pvc;
  Double_t etap, etam, ptp, ptm, etaLead, etaTrail, ptLead, ptTrail;
  Double_t mup_eff, mum_eff, mu17_eff, mu8_eff, mu_eff;
  Double_t vtx_t;
  Double_t pu_weight, weight, weight_nopu;
  Double_t et_sumPt, rap_pair, acop, num_tracks;
  TLorentzVector Pm, Pp;

  TH1D *extTrk[masses];
  TH1D *invm[masses];
  TH1D *invm_over160GeV[masses];
  TH1D *ptpair[masses], *ptpairZoom1[masses], *ptpairZoom2[masses], *ptpairZoom3[masses], *ptpairZoom4[masses];
  TH1D *pTSingleM[masses], *pTSingleP[masses];
  TH1D *etaPair[masses], *etaSingleM[masses], *etaSingleP[masses];
  TH1D *dpt[masses], *dptZoom1[masses];
  TH1D *numVtx[masses], *numVtxAfterCuts[masses], *numVtxNoPURW[masses];
  TH1D *acoplZoom1[masses], *acoplZoom2[masses], *acoplZoom3[masses], *acoplZoom4[masses];
  TH1D *etPt[masses], *etEta[masses], *etPhi[masses], *etSumPt[masses], *etDvtx[masses];
  TH1D *etPtAfter[masses], *etEtaAfter[masses];
  TH2D *acoplVsdpt[masses], *ntrkVsnvtx[masses], *acoplVsptpair[masses];
  TH2D *nvtxVspuw[masses];
  TH1D *vtxT[masses], *vtxZ[masses];
  TH1D *hCuts, *hCuts_numAfter, *hTrig;

  gSystem->Load("includes/OpenTree_C.so");
  gSystem->Load("includes/TDRStyle_C.so");
  gSystem->Load("includes/GlobalFunctions_C.so");

  //output_file = "1022-lumiSelection-0exttrk.root";
  //output_file = "1022-lumiSelection-1exttrk.root";
  //output_file = "1022-lumiSelection-1-6exttrk.root";
  //output_file = "1022-antilumiSelection-0exttrk.root";
  //output_file = "1022-antilumiSelection-0exttrk-m160.root";
  //output_file = "1022-antilumiSelection-1-6exttrk.root";
  //output_file = "1022-noSelection-0exttrk.root";
  //output_file = "1022-noSelection-0exttrk-m_gt_160.root";
  //output_file = "1022-noSelection-allexttrk.root";
  //output_file = "1022-noSelection-1-6exttrk.root";
  //output_file = "1022-signalSelection-0exttrk.root";
  //output_file = "1022-signalSelection-1-6exttrk.root";
  output_file = "test.root";

  ///////////////////////////////////////////
  //        Parameters for plotting        //
  cut_pTPair_min = 0.0;                    // in GeV/c
  cut_pTPair_max = 0.0;                    // in GeV/c
  cut_acop       = 0.1;                    // between 0.0 and 1.0 (0. = no cut)
  cut_dpt        = 1.0;                    // in GeV/c
  antisel        = false;                  // anti-selection for (acop,dpt) ?
//antisel        = true;                   //
  cut_nExtTrk    = "0";                    // 0, 1, 1-6, 0-6, 0-10, all (no cut)
  ///////////////////////////////////////////

  drawAfterwards = false;

  ofstream eventslist("eventsList.txt");

  // Adding paths to files (usage: AddFile((TString)name, (TString)path, (bool)is_monte_carlo, (Double_t)normalization_factor))

  //AddFile("signal", "CalcHEP-Signal-5kEvents.root", true, 1.);
  //AddFile("qcd", "QCD-MuEnrichedPt10-Max15ExtTrk-2806.root", true , 1.);

  AddFile("dymumu", "DYtoMuMu-powheg-1122.root", true);
  AddFile("dytautau", "DYtoTauTau-1122.root", true);
  AddFile("data", "Data-1122.root", false);
  AddFile("tt", "TT-uncomplete-1122.root", true);
  AddFile("ww-nlo", "WW-CT10-mcatnlo-1122.root", true);
  AddFile("ww-z2", "WW-TuneZ2-pythia6-1122.root", true);
  AddFile("elel", "ElEl-15GeV-MuMu-Max15ExtTrk-1121.root", true, 1.);
  AddFile("inelel", "InelEl-15GeV-MuMu-Max15ExtTrk-1121.root", true, 1.);
  AddFile("inelinel", "InelInel-15GeV-MuMu-Max15ExtTrk-1121.root", true, 1.);

  const Int_t numFiles = filename.size();

  GetMuonEfficiencyCorrections();
  GetTriggerEfficiencyCorrections();

  for (Int_t i=0; i<numFiles; ++i) {
    
    cerr << "=== " << filename.at(i) << " (" << numEvents.at(i) << " events) ===" << endl;
    cout << "================================================================================" << endl;
    cout << "Opening file \n"
         << "\t" << filename.at(i) << "\n"
         << "\t" << setprecision(7) << numEvents.at(i) << " events\n"
         << "\tMC ? " << isMC.at(i)
         << endl;
    cout << "================================================================================" << endl;
    TFile *file = new TFile(filename.at(i));
    TTree *tree = (TTree*)(file->Get("ntp1"));

    OpenTree(tree, i, isMC.at(i));

    pair_p = -1;
    pair_m = -1;
    trailing = -1;
    leading = -1;
    et_sumPt = 0.;

    // Building histograms
    for (Int_t j=0; j<masses; j++) {
      ss.str(""); ss << "numExtraTracks_" << mass_range[j] << "_" << datasetName.at(i);
      extTrk[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),6,-0.5,5.5); hList->Add(extTrk[j]);
      ss.str(""); ss << "invariantMass_ov160_" << mass_range[j] << "_" << datasetName.at(i);
      invm_over160GeV[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),70,150.,500.); hList->Add(invm_over160GeV[j]);
      ss.str(""); ss << "invariantMass_" << mass_range[j] << "_" << datasetName.at(i);
      invm[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),53,35.,300.); hList->Add(invm[j]);
      ss.str(""); ss << "pTpair_" << mass_range[j] << "_" << datasetName.at(i);
      //ptpair[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),135,0.,10.125); hList->Add(ptpair[j]);
      //ptpair[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),100,0.,50.); hList->Add(ptpair[j]);
      //ptpair[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),50,0.,20.); hList->Add(ptpair[j]);
      ptpair[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),120,0.,300.); hList->Add(ptpair[j]);
      ss.str(""); ss << "pTpairZoom1_" << mass_range[j] << "_" << datasetName.at(i);
      ptpairZoom1[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),25,0.,50.); hList->Add(ptpairZoom1[j]);
      ss.str(""); ss << "pTpairZoom2_" << mass_range[j] << "_" << datasetName.at(i);
      ptpairZoom2[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),50,0.,50.); hList->Add(ptpairZoom2[j]);
      ss.str(""); ss << "pTpairZoom3_" << mass_range[j] << "_" << datasetName.at(i);
      ptpairZoom3[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),20,0.,5.); hList->Add(ptpairZoom3[j]);
      ss.str(""); ss << "pTpairZoom4_" << mass_range[j] << "_" << datasetName.at(i);
      ptpairZoom4[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),48,0.,120.); hList->Add(ptpairZoom4[j]);
      ss.str(""); ss << "numVtx_" << mass_range[j] << "_" << datasetName.at(i);
      numVtx[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),20,-0.5,19.5); hList->Add(numVtx[j]);
      ss.str(""); ss << "numVtx_noPURW_" << mass_range[j] << "_" << datasetName.at(i);
      numVtxNoPURW[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),20,-0.5,19.5); hList->Add(numVtxNoPURW[j]);
      ss.str(""); ss << "numVtxAfterCuts_" << mass_range[j] << "_" << datasetName.at(i);
      numVtxAfterCuts[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),20,-0.5,19.5); hList->Add(numVtxAfterCuts[j]);
      ss.str(""); ss << "dpt_" << mass_range[j] << "_" << datasetName.at(i);
      dpt[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),20,0.,1.); hList->Add(dpt[j]);
      ss.str(""); ss << "dptZoom1_" << mass_range[j] << "_" << datasetName.at(i);
      dptZoom1[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),10,0.,1.); hList->Add(dptZoom1[j]);
      ss.str(""); ss << "etaSingleM_" << mass_range[j] << "_" << datasetName.at(i);
      etaSingleM[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),16,-2.4,2.4); hList->Add(etaSingleM[j]);
      ss.str(""); ss << "etaSingleP_" << mass_range[j] << "_" << datasetName.at(i);
      etaSingleP[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),16,-2.4,2.4); hList->Add(etaSingleP[j]);
      ss.str(""); ss << "pTSingleM_" << mass_range[j] << "_" << datasetName.at(i);
      pTSingleM[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),75,20.,95.); hList->Add(pTSingleM[j]);
      ss.str(""); ss << "pTSingleP_" << mass_range[j] << "_" << datasetName.at(i);
      pTSingleP[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),75,20.,95.); hList->Add(pTSingleP[j]);
      ss.str(""); ss << "etaPair_" << mass_range[j] << "_" << datasetName.at(i);
      etaPair[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),40,-4.0,4.0); hList->Add(etaPair[j]);
      ss.str(""); ss << "acoplVsdpt_" << mass_range[j] << "_" << datasetName.at(i);
      acoplVsdpt[j] = new TH2D(ss.str().c_str(),ss.str().c_str(),100,0.,1.,100,0.,5.); hList->Add(acoplVsdpt[j]);
      ss.str(""); ss << "acoplVsptpair_" << mass_range[j] << "_" << datasetName.at(i);
      acoplVsptpair[j] = new TH2D(ss.str().c_str(),ss.str().c_str(),100,0.,1.,15,0.,30.); hList->Add(acoplVsptpair[j]);
      ss.str(""); ss << "ntrkVsnvtx_" << mass_range[j] << "_" << datasetName.at(i);
      ntrkVsnvtx[j] = new TH2D(ss.str().c_str(),ss.str().c_str(),6,-0.5,5.5,15,-0.5,14.5); hList->Add(ntrkVsnvtx[j]);
      ss.str(""); ss << "nvtxVspuw_" << mass_range[j] << "_" << datasetName.at(i);
      nvtxVspuw[j] = new TH2D(ss.str().c_str(),ss.str().c_str(),15,-0.5,14.5,50,0.,2.5); hList->Add(nvtxVspuw[j]);
      ss.str(""); ss << "vtxZ_" << mass_range[j] << "_" << datasetName.at(i);
      vtxZ[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),15,-30,30); hList->Add(vtxZ[j]);
      ss.str(""); ss << "vtxT_" << mass_range[j] << "_" << datasetName.at(i);
      vtxT[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),500,0.,0.6); hList->Add(vtxT[j]);
      ss.str(""); ss << "acoplZoom1_" << mass_range[j] << "_" << datasetName.at(i);
      acoplZoom1[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),10,0.,1.); hList->Add(acoplZoom1[j]);
      ss.str(""); ss << "acoplZoom2_" << mass_range[j] << "_" << datasetName.at(i);
      acoplZoom2[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),100,0.,1.); hList->Add(acoplZoom2[j]);
      ss.str(""); ss << "acoplZoom3_" << mass_range[j] << "_" << datasetName.at(i);
      acoplZoom3[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),100,0.,0.1); hList->Add(acoplZoom3[j]);
      ss.str(""); ss << "acoplZoom4_" << mass_range[j] << "_" << datasetName.at(i);
      acoplZoom4[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),50,0.,0.1); hList->Add(acoplZoom4[j]);
      ss.str(""); ss << "etPt_" << mass_range[j] << "_" << datasetName.at(i);
      etPt[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),50,0.,20.); hList->Add(etPt[j]);
      ss.str(""); ss << "etEta_" << mass_range[j] << "_" << datasetName.at(i);
      etEta[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),60,-3.,3.); hList->Add(etEta[j]);
      ss.str(""); ss << "etPhi_" << mass_range[j] << "_" << datasetName.at(i);
      etPhi[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),16,-4.,4.); hList->Add(etPhi[j]);
      ss.str(""); ss << "etSumPt_" << mass_range[j] << "_" << datasetName.at(i);
      etSumPt[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),100,0.,200.); hList->Add(etSumPt[j]);
      ss.str(""); ss << "etDvtx_" << mass_range[j] << "_" << datasetName.at(i);
      etDvtx[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),500,0.,5.); hList->Add(etDvtx[j]);
      ss.str(""); ss << "etPtAfter_" << mass_range[j] << "_" << datasetName.at(i);
      etPtAfter[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),50,0.,20.); hList->Add(etPtAfter[j]);
      ss.str(""); ss << "etEtaAfter_" << mass_range[j] << "_" << datasetName.at(i);
      etEtaAfter[j] = new TH1D(ss.str().c_str(),ss.str().c_str(),60,-3.,3.); hList->Add(etEtaAfter[j]);
    }
    ss.str(""); ss << "triggersMatched_" << datasetName.at(i);
    hTrig = new TH1D(ss.str().c_str(),ss.str().c_str(),4,-0.5,3.5); hList->Add(hTrig);
    ss.str(""); ss << "Cuts_" << datasetName.at(i);
    hCuts = new TH1D(ss.str().c_str(),ss.str().c_str(),NUM_CUTS,0.,NUM_CUTS); hList->Add(hCuts);
    ss.str(""); ss << "Cuts_NumAfter_" << datasetName.at(i);
    hCuts_numAfter = new TH1D(ss.str().c_str(),ss.str().c_str(),NUM_CUTS,0.,NUM_CUTS); hList->Add(hCuts_numAfter);


    for (Double_t evt=0; evt<numEvents.at(i); ++evt) {

      if (fmod(evt,5000)==0) {
	cerr << "[";
	cerr << setprecision(4) << setw(5)
	     << evt/numEvents.at(i)*100 << "%";
	cerr << "] Processing event " 
	     << setprecision(7) << evt << " / " << numEvents.at(i) << endl;
      }

      tree->GetEntry(evt);

      if(var_charge[i][var_Pair[i][0]]>0){
	pair_p = var_Pair[i][0]; pair_m = var_Pair[i][1];
      }
      else { pair_p = var_Pair[i][1]; pair_m = var_Pair[i][0]; }
      if (var_pt[i][var_Pair[i][0]]>var_pt[i][var_Pair[i][1]]) {
	leading = var_Pair[i][0]; trailing = var_Pair[i][1];
      }
      else { leading = var_Pair[i][1]; trailing = var_Pair[i][0]; }

      //////////////////////////////////////////////////////////////////////////
      //                        Efficiency computation                        //
      //////////////////////////////////////////////////////////////////////////
       
      pu_weight = 1.;
      weight = 1.;

      if (isMC.at(i)) { // Muon efficiency correction only for MC
	mup_eff = 1.;
	mum_eff = 1.;
	mu17_eff = 1.;
	mu8_eff = 1.;
	mu_eff = 1.;

	weight_nopu = 1.;
      
	etap = fabs(var_eta[i][pair_p]);
	etam = fabs(var_eta[i][pair_m]);
	ptp = var_pt[i][pair_p];
	ptm = var_pt[i][pair_m];
	etaLead = fabs(var_eta[i][leading]);
	etaTrail = fabs(var_eta[i][trailing]);
	ptLead = var_pt[i][leading];
	ptTrail = var_pt[i][trailing];
	if (evt<0.453435115*numEvents.at(i)) {
	  // Muon efficiency for Run2011A
	  for (Int_t l=0; l<(Int_t)(etaLowA.size()); l++) {
	    if (etap>etaLowA.at(l) && etap<etaHighA.at(l) &&
		ptp>pTLowA.at(l) && ptp<pTHighA.at(l))
	      mup_eff = MuEffA.at(l);
	    if (etam>etaLowA.at(l) && etam<etaHighA.at(l) &&
		ptm>pTLowA.at(l) && ptm<pTHighA.at(l))
	      mum_eff = MuEffA.at(l);
	  }
	}
	else {
	  // Muon efficiency for Run2011B
	  for (Int_t l=0; l<(Int_t)(etaLowB.size()); l++) {
	    if (etap>etaLowB.at(l) && etap<etaHighB.at(l) &&
		ptp>pTLowB.at(l) && ptp<pTHighB.at(l))
	      mup_eff = MuEffB.at(l);
	    if (etam>etaLowB.at(l) && etam<etaHighB.at(l) &&
		ptm>pTLowB.at(l) && ptm<pTHighB.at(l))
	      mum_eff = MuEffB.at(l);
	  }
	}
	for (Int_t l=0; l<(Int_t)(TEetaLow17.size()); l++) {
	  if (etaLead>TEetaLow17.at(l) && etaLead<TEetaHigh17.at(l) &&
	      ptLead>TEpTLow17.at(l) && ptLead<TEpTHigh17.at(l))
	    mu17_eff = TEff17.at(l);
	}
	for (Int_t l=0; l<(Int_t)(TEetaLow8.size()); l++) {
	  if (etaTrail>TEetaLow8.at(l) && etaTrail<TEetaHigh8.at(l) &&
	      ptTrail>TEpTLow8.at(l) && ptTrail<TEpTHigh8.at(l))
	    mu8_eff = TEff8.at(l);
	}
	pu_weight = var_PUweight[i][0];
	mu_eff = mup_eff*mum_eff;
	mu_eff = mu_eff*mu17_eff*mu8_eff;
	
	weight_nopu = weight*mu_eff;
	weight = weight_nopu*pu_weight;
	//cout << weight_nopu << "\t" << weight << "\t" << pu_weight << endl;
      }
      else { // Data
	weight = 1.;
	weight_nopu = 1.;
      }

      //////////////////////////////////////////////////////////////////////////
      //                              HLT Filter                              //
      //////////////////////////////////////////////////////////////////////////

      Int_t numTriggersMatched(0);
      if (PassesTrigger(hlt_d[4][i][0],var_run[i][0]))
	numTriggersMatched += 1;
      if (PassesTrigger(hlt_d[5][i][0],var_run[i][0]))
	numTriggersMatched += 2;
      if (PassesTrigger(hlt_d[6][i][0],var_run[i][0]))
	numTriggersMatched += 4;

      if (numTriggersMatched==0) { // No trigger was matched
	hCuts->Fill(0., weight);
	continue;
      }
      hCuts_numAfter->Fill(0.,weight);


      ///////////////////////////////////////////////////////////////////////////
      //               Loop on all the primary vertex candidates               //
      ///////////////////////////////////////////////////////////////////////////
      
      num_tracks = -999;
      Int_t vtx = -999;
      
      Int_t num_vtx_with_dimuon(0);
      for (Int_t primvtx=0; primvtx<var_nvtx[i]; primvtx++) {
	if (var_vtxhasmumu[i][primvtx] != 1) continue;
	num_vtx_with_dimuon++;
	vtx = primvtx;
      }
      if (num_vtx_with_dimuon>1) {
	cout << "*** ERROR *** Ambiguosity in event " << evt 
	     << ":\n   multiple (" << num_vtx_with_dimuon 
	     << ") vertices with two muons arising" << endl;
	hCuts->Fill(1., weight);
	break;
      }
      if(fabs(var_vtxZ[i][vtx])>=24) continue;
      hCuts_numAfter->Fill(1., weight);

      num_tracks = (var_vtxTrack_is_int[i]) ? 
	var_vtxTrack_int[i][vtx] : var_vtxTrack[i][vtx];
      
      if (num_tracks-2<0) continue;

      //////////////////////////////////////////////////////////////////////////
      //                           Kinematics cuts                            //
      //////////////////////////////////////////////////////////////////////////
      
      if (var_mass[i][0]<CUT_INVM_MIN) {
	hCuts->Fill(2., weight);
	continue;
      }
      hCuts_numAfter->Fill(2., weight);

      if (var_pt[i][trailing]<CUT_PT_TRAILING || var_pt[i][leading]<CUT_PT_LEADING) {
	hCuts->Fill(3., weight);
	continue;
      }
      hCuts_numAfter->Fill(3., weight);

      if (var_nhitsTrack[i][pair_p]<8 || var_nhitsTrack[i][pair_m]<8) {
	hCuts->Fill(4., weight);
	continue;
      }
      hCuts_numAfter->Fill(4., weight);
      
      if (!var_MuCand_tightID[i][pair_p] || !var_MuCand_tightID[i][pair_m]) {
	hCuts->Fill(5);
	continue;
      }
      hCuts_numAfter->Fill(5., weight);

      if (fabs(var_eta[i][pair_p])>=2.4 || fabs(var_eta[i][pair_m])>=2.4) {
	hCuts->Fill(6);
	continue;
      }
      hCuts_numAfter->Fill(6., weight);

      if (!var_tracker[i][pair_p] || !var_global[i][pair_p]) {
	hCuts->Fill(7);
	continue;
      }
      if (!var_tracker[i][pair_m] || !var_global[i][pair_m]) {
	hCuts->Fill(7., weight);
	continue;
      }
      hCuts_numAfter->Fill(7., weight);
      
      /*if (!(bool)(isMC.at(i))) {
	cout << "candidate " << setw(6) << filterEvents << "\n"
	<< "\tRun " << var_run[i][0] << "  LS " << var_ls[i][0] << "\tEvt " << var_event[i][0] << "\n"
	<< "\tmass=" << var_mass[i][0] << " GeV" << "\n"
	<< "\tOpening angle=" << openangle << "\n"
	<< "\t\tdphi/pi=" << var_dphi[i][0]/PI <<"\n"
	<< "\t\tdpt=" << var_dpt[i][0]
	<< endl;
	}*/
      Pm.SetPtEtaPhiM(var_pt[i][pair_m], var_eta[i][pair_m], var_phi[i][pair_m], MUON_MASS);
      Pp.SetPtEtaPhiM(var_pt[i][pair_p], var_eta[i][pair_p], var_phi[i][pair_p], MUON_MASS);
      rap_pair = (Pm+Pp).Eta();
      acop = 1-fabs(Pm.DeltaPhi(Pp)/PI);

      if (cut_dpt != 0. && cut_acop != 0.) {
	// We are looking at the "lumi" (or "anti-lumi") selection
	if (!antisel) {
	  if (var_dpt[i][0]>cut_dpt) {
	    hCuts->Fill(8., weight);
	    continue;
	  }
	  hCuts_numAfter->Fill(8., weight);

	  if (acop>cut_acop) {
	    hCuts->Fill(9., weight);
	    continue;
	  }
	  hCuts_numAfter->Fill(9., weight);
	}
	else { // antiselection
	  if (var_dpt[i][0]<=cut_dpt && acop<=cut_acop) {
	    hCuts->Fill(10., weight);
	    continue;
	  }
	  hCuts_numAfter->Fill(10., weight);
	}
      }
      
      // We fill the extra tracks on vertex histogram
      if (((cut_pTPair_min != 0. && var_MuMupt[i][0]>=cut_pTPair_min) || cut_pTPair_min==0.) ||
	  ((cut_pTPair_max != 0. && var_MuMupt[i][0]<=cut_pTPair_max) || cut_pTPair_max==0.)) {

	// Tricky part - to fit Jonathan's cut flow we have to
	// ensure the pT of the pair is compatible with the selection
	// without cutting on it at this step

	num_pvc = var_nvtx[i];

	if (var_mass[i][0]>=INV_MASS_SPLIT) {
	  if (var_mass[i][0]>INV_MASS_ZPEAK_LOW && var_mass[i][0]<INV_MASS_ZPEAK_UP) {
	    extTrk[ZPEAK]->Fill(num_tracks-2-0.5, weight);
	    numVtx[ZPEAK]->Fill(num_pvc-0.5, weight);
	    numVtxNoPURW[ZPEAK]->Fill(num_pvc-0.5, weight_nopu);
	    ntrkVsnvtx[ZPEAK]->Fill(num_tracks-2-0.5, num_pvc-0.5, weight);
	    nvtxVspuw[ZPEAK]->Fill(num_pvc-0.5, pu_weight);
	  }
	  else {
	    extTrk[HIGHM]->Fill(num_tracks-2-0.5, weight);
	    numVtx[HIGHM]->Fill(num_pvc-0.5, weight);
	    numVtxNoPURW[HIGHM]->Fill(num_pvc-0.5, weight_nopu);
	    ntrkVsnvtx[HIGHM]->Fill(num_tracks-2-0.5, num_pvc-0.5, weight);
	    nvtxVspuw[HIGHM]->Fill(num_pvc-0.5, pu_weight);
	  }
	}
	extTrk[FULLM]->Fill(num_tracks-2-0.5, weight);
	numVtx[FULLM]->Fill(num_pvc-0.5, weight);
	numVtxNoPURW[FULLM]->Fill(num_pvc-0.5, weight_nopu);
	ntrkVsnvtx[FULLM]->Fill(num_tracks-2-0.5, num_pvc-0.5, weight);
	nvtxVspuw[FULLM]->Fill(num_pvc-0.5, pu_weight);
      }
      
      // Extra tracks on dimuon vertex plots
      et_sumPt = 0.;
      for (Int_t et=0; et<var_nTrack[i][0]; et++) { // loop on extra tracks on vertex
	if (var_et_vtx_dist[i][et]>0.01 || var_et_quality[i][et] != 1) continue;

	if (var_mass[i][0]>=INV_MASS_SPLIT) {
	  if (var_mass[i][0]>INV_MASS_ZPEAK_LOW && var_mass[i][0]<INV_MASS_ZPEAK_UP) {
	    etDvtx[ZPEAK]->Fill(var_et_vtx_dist[i][et], weight);
	    etPt[ZPEAK]->Fill(var_et_pt[i][et], weight);
	    etEta[ZPEAK]->Fill(var_et_eta[i][et], weight);
	    etPhi[ZPEAK]->Fill(var_et_phi[i][et], weight);
	  }
	  else {
	    etDvtx[HIGHM]->Fill(var_et_vtx_dist[i][et], weight);
	    etPt[HIGHM]->Fill(var_et_pt[i][et], weight);
	    etEta[HIGHM]->Fill(var_et_eta[i][et], weight);
	    etPhi[HIGHM]->Fill(var_et_phi[i][et], weight);
	  }
	}
	etDvtx[FULLM]->Fill(var_et_vtx_dist[i][et], weight);
	etPt[FULLM]->Fill(var_et_pt[i][et], weight);
	etEta[FULLM]->Fill(var_et_eta[i][et], weight);
	etPhi[FULLM]->Fill(var_et_phi[i][et], weight);
	et_sumPt += var_et_pt[i][et];
      }
      // Cut on number of extra tracks on mu mu vertex
      if (cut_nExtTrk != "" and cut_nExtTrk != "all") {
	if (cut_nExtTrk == "0")
	  if (num_tracks-2!=0) {
	    hCuts->Fill(11., weight);
	    continue;
	  }
	if (cut_nExtTrk == "1")
	  if (num_tracks-2!=1) {
	    hCuts->Fill(11., weight);
	    continue;
	  }
	if (cut_nExtTrk == "1-6")
	  if (num_tracks-2<1 || num_tracks-2>6) {
	    hCuts->Fill(11., weight);
	    continue;
	  }
	if (cut_nExtTrk == "0-6")
	  if (num_tracks-2>6) {
	    hCuts->Fill(11., weight);
	    continue;
	  }
	if (cut_nExtTrk == "0-10")
	  if (num_tracks-2>10) {
	    hCuts->Fill(11., weight);
	    continue;
	  }
      }
      hCuts_numAfter->Fill(11., weight);

      if (cut_pTPair_min != 0. && var_MuMupt[i][0]<cut_pTPair_min) {
	hCuts->Fill(12., weight);
	continue;
      }
      if (cut_pTPair_max != 0. && var_MuMupt[i][0]>cut_pTPair_max) {
	hCuts->Fill(12., weight);
	continue;
      }
      hCuts_numAfter->Fill(12., weight);

      vtx_t = sqrt(pow(var_vtxX[i][vtx],2)+pow(var_vtxY[i][vtx],2));

      if (var_mass[i][0]>=INV_MASS_SPLIT) {
	if (var_mass[i][0]>INV_MASS_ZPEAK_LOW && var_mass[i][0]<INV_MASS_ZPEAK_UP) {
	  invm[ZPEAK]->Fill(var_mass[i][0], weight);
	  invm_over160GeV[ZPEAK]->Fill(var_mass[i][0], weight);
	  ptpair[ZPEAK]->Fill(var_MuMupt[i][0], weight);
	  ptpairZoom1[ZPEAK]->Fill(var_MuMupt[i][0], weight);
	  ptpairZoom2[ZPEAK]->Fill(var_MuMupt[i][0], weight);
	  ptpairZoom3[ZPEAK]->Fill(var_MuMupt[i][0], weight);
	  pTSingleM[ZPEAK]->Fill(var_pt[i][pair_m], weight);
	  pTSingleP[ZPEAK]->Fill(var_pt[i][pair_p], weight);
	  etaSingleM[ZPEAK]->Fill(var_eta[i][pair_m], weight);
	  etaSingleP[ZPEAK]->Fill(var_eta[i][pair_p], weight);
	  etaPair[ZPEAK]->Fill(rap_pair, weight);
	  acoplZoom1[ZPEAK]->Fill(acop, weight);
	  acoplZoom2[ZPEAK]->Fill(acop, weight);
	  acoplZoom3[ZPEAK]->Fill(acop, weight);
	  acoplZoom4[ZPEAK]->Fill(acop, weight);
	  acoplVsdpt[ZPEAK]->Fill(acop, var_dpt[i][0], weight);
	  acoplVsptpair[ZPEAK]->Fill(acop, var_MuMupt[i][0], weight);
	  vtxZ[ZPEAK]->Fill(var_vtxZ[i][vtx], weight);
	  vtxT[ZPEAK]->Fill(vtx_t, weight);
	  dpt[ZPEAK]->Fill(var_dpt[i][0], weight);
	  dptZoom1[ZPEAK]->Fill(var_dpt[i][0], weight);
	  etSumPt[ZPEAK]->Fill(et_sumPt, weight);
	  for (Int_t et=0; et<var_nTrack[i][0]; et++) {
	    // loop on extra tracks on the dimuon vertex
	    if (var_et_vtx_dist[i][et]>0.01 ||
		var_et_quality[i][et] != 1)
	      continue;
	    etPtAfter[ZPEAK]->Fill(var_et_pt[i][et], weight);
	    etEtaAfter[ZPEAK]->Fill(var_et_eta[i][et], weight);
	  }
	  numVtxAfterCuts[ZPEAK]->Fill(num_pvc-0.5, weight);
	}
	else {
	  invm[HIGHM]->Fill(var_mass[i][0], weight);
	  invm_over160GeV[HIGHM]->Fill(var_mass[i][0], weight);
	  ptpair[HIGHM]->Fill(var_MuMupt[i][0], weight);
	  ptpairZoom1[HIGHM]->Fill(var_MuMupt[i][0], weight);
	  ptpairZoom2[HIGHM]->Fill(var_MuMupt[i][0], weight);
	  ptpairZoom3[HIGHM]->Fill(var_MuMupt[i][0], weight);
	  pTSingleM[HIGHM]->Fill(var_pt[i][pair_m], weight);
	  pTSingleP[HIGHM]->Fill(var_pt[i][pair_p], weight);
	  etaSingleM[HIGHM]->Fill(var_eta[i][pair_m], weight);
	  etaSingleP[HIGHM]->Fill(var_eta[i][pair_p], weight);
	  etaPair[HIGHM]->Fill(rap_pair, weight);
	  acoplZoom1[HIGHM]->Fill(acop, weight);
	  acoplZoom2[HIGHM]->Fill(acop, weight);
	  acoplZoom3[HIGHM]->Fill(acop, weight);
	  acoplZoom4[HIGHM]->Fill(acop, weight);
	  acoplVsdpt[HIGHM]->Fill(acop, var_dpt[i][0], weight);
	  acoplVsptpair[HIGHM]->Fill(acop, var_MuMupt[i][0], weight);
	  vtxZ[HIGHM]->Fill(var_vtxZ[i][vtx], weight);
	  vtxT[HIGHM]->Fill(vtx_t, weight);
	  dpt[HIGHM]->Fill(var_dpt[i][0], weight);
	  dptZoom1[HIGHM]->Fill(var_dpt[i][0], weight);
	  etSumPt[HIGHM]->Fill(et_sumPt, weight);
	  for (Int_t et=0; et<var_nTrack[i][0]; et++) {
	    if (var_et_vtx_dist[i][et]>0.01 ||
		var_et_quality[i][et] != 1)
	      continue;
	    etPtAfter[HIGHM]->Fill(var_et_pt[i][et], weight);
	    etEtaAfter[HIGHM]->Fill(var_et_eta[i][et], weight);
	  }
	  numVtxAfterCuts[HIGHM]->Fill(num_pvc-0.5, weight);
	}
      }
      invm[FULLM]->Fill(var_mass[i][0], weight);
      invm_over160GeV[FULLM]->Fill(var_mass[i][0], weight);
      ptpair[FULLM]->Fill(var_MuMupt[i][0], weight);
      pTSingleM[FULLM]->Fill(var_pt[i][pair_m], weight);
      pTSingleP[FULLM]->Fill(var_pt[i][pair_p], weight);
      ptpairZoom1[FULLM]->Fill(var_MuMupt[i][0], weight);
      ptpairZoom2[FULLM]->Fill(var_MuMupt[i][0], weight);
      ptpairZoom3[FULLM]->Fill(var_MuMupt[i][0], weight);
      etaSingleM[FULLM]->Fill(var_eta[i][pair_m], weight);
      etaSingleP[FULLM]->Fill(var_eta[i][pair_p], weight);
      etaPair[FULLM]->Fill(rap_pair, weight);
      acoplZoom1[FULLM]->Fill(acop, weight);
      acoplZoom2[FULLM]->Fill(acop, weight);
      acoplZoom3[FULLM]->Fill(acop, weight);
      acoplZoom4[FULLM]->Fill(acop, weight);
      acoplVsdpt[FULLM]->Fill(acop, var_dpt[i][0], weight);
      acoplVsptpair[FULLM]->Fill(acop, var_MuMupt[i][0], weight);
      vtxZ[FULLM]->Fill(var_vtxZ[i][vtx], weight);
      vtxT[FULLM]->Fill(vtx_t, weight);
      dpt[FULLM]->Fill(var_dpt[i][0], weight);
      dptZoom1[FULLM]->Fill(var_dpt[i][0], weight);
      etSumPt[FULLM]->Fill(et_sumPt, weight);
      for (Int_t et=0; et<var_nTrack[i][0]; et++) {
	if (var_et_vtx_dist[i][et]>0.01 ||
	    var_et_quality[i][et] != 1)
	  continue;
	etPtAfter[FULLM]->Fill(var_et_pt[i][et], weight);
	etEtaAfter[FULLM]->Fill(var_et_eta[i][et], weight);
      }
      numVtxAfterCuts[FULLM]->Fill(num_pvc-0.5, weight);

      if ((numTriggersMatched&&4)/4==1) {
	hTrig->Fill(2-0.5, weight);
	numTriggersMatched -= numTriggersMatched&&4;
      }
      if ((numTriggersMatched&&2)/2==1) {
	hTrig->Fill(1-0.5, weight);
	numTriggersMatched -= numTriggersMatched&&2;
      }
      if ((numTriggersMatched&&1)/1==1) {
	hTrig->Fill(0-0.5, weight);
	numTriggersMatched -= numTriggersMatched&&1;
      }
    } // end of events loop
  } // end of files loop
  TFile *out_file = new TFile(output_file, "RECREATE");
  cout << "Writing histograms data to output file..." << endl;
  hList->Write();
  cout << "Data succesfully written on " << output_file << "!" << endl;
  out_file->Close();
  
  if (drawAfterwards) {
    cout << "Now drawing the combined histograms" << endl;
    gROOT->ProcessLine(".x Draw.C");
  }
}
