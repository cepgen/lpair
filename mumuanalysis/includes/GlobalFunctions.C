#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

std::vector<bool> isMC;
std::vector<TString> datasetName;
std::vector<TString> filename;
std::vector<Double_t> numEvents;
std::vector<Double_t> datasetWeight;
std::vector<TString> quantities;
std::vector<TString> quantitiesNames;
std::vector<TString> quantitiesUnits;
std::vector<TString> quantitiesTitles;
std::vector<bool> quantitiesPlotLog;
std::vector<TString> dists;
std::vector<TString> distsTitles;
std::vector<Double_t> distsWeights;
std::vector<Int_t> distsColours;
std::vector<Int_t> distsIsMC;

// Muon efficiencies
std::vector<Double_t> etaLowA, etaHighA;
std::vector<Double_t> pTLowA, pTHighA;
std::vector<Double_t> MuEffA, errLowA, errHighA;
std::vector<Double_t> etaLowB, etaHighB;
std::vector<Double_t> pTLowB, pTHighB;
std::vector<Double_t> MuEffB, errLowB, errHighB;

// Trigger efficiencies
std::vector<Double_t> TEetaLow17, TEetaHigh17;
std::vector<Double_t> TEpTLow17, TEpTHigh17;
std::vector<Double_t> TEff17, TEerrLow17, TEerrHigh17;
std::vector<Double_t> TEetaLow8, TEetaHigh8;
std::vector<Double_t> TEpTLow8, TEpTHigh8;
std::vector<Double_t> TEff8, TEerrLow8, TEerrHigh8;

Int_t AddFile(TString name_, TString filename_, bool isMC_, Double_t datasetWeight_=1.) {
  datasetName.push_back(name_);
  filename.push_back(filename_);
  isMC.push_back(isMC_);
  TFile *file = new TFile(filename_);
  TTree *tree = (TTree*)(file->Get("ntp1"));
  numEvents.push_back(tree->GetEntries());
  datasetWeight.push_back(datasetWeight_);
  
  return 0;
}

Int_t AddQuantity(TString name_, TString title_, TString description_, TString unit_, Bool_t plotlog_=false) {
  quantities.push_back(name_);
  quantitiesNames.push_back(description_);
  quantitiesUnits.push_back(unit_);
  quantitiesTitles.push_back(title_);
  quantitiesPlotLog.push_back(plotlog_);

  return 0;
}

Int_t AddDistribution(TString name_, TString title_, Double_t weight_, Int_t colour_, Int_t isMC_=1) {
  if (false) {
    cout << "[DEBUG] : AddDistribution (" << endl;
    cout << "            name : " << name_ << endl;
    cout << "            title : " << title_ << endl;
    cout << "            weight : " << weight_ << endl;
    cout << "            colour : " << colour_ << endl;
    cout << "            isMC : " << isMC_ << endl;
    cout << "          )" << endl;
  }
  dists.push_back(name_);
  distsTitles.push_back(title_);
  distsWeights.push_back(weight_);
  distsColours.push_back(colour_);
  distsIsMC.push_back(isMC_);

  return 0;
}

bool PassesTrigger(Int_t triggerbit, Int_t runNumber=1) {
  bool pass = false;

  if(triggerbit == 1) pass = true;
  //if(runNumber>=147146) pass = true;
  if (runNumber == 163869)
    std::cout << "Run2011A --> Run2011B" << endl;
  return pass;
}


void GetMuonEfficiencyCorrections() {
  Double_t etaLowtmp, etaHightmp, pTLowtmp, pTHightmp, efftmp, errorLowtmp, errorHightmp;
  
  TString muEff_Run2011A("includes/TIGHT_nL8_2011A_ratio.txt");
  TString muEff_Run2011B("includes/TIGHT_nL8_2011B_ratio.txt");

  std::ifstream muEffA(muEff_Run2011A);
  std::cout << "===== Run2011A =====" << std::endl;
  while(!muEffA.eof()) {
    muEffA >> etaLowtmp >> etaHightmp 
	   >> pTLowtmp >> pTHightmp 
	   >> efftmp >> errorLowtmp >> errorHightmp;

    // Run2011A
    etaLowA.push_back(etaLowtmp); etaHighA.push_back(etaHightmp);
    pTLowA.push_back(pTLowtmp); pTHighA.push_back(pTHightmp);
    MuEffA.push_back(efftmp);
    errLowA.push_back(etaLowtmp); errHighA.push_back(errorHightmp);
    std::cout << "Efficiency for #eta in [" << etaLowtmp << "; " << etaHightmp << "] "
	      << "and pT in [" << pTLowtmp << "; " << pTHightmp << "] : "
	      << efftmp << std::endl;
  }
  std::ifstream muEffB(muEff_Run2011B);
  std::cout << "===== Run2011B =====" << std::endl;
  while(!muEffB.eof()) {
    muEffB >> etaLowtmp >> etaHightmp 
	   >> pTLowtmp >> pTHightmp 
	   >> efftmp >> errorLowtmp >> errorHightmp;

    // Run2011B
    etaLowB.push_back(etaLowtmp); etaHighB.push_back(etaHightmp);
    pTLowB.push_back(pTLowtmp); pTHighB.push_back(pTHightmp);
    MuEffB.push_back(efftmp);
    errLowA.push_back(etaLowtmp); errHighB.push_back(errorHightmp);
    std::cout << "Efficiency for #eta in [" << etaLowtmp << "; " << etaHightmp << "] "
	      << "and pT in [" << pTLowtmp << "; " << pTHightmp << "] : "
	      << efftmp << std::endl;
  }
}

void GetTriggerEfficiencyCorrections() {
  Double_t etaLowtmp, etaHightmp, pTLowtmp, pTHightmp, efftmp, errorLowtmp, errorHightmp;

  TString mu17_file("includes/Mu17_2011_ratio.txt");
  TString mu8_file("includes/Mu8_2011_ratio.txt");

  std::cout << "===== Trigger Efficiencies (Mu17) =====" << std::endl;
  std::ifstream Eff17(mu17_file);
  while(!Eff17.eof()) {
    Eff17 >> etaLowtmp >> etaHightmp
	  >> pTLowtmp >> pTHightmp
	  >> efftmp >> errorLowtmp >> errorHightmp;

    TEetaLow17.push_back(etaLowtmp); TEetaHigh17.push_back(etaHightmp);
    TEpTLow17.push_back(pTLowtmp); TEpTHigh17.push_back(pTHightmp);
    TEff17.push_back(efftmp);
    TEerrLow17.push_back(etaLowtmp); TEerrHigh17.push_back(errorHightmp);
    std::cout << "Efficiency for #eta in [" << etaLowtmp << "; " << etaHightmp << "] "
	      << "and pT in [" << pTLowtmp << "; " << pTHightmp << "] : "
	      << efftmp << std::endl;
  }
  std::cout << "===== Trigger Efficiencies (Mu8) =====" << std::endl;
  std::ifstream Eff8(mu8_file);
  while(!Eff8.eof()) {
    Eff8 >> etaLowtmp >> etaHightmp
	  >> pTLowtmp >> pTHightmp
	  >> efftmp >> errorLowtmp >> errorHightmp;

    TEetaLow8.push_back(etaLowtmp); TEetaHigh8.push_back(etaHightmp);
    TEpTLow8.push_back(pTLowtmp); TEpTHigh8.push_back(pTHightmp);
    TEff8.push_back(efftmp);
    TEerrLow8.push_back(etaLowtmp); TEerrHigh8.push_back(errorHightmp);
    std::cout << "Efficiency for #eta in [" << etaLowtmp << "; " << etaHightmp << "] "
	      << "and pT in [" << pTLowtmp << "; " << pTHightmp << "] : "
	      << efftmp << std::endl;
  }
}

Double_t VertexSeparation(Int_t nvtx, Double_t *vtxtrks, Double_t *vtxz, Int_t *ismumuvtx) {
  Double_t closestvtx = 9999.;
  std::cout << "Tracks on vertex : " << vtxtrks << std::endl;
  if(nvtx == 1) return 9999.;

  for(Int_t i=0; i<nvtx; ++i) {
    if (ismumuvtx[i] == 0) continue;
	
    for(Int_t j=0; j<nvtx && (j!=i); j++) {
    //for (Int_t j=0; j<i; ++j) {
      if (fabs(vtxz[i]-vtxz[j])<fabs(closestvtx))
        closestvtx = vtxz[i]-vtxz[j];
    }
  }
  return closestvtx;
}
Double_t VertexSeparation(Int_t nvtx, Int_t *vtxtrks, Double_t *vtxz, Int_t *ismumuvtx) {
  Double_t closestvtx = 9999.;
  std::cout << "Tracks on vertex : " << vtxtrks << std::endl;
  if(nvtx == 1) return 9999.;

  for(Int_t i=0; i<nvtx; ++i) {
    if (ismumuvtx[i] == 0) continue;
	
    for(Int_t j=0; j<nvtx && (j!=i); j++) {
      if (fabs(vtxz[i]-vtxz[j])<fabs(closestvtx))
        closestvtx = vtxz[i]-vtxz[j];
    }
  }
  return closestvtx;
}

typedef enum {
  NO_CUT,
  TRIGGER,
  INVM,
  PT,
  ETA,
  CHISQ_VTX,
  SEPARATION_VTX,
  DIMUON_VTX,
  Z_VTX,
  NUM_VTX_TRACKS,
  TRANSV_VTX,
  TRACK_HITS,
  MU_HITS,
  CHISQ,
  GLOB_TRACK,
  DPT,
  DPHI,
  OPEN_ANGLE,
  EXCLU_TRACK,
  ZDC
} cuts;
