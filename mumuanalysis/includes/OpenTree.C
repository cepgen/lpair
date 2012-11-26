#include "TTree.h"
#include <iostream>
using namespace std;
#define MAXNUM 12 // the maximum of sources for TTrees
#define NUMTRIG 7

// Bunch crossing
Int_t var_bx[MAXNUM][1], var_run[MAXNUM][1], var_ls[MAXNUM][1], var_event[MAXNUM][1];
//FIXME -> Double_t? BUG?
Int_t hlt_d[NUMTRIG][MAXNUM][1];
//Int_t hlt_d0[MAXNUM][1], hlt_d1[MAXNUM][1], hlt_d2[MAXNUM][1], hlt_d3[MAXNUM][1], hlt_d4[MAXNUM][1], hlt_d5[MAXNUM][1];
Int_t techBit[MAXNUM][128];
//Int_t var_idA[MAXNUM][100],var_idB[MAXNUM][100], var_idC[MAXNUM][100], var_idD[MAXNUM][100], var_idE[MAXNUM][100];
// RecoTrack
Int_t var_nTrack[MAXNUM][1];
Int_t var_nTrackQual[MAXNUM][1];
Double_t var_et_p[MAXNUM][2000];
Double_t var_et_pt[MAXNUM][2000]; // As many as track candidates
Double_t var_et_eta[MAXNUM][2000];
Double_t var_et_phi[MAXNUM][2000];
Double_t var_et_vtx_z[MAXNUM][2000];
Double_t var_et_vtx_dist[MAXNUM][2000];
Double_t var_et_quality[MAXNUM][2000];
// Primary Vtx
Int_t var_nvtx[MAXNUM];
Double_t var_vtxX[MAXNUM][100], var_vtxY[MAXNUM][100], var_vtxZ[MAXNUM][100]; // As many as prim vtx candidates
Double_t var_vertexChi2[MAXNUM][100], var_vertexNdf[MAXNUM][100];
Double_t var_vtxTrack[MAXNUM][100];
Int_t var_vtxTrack_int[MAXNUM][100];
Bool_t var_vtxTrack_is_int[MAXNUM];
Int_t var_vtxmumu[MAXNUM][100];
Int_t var_vtxhasmumu[MAXNUM][100];

// CaloTowers 
/*Int_t var_ncalo[MAXNUM][1];
  Int_t var_caloId[MAXNUM][2000];
  Double_t var_caloEn[MAXNUM][2000], var_caloTime[MAXNUM][2000], var_calodR[MAXNUM][2000];
  Double_t var_caloEta[MAXNUM][2000], var_caloPhi[MAXNUM][2000];
  Double_t var_etmiss[MAXNUM][1];
  Int_t var_tower[MAXNUM][1];*/
/*Double_t var_caloEmE[MAXNUM][2000];
  Double_t var_caloHadE[MAXNUM][2000];
  Double_t var_caloZ[MAXNUM][2000];*/

// MuMu kinematics
Double_t var_mass[MAXNUM][5], var_dpt[MAXNUM][5], var_dphi[MAXNUM][5];
Double_t var_MuMupt[MAXNUM][5];
Double_t var_MuMuvtxX[MAXNUM][100], var_MuMuvtxY[MAXNUM][100], var_MuMuvtxZ[MAXNUM][100];
Int_t var_MuMuvtxValid[MAXNUM][100];

// Single muon properties
Int_t var_global[MAXNUM][100], var_tracker[MAXNUM][100],var_standalone[MAXNUM][100]; // Muon reconstruction algorithm
Int_t var_charge[MAXNUM][100]; // muon charge
Double_t var_normChi[MAXNUM][100];
Int_t var_matches[MAXNUM][100];
Int_t var_nhitsMuon[MAXNUM][100]; // number of hits in the muon chambers
Int_t var_nhitsTrack[MAXNUM][100]; // number of hits in the tracker
Int_t var_nhitsPixel[MAXNUM][100]; // number of hits in the pixel chambers
Double_t var_pt[MAXNUM][100], var_phi[MAXNUM][100], var_eta[MAXNUM][100];
Double_t var_p[MAXNUM][100];
Double_t var_px[MAXNUM][100], var_py[MAXNUM][100], var_pz[MAXNUM][100];
Int_t var_MuCand_tightID[MAXNUM][100];

Int_t var_Pair[MAXNUM][2];

// Efficiency
Double_t var_eff[MAXNUM][100];

// ZDC
/*Int_t var_nZDC[MAXNUM][1];
  Double_t var_zdcEmMinus[MAXNUM][1], var_zdcHadMinus[MAXNUM][1];
  Double_t var_zdcEmPlus[MAXNUM][1], var_zdcHadPlus[MAXNUM][1];
  Double_t var_zdcTime[MAXNUM][5000];
  Double_t var_zdcE[MAXNUM][5000];
  Int_t var_zdcsection[MAXNUM][5000];*/

// Castor
/*Int_t var_nCastor[MAXNUM][1];
  Double_t var_CastorE[MAXNUM][1000];
  Double_t var_CastorEta[MAXNUM][1000];
  Double_t var_CastorPhi[MAXNUM][1000];
  Double_t var_CastorRecHit[MAXNUM][1];*/

// Pileup reweighting
Double_t var_PUweight[MAXNUM][1];
Double_t var_PUweightA[MAXNUM][1];
Bool_t var_issetPUweightA[MAXNUM][1];
Double_t var_PUweightB[MAXNUM][1];
Bool_t var_issetPUweightB[MAXNUM][1];

int OpenTree(TTree *tree, int index, bool isMC) {
  if (!isMC) {
    tree->SetBranchAddress("BX",var_bx[index]);
    tree->SetBranchAddress("LumiSection",var_ls[index]);
    tree->SetBranchAddress("EventNum",var_event[index]);
  }
  tree->SetBranchAddress("Run",var_run[index]);

  //////////////////////////////////////////////////////////////////////////////
  //                       High Level Trigger Decisions                       //
  //////////////////////////////////////////////////////////////////////////////
  
  tree->SetBranchAddress("HLT_DoubleMu4Acoplanarity",hlt_d[0][index]);
  tree->SetBranchAddress("HLT_DoubleMu5Acoplanarity",hlt_d[1][index]);
  tree->SetBranchAddress("HLT_DoubleMu6Acoplanarity",hlt_d[2][index]);
  tree->SetBranchAddress("HLT_DoubleMu7Acoplanarity",hlt_d[3][index]);
  tree->SetBranchAddress("HLT_Mu13Mu8",hlt_d[4][index]);
  tree->SetBranchAddress("HLT_Mu17Mu8",hlt_d[5][index]);
  tree->SetBranchAddress("HLT_DoubleMu7",hlt_d[6][index]);
  
  //////////////////////////////////////////////////////////////////////////////
  //                           L1 Trigger Decisions                           //
  //////////////////////////////////////////////////////////////////////////////
  
  tree->SetBranchAddress("L1TechnicalTriggers",techBit[index]);
  //cout << “Hahaha!” << endl;
  /*
    tree->SetBranchAddress("MuonCand_tmlsOptLowPtloosemuonid",var_idA[index]);
    tree->SetBranchAddress("MuonCand_tmlsAngloosemuonid",var_idB[index]);
    tree->SetBranchAddress("MuonCand_tmlsAngtightmuonid",var_idC[index]);
    tree->SetBranchAddress("MuonCand_tmosAngloosemuonid",var_idD[index]);
    tree->SetBranchAddress("MuonCand_tmosAngtightmuonid",var_idE[index]);
  */
  tree->SetBranchAddress("MuonCand_tightID",var_MuCand_tightID[index]);

  tree->SetBranchAddress("nTrackCand",var_nTrack[index]);
  tree->SetBranchAddress("nQualityTrackCand",var_nTrackQual[index]);
  tree->SetBranchAddress("TrackCand_p",var_et_p[index]);
  tree->SetBranchAddress("TrackCand_pt",var_et_pt[index]);
  tree->SetBranchAddress("TrackCand_eta",var_et_eta[index]);
  tree->SetBranchAddress("TrackCand_phi",var_et_phi[index]);
  tree->SetBranchAddress("TrackCand_vtxZ",var_et_vtx_z[index]);
  tree->SetBranchAddress("TrackCand_vtxdxyz",var_et_vtx_dist[index]);
  tree->SetBranchAddress("TrackCand_purity",var_et_quality[index]);

  //////////////////////////////////////////////////////////////////////////////
  //                        Primary vertex information                        //
  //////////////////////////////////////////////////////////////////////////////
  
  tree->SetBranchAddress("nPrimVertexCand",&(var_nvtx[index]));
  tree->SetBranchAddress("PrimVertexCand_x",var_vtxX[index]);
  tree->SetBranchAddress("PrimVertexCand_y",var_vtxY[index]);
  tree->SetBranchAddress("PrimVertexCand_z",var_vtxZ[index]);
  tree->SetBranchAddress("PrimVertexCand_chi2",var_vertexChi2[index]);
  tree->SetBranchAddress("PrimVertexCand_ndof",var_vertexNdf[index]);

  if (tree->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack[index]) == -2) {
    var_vtxTrack_is_int[index] = true;
    cout << "[DEBUG] ----> PrimVertexCand_tracks leaf is an Int_t" << endl;
    tree->SetBranchAddress("PrimVertexCand_tracks",var_vtxTrack_int[index]);
    /*cout << "----> " << var_vtxTrack_int[index][0] << endl;
      memcpy(var_vtxTrack[index],var_vtxTrack_int[index],sizeof(var_vtxTrack_int)+1);
      cout << "----> " << var_vtxTrack[index][0] << endl;*/
  }
  //else memcpy(var_vtxTrack[index],var_vtxTrack_dbl[index],sizeof(var_vtxTrack_dbl)+1);
  //cout << "----> [FINAL] " << var_vtxTrack[index][0] << endl;
  tree->SetBranchAddress("PrimVertexCand_mumuExactlyTwoTracks",var_vtxmumu[index]); 
  tree->SetBranchAddress("PrimVertexCand_mumuTwoTracks",var_vtxhasmumu[index]);

  tree->SetBranchAddress("MuMu_Kalmanvtxx",var_MuMuvtxX[index]);
  tree->SetBranchAddress("MuMu_Kalmanvtxy",var_MuMuvtxY[index]);
  tree->SetBranchAddress("MuMu_Kalmanvtxz",var_MuMuvtxZ[index]);
  tree->SetBranchAddress("MuMu_Kalmanvtxisvalid",var_MuMuvtxValid[index]);
  
  //////////////////////////////////////////////////////////////////////////////
  //                        Calorimeters information                          //
  //////////////////////////////////////////////////////////////////////////////
  /*  
      tree->SetBranchAddress("nCaloCand",var_ncalo[index]);
      tree->SetBranchAddress("CaloTower_ID",var_caloId[index]);
      tree->SetBranchAddress("CaloTower_e",var_caloEn[index]);
      tree->SetBranchAddress("CaloTower_t",var_caloTime[index]);
      tree->SetBranchAddress("CaloTower_dr",var_calodR[index]);
      tree->SetBranchAddress("CaloTower_eta",var_caloEta[index]);
      tree->SetBranchAddress("CaloTower_phi",var_caloPhi[index]);
      tree->SetBranchAddress("Etmiss",var_etmiss[index]);
      tree->SetBranchAddress("nExtraCaloTowersE5",var_tower[index]);
  */
  //////////////////////////////////////////////////////////////////////////////
  //                            Dimuon kinematics                             //
  //////////////////////////////////////////////////////////////////////////////
  
  tree->SetBranchAddress("MuMu_mass",var_mass[index]);
  tree->SetBranchAddress("MuMu_pt",var_MuMupt[index]);
  tree->SetBranchAddress("MuMu_dpt",var_dpt[index]);
  tree->SetBranchAddress("MuMu_dphi",var_dphi[index]);

  tree->SetBranchAddress("MuonCand_isglobal",var_global[index]);
  tree->SetBranchAddress("MuonCand_istracker",var_tracker[index]);
  tree->SetBranchAddress("MuonCand_isstandalone",var_standalone[index]);
  tree->SetBranchAddress("MuonCand_validtrackhits",var_nhitsTrack[index]);
  tree->SetBranchAddress("MuonCand_charge",var_charge[index]);
  tree->SetBranchAddress("MuonCand_normchi2",var_normChi[index]);
  tree->SetBranchAddress("MuonCand_matches",var_matches[index]);
  tree->SetBranchAddress("MuonCand_validmuonhits",var_nhitsMuon[index]);
  tree->SetBranchAddress("MuonCand_validpixelhits",var_nhitsPixel[index]);

  tree->SetBranchAddress("MuonCand_pt",var_pt[index]);
  tree->SetBranchAddress("MuonCand_pz",var_pz[index]);
  tree->SetBranchAddress("MuonCand_phi",var_phi[index]);
  tree->SetBranchAddress("MuonCand_eta",var_eta[index]);
  tree->SetBranchAddress("MuonCand_p",var_p[index]);
  tree->SetBranchAddress("MuonCand_px",var_px[index]);
  tree->SetBranchAddress("MuonCand_py",var_py[index]);
  tree->SetBranchAddress("MuonCand_pz",var_pz[index]);

  tree->SetBranchAddress("MuonPairCand",var_Pair[index]);
  if (isMC)
    tree->SetBranchAddress("MuonCand_efficiency",var_eff[index]);

  //////////////////////////////////////////////////////////////////////////////
  //                             ZDC information                              //
  //////////////////////////////////////////////////////////////////////////////
  
  /*
    tree->SetBranchAddress("nZDChitCand",var_nZDC[index]);
    tree->SetBranchAddress("ZDCsumEMminus",var_zdcEmMinus[index]);
    tree->SetBranchAddress("ZDCsumHADminus",var_zdcHadMinus[index]);
    tree->SetBranchAddress("ZDCsumEMplus",var_zdcEmPlus[index]);
    tree->SetBranchAddress("ZDCsumHADplus",var_zdcHadPlus[index]);
    tree->SetBranchAddress("ZDChit_time",var_zdcTime[index]);
    tree->SetBranchAddress("ZDChit_energy",var_zdcE[index]);
    tree->SetBranchAddress("ZDChit_section",var_zdcsection[index]);
  */

  //////////////////////////////////////////////////////////////////////////////
  //                           Castor information                             //
  //////////////////////////////////////////////////////////////////////////////
  
  /*
    tree->SetBranchAddress("nCastorTowerCand",var_nCastor[index]);
    tree->SetBranchAddress("CastorTower_e",var_CastorE[index]);
    tree->SetBranchAddress("CastorTower_eta",var_CastorEta[index]);
    tree->SetBranchAddress("CastorTower_phi",var_CastorPhi[index]);
    tree->SetBranchAddress("CASTORsumRecHitsE",var_CastorRecHit[index]);
  */

  //////////////////////////////////////////////////////////////////////////////
  //                             Other information                            //
  //////////////////////////////////////////////////////////////////////////////

  if (tree->SetBranchAddress("Weight3D",var_PUweight[index])<0) {
    cout << "[DEBUG] ----> Weight3D leaf is not set : replacing it with value 1." << endl;
    var_PUweight[index][0] = 1.;
  }

  if (tree->SetBranchAddress("Weight3D_2011A",var_PUweightA[index])<0) {
    cout << "[DEBUG] ----> [2011A] Weight3D leaf is not set : replacing it with value 1." << endl;
    var_PUweightA[index][0] = 1.;
    var_issetPUweightA[index][0] = false;
  }
  else var_issetPUweightA[index][0] = true;

  if (tree->SetBranchAddress("Weight3D_2011B",var_PUweightB[index])<0) {
    cout << "[DEBUG] ----> [2011B] Weight3D leaf is not set : replacing it with value 1." << endl;
    var_PUweightB[index][0] = 1.;
    var_issetPUweightB[index][0] = false;
  }
  else var_issetPUweightB[index][0] = true;

  return 0;

}
