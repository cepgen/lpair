// -*- C++ -*-
//
// Package:    KinAnalysis
// Class:      KinAnalysis
// 
/**\class KinAnalysis KinAnalysis.cc KinAnalysis/KinAnalysis/src/KinAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Sat Apr 16 13:18:40 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// Reconstruction
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// Colins-Soper reference frame
//#include "PhysicsTools/TagAndProbe/interface/ColinsSoperVariables.h"
#include "../interface/ColinsSoperVariables.h" //FIXME

// ROOT
#include "TH1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


//
// class declaration
//

class KinAnalysis : public edm::EDAnalyzer {
public:
  explicit KinAnalysis(const edm::ParameterSet&);
  ~KinAnalysis();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  TFile* file;
  TH1D *hist_etam, *hist_etap, *hist_eta, *hist_mass, *hist_asy, *hist_ell, *hist_cos_theta_cs;
  TH1D *hist_ptm, *hist_ptp, *hist_pt;
  float M;
  int ch1,ch2;
  double x, etamax;
  double sum_ell_inf, sum_ell_med, sum_ell_sup;
  float ell;
  double tot_ell;
  double zstar, zmax;
  double sum_asy_inf, sum_asy_med, sum_asy_sup;
  double sum_asy_cs_inf, sum_asy_cs_med, sum_asy_cs_sup;
  float asy;
  double tot_asy;
  float asy_cs;
  double tot_asy_cs;
  float z, z_cs;
  
  bool workOnGenLevel;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KinAnalysis::KinAnalysis(const edm::ParameterSet& iConfig)
  : x(1.), etamax(2.5), zstar(.5), zmax(1.)
{
   //now do what ever initialization is needed

  //workOnGenLevel = iConfig.getParameter<bool>('workOnGenLevel');

  edm::Service<TFileService> fs;
  hist_etam=fs->make<TH1D>("eta_m","#eta for #mu^{-}",100,-4.,4.);
  hist_etap=fs->make<TH1D>("eta_p","#eta for #mu^{+}",100,-4.,4.);
  hist_eta=fs->make<TH1D>("eta","#eta for #mu^{#pm}",100,-4.,4.);
  hist_ptm=fs->make<TH1D>("pt_m","p_{T} for #mu^{-}",400,0.,400.);
  hist_ptp=fs->make<TH1D>("pt_p","p_{T} for #mu^{+}",400,0.,400.);
  hist_pt=fs->make<TH1D>("pt","p_{T} for #mu^{#pm}",400,0.,400.);
  hist_mass=fs->make<TH1D>("mass","M_{#mu#mu}",500,0.,500.);
  hist_asy=fs->make<TH1D>("asym","asym",1000,-10.,10.);
  hist_ell=fs->make<TH1D>("ellip","ellip",1000,-10.,10.);
  hist_cos_theta_cs=fs->make<TH1D>("cos_theta","cos_theta",100,-1.,1.);
}


KinAnalysis::~KinAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
KinAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


  //else {
    edm::Handle<reco::MuonCollection> muHandle;
    iEvent.getByLabel("muons",muHandle);
    reco::MuonCollection muColl = *muHandle.product();
    unsigned int NMuons = muColl.size();

    if (NMuons >= 2) {
      int ind[2]={-1,-1};

      // We only keep the highest invariant mass pair
      M = 0;
      float Mtemp=0;
      for(unsigned int i=0; i<NMuons; i++) {
        reco::MuonRef m1 = reco::MuonRef(muHandle,i);
        TLorentzVector P1(m1->momentum().x(),m1->momentum().y(),m1->momentum().z(),sqrt(m1->momentum().mag2()));
        for(unsigned int j=0; j<i; j++) {
      	 reco::MuonRef m2 = reco::MuonRef(muHandle,j);
      	 TLorentzVector P2(m2->momentum().x(),m2->momentum().y(),m2->momentum().z(),sqrt(m2->momentum().mag2()));
      	 Mtemp = (P1+P2).M();
      	 ch1 = m1->charge();
      	 ch2 = m2->charge();
      	 //std::cout << "muon 1: " << ch1 << ", muon 2: " << ch2 << std::endl;
      	 if ((Mtemp>M) && (ch1*ch2==-1)) {
      	   M = Mtemp;
      	   if(ch1==-1) {
      	     ind[0] = i;
      	     ind[1] = j;
      	   }
      	   else {
      	     ind[0] = j;
      	     ind[1] = i;
      	   }
      	 }
        }
      }
      if (ind[0]!=-1 && ind[1]!=-1) { // We have found the opposite charge/highest invariant mass pair
        //std::cout << "mu m is " << ind[0] << ", mu p is " << ind[1] << std::endl;
        
        reco::MuonRef mum = reco::MuonRef(muHandle,ind[0]);
        reco::MuonRef mup = reco::MuonRef(muHandle,ind[1]);

        if ((mum->eta())>=(-etamax)) {
      	 if (mum->eta()<(-x)) sum_ell_inf++;
      	 else if (mum->eta()<x) sum_ell_med++;
      	 else if (mum->eta()<=etamax) sum_ell_sup++;
        }

        hist_ptm->Fill(mum->pt());
        hist_ptp->Fill(mup->pt());
        hist_pt->Fill(mum->pt());
        hist_pt->Fill(mup->pt());

        hist_etam->Fill(mum->eta());
        hist_etap->Fill(mup->eta());
        hist_eta->Fill(mum->eta());
        hist_eta->Fill(mup->eta());

        TLorentzVector Pm(mum->momentum().x(),mum->momentum().y(),mum->momentum().z(),sqrt(mum->momentum().mag2()));
        TLorentzVector Pp(mup->momentum().x(),mup->momentum().y(),mup->momentum().z(),sqrt(mup->momentum().mag2()));
        double result_cs[3];
        calCSVariables(Pm, Pp, result_cs, false);
        hist_mass->Fill((Pm+Pp).M());
        hist_cos_theta_cs->Fill(result_cs[0]);
        z = tanh((mum->eta()-mup->eta())/2);
        if (z>=(-zmax)) {
      	 if (z<(-zstar)) sum_asy_inf++;
      	 else if (z<zstar) sum_asy_med++;
      	 else if (z<zmax) sum_asy_sup++;
        }
        z_cs = result_cs[0];
        if (z_cs>=(-zmax)) {
      	 if (z_cs<(-zstar)) sum_asy_cs_inf++;
      	 else if (z_cs<zstar) sum_asy_cs_med++;
      	 else if (z_cs<zmax) sum_asy_cs_sup++;
        }
      }
    }
  //}
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
KinAnalysis::beginJob()
{
  sum_ell_inf=0;
  sum_ell_med=0;
  sum_ell_sup=0;
  sum_asy_inf=0;
  sum_asy_med=0;
  sum_asy_sup=0;
  sum_asy_cs_inf=0;
  sum_asy_cs_med=0;
  sum_asy_cs_sup=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KinAnalysis::endJob() {
  tot_ell = sum_ell_inf+sum_ell_med+sum_ell_sup;
  ell = (sum_ell_med-sum_ell_inf-sum_ell_sup)/tot_ell;
  std::cout << "Total for ellipcity: "<< tot_ell << std::endl;
  std::cout << "Ellipcity: " << ell << std::endl;
  std::cout << "==================================" << std::endl;
  tot_asy = sum_asy_inf+sum_asy_med+sum_asy_sup;
  asy = (sum_asy_med-sum_asy_inf-sum_asy_sup)/tot_asy;
  std::cout << "Total for asymmetry (using eta difference): "<< tot_asy << std::endl;
  std::cout << "Center-edge asymmetry (using eta difference): " << asy << std::endl;
  tot_asy_cs = sum_asy_cs_inf+sum_asy_cs_med+sum_asy_cs_sup;
  asy_cs = (sum_asy_cs_med-sum_asy_cs_inf-sum_asy_cs_sup)/tot_asy_cs;
  std::cout << "Total for asymmetry (using CS frame): "<< tot_asy_cs << std::endl;
  std::cout << "Center-edge asymmetry (using CS frame): " << asy_cs << std::endl;
  hist_ell->Fill(ell);
  hist_asy->Fill(asy);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinAnalysis);
