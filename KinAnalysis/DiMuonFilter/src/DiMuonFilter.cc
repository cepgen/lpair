// -*- C++ -*-
//
// Package:    DiMuonFilter
// Class:      DiMuonFilter
// 
/**\class DiMuonFilter DiMuonFilter.cc Analysis/DiMuonFilter/src/DiMuonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Suzan Basegmez
//         Created:  Tue Apr 19 15:34:21 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//
// class declaration
//

class DiMuonFilter : public edm::EDFilter {
   public:
      explicit DiMuonFilter(const edm::ParameterSet&);
      ~DiMuonFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  
      bool   isHltMatched(const edm::Event&, const edm::EventSetup&, const std::string&, const reco::Muon&);
      bool   isHltPassed (const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
      double DR(double, double, double, double);


      // ----------member data ---------------------------

  edm::Handle<reco::BeamSpot> bsHandle;
  edm::ESHandle<TransientTrackBuilder> ttkb;
  edm::Handle<std::vector<reco::Muon> > muonCollectionHandle;
  
  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;
  std::string   processName_;
  std::vector<std::string> triggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  std::string fname;
  std::string triggerPath;

  bool checkHLT;
 
  double muonPTcut;
  double invMcut;


};

DiMuonFilter::DiMuonFilter(const edm::ParameterSet& iConfig)
{
  
   
  processName_       = iConfig.getParameter<std::string>("processName");
  triggerName_       = iConfig.getParameter<std::vector<std::string> >("triggerName");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  triggerEventTag_   = iConfig.getParameter<edm::InputTag>("triggerEvent");
  checkHLT           = iConfig.getParameter<bool>( "CheckHLT" );
  muonPTcut          = iConfig.getParameter<double>( "MuonPtCut" );
  invMcut            = iConfig.getParameter<double>( "InvMassCut" );
  

}


DiMuonFilter::~DiMuonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiMuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   bool KeepIT = false;

   if(checkHLT) {
     if (triggerName_.size() > 0) {
       iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
       if (!triggerResultsHandle_.isValid()) {
	 cout << "ERROR in getting TriggerResults product from Event!" << endl;
	 return false;
       }
       iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
       if (!triggerEventHandle_.isValid()) {
	 cout << "ERROR in getting TriggerEvent product from Event!" << endl;
	 return false;
       }
       cout << "triggerResultsHandle: " << (triggerResultsHandle_->size()) << endl;
       cout << "hltConfig: " << hltConfig_.size() << endl;
       assert(triggerResultsHandle_->size()==hltConfig_.size());
       
       bool hltEvOk = false;
       for (unsigned int i = 0; i < triggerName_.size(); i++) {
	 bool path = isHltPassed(iEvent,iSetup,triggerName_.at(i));
	 if (path && !hltEvOk) {
	   hltEvOk = true;
	   triggerPath = triggerName_.at(i);
	 }
       }
       if (!hltEvOk) return false;
     }
   }


   edm::InputTag bSLabel ("offlineBeamSpot");
   iEvent.getByLabel(bSLabel, bsHandle);
   reco::BeamSpot bSpot = *bsHandle;

   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);

   edm::InputTag muonLabel ("muons");
   iEvent.getByLabel(muonLabel, muonCollectionHandle);
   std::vector<reco::Muon> const & muon = *muonCollectionHandle;
   int nmuon = muon.size();
   int invM = 0;
   double angle3d_ = 0;
   
   // At least two muons in the event
   if (nmuon < 2 ) return false;
   
   int muon1 = 0; int muon2 = 0;
   for ( int j = 0; j < nmuon;  ++j ) {
     for ( int k = nmuon-1; k >= 0; --k) {
       
       // Opposite charge
       if (muon[j].charge()*muon[k].charge() >= 0) continue;
       
       // 2 Global muon
       if (!muon[j].isGlobalMuon() || !muon[k].isGlobalMuon()) continue;
       
       // both muons should have at least 10 tracker hits
       if (muon[j].innerTrack()->numberOfValidHits() < 10 || muon[k].innerTrack()->numberOfValidHits() < 10) continue;
       
       // both muons with pT > 20 GeV/c
       //  if (muon[j].pt() <= 20 || muon[k].pt() <= 20) continue;
        if (muon[j].pt() <= muonPTcut  || muon[k].pt() <= muonPTcut) continue;
	
       double muon1_d0 = muon[j].innerTrack()->dxy(bSpot.position());
       double muon2_d0 = muon[k].innerTrack()->dxy(bSpot.position());
       
       // at least one muon has to pass the following cuts
       // - dxy < 0.2 cm (dB < 0.2)
       // - muon global track chi2/ndof < 10 (globalTrack.normalizedChi2 < 10)
       // - at least one pixel hit (innerTrack.hitPattern.numberOfValidPixelHits >= 1)
       // - at least two muon stations in the fit (globalTrack.hitPattern.muonStationsWithValidHits >= 2)
       // - must be a tracker muon
       // - trigger matching to the single muon HLT path
       // muon[j].combinedMuon()->hitPattern().muonStationsWithValidHits() >= 2

       bool ms1 = false; bool ms2 = false;
       
       reco::TrackRef const gT1 = muon[j].globalTrack();
       reco::TrackRef const gT2 = muon[k].globalTrack();


       if(checkHLT) {
	 
       if (muon[j].isTrackerMuon() && muon[j].innerTrack()->hitPattern().numberOfValidPixelHits() >= 1 && muon[j].combinedMuon()->hitPattern().muonStationsWithValidHits() >= 2 && fabs(muon1_d0) < 0.2 && muon[j].globalTrack()->normalizedChi2() < 10 && isHltMatched(iEvent, iSetup, triggerPath, muon[j])) ms1 = true;

       if (muon[k].isTrackerMuon() && muon[k].innerTrack()->hitPattern().numberOfValidPixelHits() >= 1 && muon[k].combinedMuon()->hitPattern().muonStationsWithValidHits() >= 2 && fabs(muon2_d0) < 0.2 && muon[k].globalTrack()->normalizedChi2() < 10 && isHltMatched(iEvent, iSetup, triggerPath, muon[k])) ms2 = true;
       
       } else { 
	 
	 if ( muon[j].isTrackerMuon() && muon[j].innerTrack()->hitPattern().numberOfValidPixelHits() >= 1 && muon[j].combinedMuon()->hitPattern().muonStationsWithValidHits() >= 2 && fabs(muon1_d0) < 0.2 && muon[j].globalTrack()->normalizedChi2() < 10 ) ms1 = true;

	 if (muon[k].isTrackerMuon() && muon[k].innerTrack()->hitPattern().numberOfValidPixelHits() >= 1 && muon[k].combinedMuon()->hitPattern().muonStationsWithValidHits() >= 2 && fabs(muon2_d0) < 0.2 && muon[k].globalTrack()->normalizedChi2() < 10 ) ms2 = true;

       }
       
       if (!ms1 && !ms2) continue;


       // Cosmic rejection based on 3d angle
       double angle3d = acos(-muon[j].momentum().Dot(muon[k].momentum())/sqrt(muon[j].momentum().mag2())/sqrt(muon[k].momentum().mag2()));
       if (angle3d < 0.02) continue;

       
       // Isolation of Muons 
       if ((muon[j].isolationR03().sumPt/muon[j].pt()) >= 0.1 || (muon[k].isolationR03().sumPt/muon[k].pt()) >= 0.1) continue;

      
       // Common Vertex
       const reco::TrackRef& tk0 = muon[j].combinedMuon();
       const reco::TrackRef& tk1 = muon[k].combinedMuon();
       std::vector<reco::TransientTrack> ttv;
       ttv.push_back(ttkb->build(tk0));
       ttv.push_back(ttkb->build(tk1));
       KalmanVertexFitter kvf(true);
       TransientVertex tv = kvf.vertex(ttv);
       if(!tv.isValid() || tv.normalisedChiSquared() >= 10) continue;

      
       // Invariant Mass
       TLorentzVector P1(muon[j].momentum().x(), muon[j].momentum().y(), muon[j].momentum().z(), sqrt(muon[j].momentum().mag2()));
       TLorentzVector P2(muon[k].momentum().x(), muon[k].momentum().y(), muon[k].momentum().z(), sqrt(muon[k].momentum().mag2()));
       
       double Mz = (P1+P2).M();

       if (Mz > invM) {
	 invM = Mz;
	 muon1 = j;
	 muon2 = k;
	 angle3d_ = angle3d;
       }
     
     }  // for innerIter
   }  // for outerIter


   if(invM> invMcut ) KeepIT = true;

 
   return KeepIT;

}

// same check for isHltPassed +
// check if the muon is the one firing the HLT path
bool
DiMuonFilter::isHltMatched(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             const std::string& triggerName,
                             const reco::Muon& muon){ //const reco::Track& muon){

  if (triggerName == "") return true;

  bool isMatched=false;
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
  if (triggerIndex>=n) {
    cout << "DiMuonFilter::isHltMatched: path "
         << triggerName << " - not found!" << endl;
    return isMatched;
  }
  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error  = triggerResultsHandle_->error (triggerIndex);
  if (!wasRun) return isMatched;
  if (!accept) return isMatched;
  if ( error ) return isMatched;
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  assert (moduleIndex<m);
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
      const Vids& VIDS (triggerEventHandle_->filterIds (filterIndex));
      const Keys& KEYS (triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
      for (size_type i=0; i!=n; ++i) {
        const TriggerObject& TO(TOC[KEYS[i]]);
        double DPT=fabs(TO.pt()-muon.pt());
        if ( DR( TO.eta(),muon.eta(), TO.phi(),muon.phi() ) < 0.2
	     && fabs(TO.eta())<2.1
             && (DPT/muon.pt() < 1) ) isMatched=true;
      }
    }
  }
  return isMatched;
}



// ------------ method called once each job just before starting event loop  ------------
void 
DiMuonFilter::beginJob()
{
}


// Checking for Trigger Configuration 
void 
DiMuonFilter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  
  using namespace std;
  using namespace edm;

  if (triggerName_.size() > 0) {
    bool changed(true);
    if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
      //    if (changed) {

      const unsigned int n(hltConfig_.size());
      const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_.at(0)));
      if (triggerIndex>=n) {
        std::cout << "\n\nHLTEventAnalyzerAOD::analyze:"
                  << " HEY DUMBASS!!! TriggerName \"" << triggerPath
                  << "\" is NOT available in (new) config!" << std::endl << std::endl;
        std::cout << "Are you working in CMS? Mmmh I doubt it. Btw the "
                  << "available TriggerNames are: " << std::endl;
        hltConfig_.dump("Triggers");
        throw cms::Exception("DiMuonFilter")<<"I need to throw an exception to let you realize "
					    << "the trigger path name you want to check DOES NOT "
					    << "EXIST GODDAMMIT!!!! ";

        // }
        //hltConfig_.dump("Streams");
        //hltConfig_.dump("Datasets");
        //hltConfig_.dump("PrescaleTable");
        //hltConfig_.dump("ProcessPSet");
      }
    } else {
      cout << "HLTEventAnalyzerAOD::analyze:"
           << " config extraction failure with process name "
           << processName_ << endl;
      throw cms::Exception("DiMuonFilter")<<"Wrong processName_(\""<<processName_
					  <<"\"): please double check what you passed "
					  <<"in the python file... ";
    }
  }

}


bool
DiMuonFilter::isHltPassed(const edm::Event& iEvent,
                          const edm::EventSetup& iSetup,
                          const std::string& triggerName) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
  if (triggerIndex>=n) {
    //    cout << "DiMuonFilter::isHltPassed: path "
    //   << triggerName << " - not found!" << endl;
    return false;
  }
  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error  = triggerResultsHandle_->error (triggerIndex);
  bool isPassed=true;
  if (!wasRun) isPassed=false;
  if (!accept) isPassed=false;
  if ( error ) isPassed=false;
  return isPassed;
}


double
DiMuonFilter::DR(double eta1, double eta2,
                   double phi1, double phi2) {
  double diffEta = eta1 - eta2;
  double diffPhi = phi1 - phi2;
  double dr = sqrt(diffEta*diffEta + diffPhi*diffPhi);
  return dr;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiMuonFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonFilter);
