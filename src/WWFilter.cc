// -*- C++ -*-
//
// Faster dilepton filter based on AOD objects to pre-filter dilepton
// events for WW analysis needs
//

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

class WWFilter : public edm::EDFilter {
public:
  explicit WWFilter(const edm::ParameterSet&);
  ~WWFilter(){}
  
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  bool passedIso(const reco::GsfElectron* el);
  bool passedIso(const reco::Muon* mu);
  bool passedId(const reco::GsfElectron* el);
  bool passedId(const reco::Muon* mu);
  
  edm::InputTag muons_;
  edm::InputTag electrons_;
  std::vector<edm::InputTag> mets_;
  double minMuPt_;
  double minElePt_;
  double minMET_;
  double minMass_;
  bool applyEleId_;
  bool applyEleIso_;
  bool applyMuId_;
  bool applyMuIso_;
  int  prescale_;
};
bool WWFilter::passedIso(const reco::GsfElectron* el){
  if (!applyEleIso_) return true;
  return el->dr03TkSumPt()/el->pt()<0.2  && 
    el->dr03EcalRecHitSumEt()/el->pt()<0.2 && 
    el->dr03HcalTowerSumEt()/el->pt()<0.2;
}
bool WWFilter::passedIso(const reco::Muon* mu){
  if (!applyMuIso_) return true;
  return (mu->isolationR03().sumPt + mu->isolationR03().emEt +
	  mu->isolationR03().hadEt)/mu->pt()<1.0;
}
bool WWFilter::passedId(const reco::GsfElectron* el){
  if (!applyEleId_) return true; 
  // VBTF90
  if (el->isEB())
    return el->sigmaIetaIeta()<0.01 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.08 &&
      fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.007;
  else
    return el->sigmaIetaIeta()<0.03 && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.07 &&
      fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.009;
}

bool WWFilter::passedId(const reco::Muon* mu){
  if (!applyMuId_) return true; 
  return mu->isGlobalMuon() && mu->isTrackerMuon();
}

WWFilter::WWFilter(const edm::ParameterSet& iConfig):
  muons_          ( iConfig.getParameter<edm::InputTag>("muons") ),
  electrons_      ( iConfig.getParameter<edm::InputTag>("electrons") ),
  mets_           ( iConfig.getParameter<std::vector<edm::InputTag> >("mets") ),
  minMuPt_        ( iConfig.getParameter<double>("minMuPt") ),
  minElePt_       ( iConfig.getParameter<double>("minElePt") ),
  minMET_         ( iConfig.getParameter<double>("minMET") ),
  minMass_        ( iConfig.getParameter<double>("minMass") ),
  applyEleId_     ( iConfig.getParameter<bool>("applyEleId") ),
  applyEleIso_    ( iConfig.getParameter<bool>("applyEleIso") ),
  applyMuId_      ( iConfig.getParameter<bool>("applyMuId") ),
  applyMuIso_     ( iConfig.getParameter<bool>("applyMuIso") ),
  prescale_       ( iConfig.getParameter<int>("prescale") )
{
  std::string label = iConfig.getParameter<std::string>( "@module_label" );
  produces<bool>("passed").setBranchAlias("filter_"+label+"_passed");
  produces<bool>("run").setBranchAlias("filter_"+label+"_run");
  produces<int>("prescale").setBranchAlias("filter_"+label+"_prescale");
}

bool
WWFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<int>  filter_prescale(new int(prescale_));
  std::auto_ptr<bool> filter_passed(new bool(false));
  std::auto_ptr<bool> filter_run(new bool(false));
  if (prescale_<=0){
    iEvent.put(filter_prescale, "prescale");
    iEvent.put(filter_passed,   "passed");
    iEvent.put(filter_run,      "run");
    return false;
  }
  if ( iEvent.eventAuxiliary().event() % prescale_ != 0 ) {
    iEvent.put(filter_prescale, "prescale");
    iEvent.put(filter_passed,   "passed");
    iEvent.put(filter_run,      "run");
    return false;
  }
  *filter_run = true;
  if ( minMET_ > 0 ){
    bool passedMET = false;
    for (std::vector<edm::InputTag>::const_iterator tag = mets_.begin();
	 tag != mets_.end(); ++tag){
      edm::Handle<edm::View<reco::Candidate> > met;
      iEvent.getByLabel(*tag, met);
      if (met->front().pt()>minMET_){
	passedMET = true;
	break;
      }
    }
    if ( !passedMET ) {
      iEvent.put(filter_prescale, "prescale");
      iEvent.put(filter_passed,   "passed");
      iEvent.put(filter_run,      "run");
      return false;
    }
  }

  std::vector<const reco::RecoCandidate*> cands;
  
  edm::Handle<std::vector<reco::Muon> > muons;
  iEvent.getByLabel(muons_, muons);
  for ( std::vector<reco::Muon>::const_iterator muon = muons->begin();
	muon != muons->end(); ++muon )
    if ( muon->pt()>minMuPt_ && passedId(&*muon) && passedIso(&*muon) )
      cands.push_back(&*muon);
  
  edm::Handle<std::vector<reco::GsfElectron> > eles;
  iEvent.getByLabel(electrons_, eles);
  for ( std::vector<reco::GsfElectron>::const_iterator ele = eles->begin();
	ele != eles->end(); ++ele )
    if ( ele->pt()>minElePt_ && passedId(&*ele) && passedIso(&*ele) )
      cands.push_back(&*ele);

  if (cands.size()<2) {
    iEvent.put(filter_prescale, "prescale");
    iEvent.put(filter_passed,   "passed");
    iEvent.put(filter_run,      "run");
    return false;
  }

  if ( minMass_ > 0 ) {
    bool passedMass = false;
    for (unsigned int i=0; i < cands.size()-1; ++i)
      for (unsigned int j=i+1; j < cands.size(); ++j){
	if ((cands.at(i)->p4()+cands.at(j)->p4()).mass2() <= 0 ) continue;
	if ((cands.at(i)->p4()+cands.at(j)->p4()).mass() < minMass_ ) continue;
	passedMass = true;
      }
    if (!passedMass){
      iEvent.put(filter_prescale, "prescale");
      iEvent.put(filter_passed,   "passed");
      iEvent.put(filter_run,      "run");
      return false;
    }
  }
  *filter_passed = true;
  iEvent.put(filter_prescale, "prescale");
  iEvent.put(filter_passed,   "passed");
  iEvent.put(filter_run,      "run");
  return true;
}
//define this as a plug-in
DEFINE_FWK_MODULE(WWFilter);
