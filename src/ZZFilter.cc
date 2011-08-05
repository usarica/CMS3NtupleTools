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

class ZZFilter : public edm::EDFilter {
public:
  explicit ZZFilter(const edm::ParameterSet&);
  ~ZZFilter(){}
  
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  edm::InputTag muons_;
  edm::InputTag electrons_;
  double minMuPt_;
  double minElePt_;
  double minLeadingLepPt_;
  double minNLeptons_;
  int  prescale_;
};
ZZFilter::ZZFilter(const edm::ParameterSet& iConfig):
  muons_          ( iConfig.getParameter<edm::InputTag>("muons") ),
  electrons_      ( iConfig.getParameter<edm::InputTag>("electrons") ),
  minMuPt_        ( iConfig.getParameter<double>("minMuPt") ),
  minElePt_       ( iConfig.getParameter<double>("minElePt") ),
  minLeadingLepPt_( iConfig.getParameter<double>("minLeadingLepPt") ),
  minNLeptons_    ( iConfig.getParameter<int>("minNLeptons") ),
  prescale_       ( iConfig.getParameter<int>("prescale") )
{
  std::string label = iConfig.getParameter<std::string>( "@module_label" );
  produces<bool>("passed").setBranchAlias("filter_"+label+"_passed");
  produces<bool>("run").setBranchAlias("filter_"+label+"_run");
  produces<int>("prescale").setBranchAlias("filter_"+label+"_prescale");
}

bool
ZZFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  std::vector<const reco::RecoCandidate*> cands;
  
  edm::Handle<std::vector<reco::Muon> > muons;
  iEvent.getByLabel(muons_, muons);
  double leadingLeptonPt(0);
  for ( std::vector<reco::Muon>::const_iterator muon = muons->begin(); muon != muons->end(); ++muon )
    {
      if ( muon->pt()>minMuPt_ ) cands.push_back(&*muon);
      if ( muon->pt()>leadingLeptonPt ) leadingLeptonPt = muon->pt();
    }  

  edm::Handle<std::vector<reco::GsfElectron> > eles;
  iEvent.getByLabel(electrons_, eles);
  for ( std::vector<reco::GsfElectron>::const_iterator ele = eles->begin(); ele != eles->end(); ++ele )
    {
      if ( ele->pt()>minElePt_ ) cands.push_back(&*ele);
      if ( ele->pt()>leadingLeptonPt ) leadingLeptonPt = ele->pt();
    }
  bool  passed = (cands.size()>=minNLeptons_) && (leadingLeptonPt >= minLeadingLepPt_);
  *filter_passed = passed;
  iEvent.put(filter_prescale, "prescale");
  iEvent.put(filter_passed,   "passed");
  iEvent.put(filter_run,      "run");
  return passed;
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZZFilter);
