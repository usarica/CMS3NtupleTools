// -*- C++ -*-
//
// Filter to select two events for two different cases:
// 1) dilepton hypothesis
// 2) dilepton + HT (sumJetPt) - SUSY
//
// The first one nominally has tight lepton pt requirements
//

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
typedef math::XYZTLorentzVectorF LorentzVector;
class HypDilepPruner : public edm::EDFilter {
public:
  explicit HypDilepPruner(const edm::ParameterSet&);
  ~HypDilepPruner(){}
  
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  edm::InputTag hyp_lt_p4_;
  edm::InputTag hyp_ll_p4_;
  edm::InputTag hyp_sumJetPt_;
  double minNominalLTPt_;
  double minNominalLLPt_;
  double minSusyLTPt_;
  double minSusyLLPt_;
  double minSusySumJetPt_;
  
};
HypDilepPruner::HypDilepPruner(const edm::ParameterSet& iConfig):
  hyp_lt_p4_      ( iConfig.getParameter<edm::InputTag>("hyp_lt_p4") ),
  hyp_ll_p4_      ( iConfig.getParameter<edm::InputTag>("hyp_ll_p4") ),
  hyp_sumJetPt_   ( iConfig.getParameter<edm::InputTag>("hyp_sumJetPt") ),
  minNominalLTPt_ ( iConfig.getParameter<double>("minNominalLTPt") ),
  minNominalLLPt_ ( iConfig.getParameter<double>("minNominalLLPt") ),
  minSusyLTPt_    ( iConfig.getParameter<double>("minSusyLTPt") ),
  minSusyLLPt_    ( iConfig.getParameter<double>("minSusyLLPt") ),
  minSusySumJetPt_( iConfig.getParameter<double>("minSusySumJetPt") )
{
}
bool
HypDilepPruner::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<LorentzVector> > hyp_lt_p4;
  iEvent.getByLabel(hyp_lt_p4_, hyp_lt_p4);
  edm::Handle<std::vector<LorentzVector> > hyp_ll_p4;
  iEvent.getByLabel(hyp_ll_p4_, hyp_ll_p4);
  edm::Handle<std::vector<float> > hyp_sumJetPt;
  iEvent.getByLabel(hyp_sumJetPt_, hyp_sumJetPt);
  for ( unsigned int i=0; i<hyp_lt_p4->size(); ++i )
    {
      if (std::max(hyp_lt_p4->at(i).pt(),hyp_ll_p4->at(i).pt()) >= minNominalLTPt_ &&
	  std::min(hyp_lt_p4->at(i).pt(),hyp_ll_p4->at(i).pt()) >= minNominalLLPt_ )
	return true;
      if (std::max(hyp_lt_p4->at(i).pt(),hyp_ll_p4->at(i).pt()) >= minSusyLTPt_ &&
	  std::min(hyp_lt_p4->at(i).pt(),hyp_ll_p4->at(i).pt()) >= minSusyLLPt_ &&
	  hyp_sumJetPt->at(i) > minSusySumJetPt_)
	return true;
    }
  return false;
}
//define this as a plug-in
DEFINE_FWK_MODULE(HypDilepPruner);
