// -*- C++ -*-
//
// Package:    PATMETMaker
// Class:      PATMETMaker
// 
/**\class PATMETMaker PATMETMaker.cc CMS2/NtupleMaker/src/PATMETMaker.cc

Description: copy additional PAT met variables in simple data structures into the EDM event tree

 Implementation:
     - take PAT mets
     - extract and fill variables
*/
//
// Original Author:  pts/4
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATMETMaker.cc,v 1.3 2009/05/18 22:43:06 fgolf Exp $
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/PATMETMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
//
// class decleration
//

//
// constructors and destructor
//
PATMETMaker::PATMETMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<float> ("metpatmetCor"         ).setBranchAlias("met_pat_metCor"         ); // MET corrected for muons+JES
  produces<float> ("metpatmetPhiCor"      ).setBranchAlias("met_pat_metPhiCor"      ); // phi of MET corrected for muons+JES
  produces<float> ("metpatmetUncor"       ).setBranchAlias("met_pat_metUncor"       ); // MET with no corrections 
  produces<float> ("metpatmetPhiUncor"    ).setBranchAlias("met_pat_metPhiUncor"    ); // phi of MET with no corrections 
  produces<float> ("metpatmetUncorJES"    ).setBranchAlias("met_pat_metUncorJES"    ); // MET corrected for muons only
  produces<float> ("metpatmetPhiUncorJES" ).setBranchAlias("met_pat_metPhiUncorJES" ); // phi of MET corrected for muons only
  produces<float> ("metpatmetUncorMuon"   ).setBranchAlias("met_pat_metUncorMuon"   ); // MET corrected for JES only
  produces<float> ("metpatmetPhiUncorMuon").setBranchAlias("met_pat_metPhiUncorMuon"); // phi of MET corrected for JES only
  
   // parameters from configuration
  patMETInputTag = iConfig.getParameter<edm::InputTag>("patMETsInputTag");
}


PATMETMaker::~PATMETMaker()
{
 
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
PATMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  edm::Handle<edm::View<pat::MET> > patMETView;
  iEvent.getByLabel(patMETInputTag, patMETView);
  const pat::MET met = (*patMETView).at(0);  
    
  // create containers
  auto_ptr<float> met_pat_metCor         (new float);
  auto_ptr<float> met_pat_metPhiCor      (new float);
  auto_ptr<float> met_pat_metUncor       (new float);
  auto_ptr<float> met_pat_metPhiUncor    (new float);
  auto_ptr<float> met_pat_metUncorJES    (new float);
  auto_ptr<float> met_pat_metPhiUncorJES (new float);
  auto_ptr<float> met_pat_metUncorMuon   (new float);
  auto_ptr<float> met_pat_metPhiUncorMuon(new float);


  *met_pat_metCor          = met.pt();
  *met_pat_metPhiCor       = met.phi();
  
  *met_pat_metUncor        = met.uncorrectedPt();
  *met_pat_metPhiUncor     = met.uncorrectedPhi();
  
  *met_pat_metUncorJES     = met.uncorrectedPt (pat::MET::uncorrJES);
  *met_pat_metPhiUncorJES  = met.uncorrectedPhi(pat::MET::uncorrJES);
  
  *met_pat_metUncorMuon    = met.uncorrectedPt (pat::MET::uncorrMUON);
  *met_pat_metPhiUncorMuon = met.uncorrectedPhi(pat::MET::uncorrMUON);
  
  
  // put containers into event
  iEvent.put(met_pat_metCor         , "metpatmetCor"          );
  iEvent.put(met_pat_metPhiCor      , "metpatmetPhiCor"       );
  iEvent.put(met_pat_metUncor       , "metpatmetUncor"        );
  iEvent.put(met_pat_metPhiUncor    , "metpatmetPhiUncor"     );
  iEvent.put(met_pat_metUncorJES    , "metpatmetUncorJES"     );
  iEvent.put(met_pat_metPhiUncorJES , "metpatmetPhiUncorJES"  );
  iEvent.put(met_pat_metUncorMuon   , "metpatmetUncorMuon"    );
  iEvent.put(met_pat_metPhiUncorMuon, "metpatmetPhiUncorMuon" );
}
// ------------ method called once each job just before starting event loop  ------------
void 
PATMETMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATMETMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATMETMaker);
