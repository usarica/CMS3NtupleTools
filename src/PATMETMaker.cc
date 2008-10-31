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
// $Id: PATMETMaker.cc,v 1.1 2008/10/31 23:14:18 kalavase Exp $
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
  produces<float> ("metpatmetCor").setBranchAlias("met_pat_metCor"); // fully corrected MET
  produces<float> ("metpatmetPhiCor").setBranchAlias("met_pat_metPhiCor"); // fully corrected METPhi
  produces<float> ("metpatmetUncor").setBranchAlias("met_pat_metUncor"); // MET with no corrections 
  produces<float> ("metpatmetPhiUncor").setBranchAlias("met_pat_metPhiUncor"); // MET with no corrections 
  produces<float> ("metpatmetUncorJES").setBranchAlias("met_pat_metUncorJES"); // PAT gen parton ID
  produces<float> ("metpatmetPhiUncorJES").setBranchAlias("met_pat_metPhiUncorJES"); // PAT gen parton ID
  produces<float> ("metpatmetUncorMuon").setBranchAlias("met_pat_metUncorMuon"); // PAT gen parton ID
  produces<float> ("metpatmetPhiUncorMuon").setBranchAlias("met_pat_metPhiUncorMuon"); // PAT gen parton ID
  
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
  
  // get jet collection
  edm::Handle<pat::MET> patMETHandle;
  iEvent.getByLabel(patMETInputTag, patMETHandle);
  
  // create containers
  auto_ptr<float> met_pat_metCor(new float);   //MET with all corrections
  auto_ptr<float> met_pat_metPhiCor(new float);   //MET with all corrections
  auto_ptr<float> met_pat_metUncor(new float); // MET with no corrections 
  auto_ptr<float> met_pat_metPhiUncor(new float); // MET with no corrections 
  auto_ptr<float> met_pat_metUncorJES(new float); //MET with no JES corrections
  auto_ptr<float> met_pat_metPhiUncorJES(new float); //MET with no JES corrections
  auto_ptr<float> met_pat_metUncorMuon(new float); //MET with no muon correction
  auto_ptr<float> met_pat_metPhiUncorMuon(new float); //MET with no muon correction
  
  *met_pat_metCor    = patMETHandle->pt();
  *met_pat_metPhiCor = patMETHandle->phi();
  
  *met_pat_metUncor    = patMETHandle->corSumEt(pat::MET::uncorrALL);
  *met_pat_metPhiUncor = (patMETHandle->corEx(pat::MET::uncorrALL)==0 && patMETHandle->corEy(pat::MET::uncorrALL)==0) 
    ? 0 : atan2(patMETHandle->corEy(pat::MET::uncorrALL), patMETHandle->corEx(pat::MET::uncorrALL) );
  
  *met_pat_metUncorJES    = patMETHandle->corSumEt(pat::MET::uncorrJES);
  *met_pat_metPhiUncorJES = (patMETHandle->corEx(pat::MET::uncorrJES)==0 && patMETHandle->corEy(pat::MET::uncorrJES)==0) 
    ? 0 : atan2(patMETHandle->corEy(pat::MET::uncorrJES), patMETHandle->corEx(pat::MET::uncorrJES) );
  
  *met_pat_metUncorMuon = patMETHandle->corSumEt(pat::MET::uncorrMUON);
  *met_pat_metPhiUncorMuon = (patMETHandle->corEx(pat::MET::uncorrMUON)==0 && patMETHandle->corEy(pat::MET::uncorrMUON)==0) 
    ? 0 : atan2(patMETHandle->corEy(pat::MET::uncorrMUON), patMETHandle->corEx(pat::MET::uncorrMUON) );

 
  
  // put containers into event
  iEvent.put(met_pat_metCor,"met_pat_metCor");
  iEvent.put(met_pat_metPhiCor,"met_pat_metPhiCor");
  iEvent.put(met_pat_metUncor, "met_pat_metUncor");
  iEvent.put(met_pat_metPhiUncor, "met_pat_metPhiUncor");
  iEvent.put(met_pat_metUncorJES,"met_pat_metUncorJES");
  iEvent.put(met_pat_metPhiUncorJES,"met_pat_metPhiUncorJES");
  iEvent.put(met_pat_metUncorMuon,"met_pat_metUncorMuon");
  iEvent.put(met_pat_metPhiUncorMuon,"met_pat_metPhiUncorMuon");


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
