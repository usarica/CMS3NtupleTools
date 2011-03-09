//-*- C++ -*-
//
// Package:    FastJetMaker
// Class:      FastJetMaker
// 
/**\class FastJetMaker FastJetMaker.cc CMS2/FastJetMaker/src/FastJetMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: FastJetMaker.cc,v 1.1 2011/03/09 16:49:14 benhoob Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/FastJetMaker.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/Common/interface/ValueMap.h" 

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

FastJetMaker::FastJetMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
       branchprefix.replace(branchprefix.find("_"),1,"");

  produces<float>         (branchprefix+"rho"       ).setBranchAlias(aliasprefix_+"_rho"       );

  // input tags
  rho_tag     = iConfig.getParameter<edm::InputTag>("rho_tag");
}


FastJetMaker::~FastJetMaker()
{
}

void  FastJetMaker::beginJob()
{
}

void FastJetMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void FastJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float> evt_rho(new float);
 
  edm::Handle<double> rhoH;
  iEvent.getByLabel( rho_tag , rhoH);
  
  *evt_rho = *rhoH; 

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
       branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_rho            , branchprefix+"rho"           );

}

//define this as a plug-in
DEFINE_FWK_MODULE(FastJetMaker);





  
