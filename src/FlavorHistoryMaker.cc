//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      FlavorHistoryMaker
// 
/**\class FlavorHistoryMaker FlavorHistoryMaker.cc CMS2/NtupleMaker/src/FlavorHistoryMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Tues Sep  1 11:07:38 CDT 2009
// $Id: FlavorHistoryMaker.cc,v 1.3 2010/03/18 02:12:11 kalavase Exp $
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/FlavorHistoryMaker.h" 


//
// constructors and destructor
//

FlavorHistoryMaker::FlavorHistoryMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
  
  using namespace edm;
  using namespace std;

  produces<unsigned int>    (branchprefix+"flavorHistoryFilterResult"        ).setBranchAlias(aliasprefix_+"_flavorHistoryFilterResult");
  
  flavorHistoryFilterTag_ = iConfig.getParameter<InputTag>    ("flavorHistoryFilterTag"                  );
  
}

FlavorHistoryMaker::~FlavorHistoryMaker()
{
}

void  FlavorHistoryMaker::beginJob()
{
}

void FlavorHistoryMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void FlavorHistoryMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  
  auto_ptr<unsigned int>    genps_flavorHistoryFilterResult     (new unsigned int);

  //get the flavor history result
  edm::Handle<unsigned int> flavorHistoryFilterH;
  iEvent.getByLabel(flavorHistoryFilterTag_, flavorHistoryFilterH);

  if( !flavorHistoryFilterH.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve gen particle collection";
    edm::LogInfo("OutputInfo") << " FlavorHistoryMaker cannot continue...!";
    throw cms::Exception("FlavorHistoryMaker: FlavorHistoryFilter Handle is not valid. Exiting");
  }
  

  *genps_flavorHistoryFilterResult = *flavorHistoryFilterH.product();
  

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(genps_flavorHistoryFilterResult             , branchprefix+"flavorHistoryFilterResult");

}

//define this as a plug-in
DEFINE_FWK_MODULE(FlavorHistoryMaker);
