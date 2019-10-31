//-*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "CMS3/NtupleMaker/interface/PUSummaryInfoMaker.h" 

using namespace edm;
using namespace std;

PUSummaryInfoMaker::PUSummaryInfoMaker(const edm::ParameterSet& iConfig){

  aliasprefix_		= iConfig.getUntrackedParameter	<std::string>	("aliasPrefix"		);
  PUInfoToken = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfoInputTag"));
  
  produces<vector<int> >	        (aliasprefix_ + "nPUvertices"	     ).setBranchAlias(aliasprefix_ + "_nPUvertices"	       );
  produces<vector<int> >	        (aliasprefix_ + "bunchCrossing"	     ).setBranchAlias(aliasprefix_ + "_bunchCrossing"      );
  produces<vector<vector<float> > >	(aliasprefix_ + "instLumi"	         ).setBranchAlias(aliasprefix_ + "_instLumi"           );
  produces<vector<float> >	        (aliasprefix_ + "trueNumInteractions").setBranchAlias(aliasprefix_ + "_trueNumInteractions");

}

PUSummaryInfoMaker::~PUSummaryInfoMaker(){
}

void  PUSummaryInfoMaker::beginJob(){
}

void PUSummaryInfoMaker::endJob(){
}

// ------------ method called to produce the data  ------------
void PUSummaryInfoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  unique_ptr<vector<int> >	        puInfo_nPUvertices	     ( new vector<int> );
  unique_ptr<vector<int> >	        puInfo_bunchCrossing	 ( new vector<int> );
  unique_ptr<vector<vector<float> > >	puInfo_instLumi          ( new vector<vector<float> > );
  unique_ptr<vector<float> >	        puInfo_trueninteractions ( new vector<float> );

  Handle<vector<PileupSummaryInfo> > puInfoH;
  bool bPuInfo=iEvent.getByToken(PUInfoToken, puInfoH); 

  if(bPuInfo){
	for (vector<PileupSummaryInfo>::const_iterator itr = puInfoH->begin(); itr != puInfoH->end(); ++itr){
      puInfo_nPUvertices->      push_back( itr->getPU_NumInteractions()  );
      puInfo_bunchCrossing->    push_back( itr->getBunchCrossing()       );
      puInfo_trueninteractions->push_back( itr->getTrueNumInteractions() );
    }
  }
	
  iEvent.put(std::move(puInfo_nPUvertices),	    aliasprefix_ + "nPUvertices"	     );
  iEvent.put(std::move(puInfo_bunchCrossing),	    aliasprefix_ + "bunchCrossing"       );
  iEvent.put(std::move(puInfo_instLumi),		    aliasprefix_ + "instLumi"	         );
  iEvent.put(std::move(puInfo_trueninteractions),	aliasprefix_ + "trueNumInteractions" );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PUSummaryInfoMaker);
