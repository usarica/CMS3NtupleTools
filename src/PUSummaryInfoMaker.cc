//-*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "CMS3/NtupleMaker/interface/PUSummaryInfoMaker.h" 

using namespace edm;
using namespace std;

//
// constructors and destructor
//

PUSummaryInfoMaker::PUSummaryInfoMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_		= iConfig.getUntrackedParameter	<std::string>	("aliasPrefix"		);
  PUInfoInputTag_	= iConfig.getParameter		<edm::InputTag>	("PUInfoInputTag"	);
  
  produces<vector<int> >	        (aliasprefix_ + "nPUvertices"	     ).setBranchAlias(aliasprefix_ + "_nPUvertices"	       );
  produces<vector<int> >	        (aliasprefix_ + "bunchCrossing"	     ).setBranchAlias(aliasprefix_ + "_bunchCrossing"      );
  produces<vector<vector<float> > >	(aliasprefix_ + "zpositions"	     ).setBranchAlias(aliasprefix_ + "_zpositions"	       );
  produces<vector<vector<float> > >	(aliasprefix_ + "sumptlowpt"	     ).setBranchAlias(aliasprefix_ + "_sumpt_lowpt"	       );
  produces<vector<vector<float> > >	(aliasprefix_ + "sumpthighpt"	     ).setBranchAlias(aliasprefix_ + "_sump_highpt"	       );
  produces<vector<vector<float> > >	(aliasprefix_ + "instLumi"	         ).setBranchAlias(aliasprefix_ + "_instLumi"           );
  produces<vector<vector<int> > >	(aliasprefix_ + "ntrkslowpt"	     ).setBranchAlias(aliasprefix_ + "_ntrks_lowpt"	       );
  produces<vector<vector<int> > >	(aliasprefix_ + "ntrkshighpt"	     ).setBranchAlias(aliasprefix_ + "_ntrks_highpt"	   );
  produces<vector<float> >	        (aliasprefix_ + "trueNumInteractions").setBranchAlias(aliasprefix_ + "_trueNumInteractions");

}

PUSummaryInfoMaker::~PUSummaryInfoMaker()
{
}

void  PUSummaryInfoMaker::beginJob()
{
}

void PUSummaryInfoMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void PUSummaryInfoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<int> >	        puInfo_nPUvertices	     ( new vector<int> );
  auto_ptr<vector<int> >	        puInfo_bunchCrossing	 ( new vector<int> );
  auto_ptr<vector<vector<float> > >	puInfo_zpositions	     ( new vector<vector<float> > );
  auto_ptr<vector<vector<float> > > puInfo_sumptlowpt	     ( new vector<vector<float> > );
  auto_ptr<vector<vector<float> > >	puInfo_sumpthighpt	     ( new vector<vector<float> > );
  auto_ptr<vector<vector<float> > >	puInfo_instLumi          ( new vector<vector<float> > );
  auto_ptr<vector<vector<int> > >	puInfo_ntrkslowpt	     ( new vector<vector<int> > );
  auto_ptr<vector<vector<int> > >	puInfo_ntrkshighpt	     ( new vector<vector<int> > );
  auto_ptr<vector<float> >	        puInfo_trueninteractions ( new vector<float> );

  Handle<vector<PileupSummaryInfo> > puInfoH;
  bool bPuInfo=iEvent.getByLabel(PUInfoInputTag_, puInfoH); 

  if(bPuInfo) {
	for ( vector<PileupSummaryInfo>::const_iterator itr = puInfoH->begin(); 
	  itr != puInfoH->end(); ++itr )
      {
		puInfo_nPUvertices->      push_back( itr->getPU_NumInteractions()  );
		puInfo_bunchCrossing->    push_back( itr->getBunchCrossing()       );
		puInfo_zpositions->       push_back( itr->getPU_zpositions()       );
		puInfo_sumptlowpt->       push_back( itr->getPU_sumpT_lowpT()      );
		puInfo_sumpthighpt->      push_back( itr->getPU_sumpT_highpT()     );
		puInfo_ntrkslowpt->       push_back( itr->getPU_ntrks_lowpT()      );
		puInfo_ntrkshighpt->      push_back( itr->getPU_ntrks_highpT()     );
		puInfo_trueninteractions->push_back( itr->getTrueNumInteractions() );
      }
  }
	
  iEvent.put(puInfo_nPUvertices,	    aliasprefix_ + "nPUvertices"	     );
  iEvent.put(puInfo_bunchCrossing,	    aliasprefix_ + "bunchCrossing"       );
  iEvent.put(puInfo_zpositions,		    aliasprefix_ + "zpositions"	         );
  iEvent.put(puInfo_sumptlowpt,		    aliasprefix_ + "sumptlowpt"	         );
  iEvent.put(puInfo_sumpthighpt,	    aliasprefix_ + "sumpthighpt"	     );
  iEvent.put(puInfo_instLumi,		    aliasprefix_ + "instLumi"	         );
  iEvent.put(puInfo_ntrkslowpt,		    aliasprefix_ + "ntrkslowpt"	         );
  iEvent.put(puInfo_ntrkshighpt,	    aliasprefix_ + "ntrkshighpt"	     );
  iEvent.put(puInfo_trueninteractions,	aliasprefix_ + "trueNumInteractions" );

  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PUSummaryInfoMaker);
