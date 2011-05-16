//-*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "CMS2/NtupleMaker/interface/PUSummaryInfoMaker.h" 

using namespace edm;
using namespace std;

//
// constructors and destructor
//

PUSummaryInfoMaker::PUSummaryInfoMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_		= iConfig.getUntrackedParameter	<std::string>	("aliasPrefix"		);
  PUInfoInputTag_	= iConfig.getParameter		<edm::InputTag>	("PUInfoInputTag"	);
  
  produces<vector<int> >	        (aliasprefix_ + "nPUvertices"	).setBranchAlias(aliasprefix_ + "_nPUvertices"	);
  produces<vector<int> >	        (aliasprefix_ + "bunchCrossing"	).setBranchAlias(aliasprefix_ + "_bunchCrossing");
  produces<vector<vector<float> > >	(aliasprefix_ + "zpositions"	).setBranchAlias(aliasprefix_ + "_zpositions"	);
  produces<vector<vector<float> > >	(aliasprefix_ + "sumptlowpt"	).setBranchAlias(aliasprefix_ + "_sumpt_lowpt"	);
  produces<vector<vector<float> > >	(aliasprefix_ + "sumpthighpt"	).setBranchAlias(aliasprefix_ + "_sump_highpt"	);
  produces<vector<vector<float> > >	(aliasprefix_ + "instLumi"	).setBranchAlias(aliasprefix_ + "_instLumi"     );
  produces<vector<vector<int> > >	(aliasprefix_ + "ntrkslowpt"	).setBranchAlias(aliasprefix_ + "_ntrks_lowpt"	);
  produces<vector<vector<int> > >	(aliasprefix_ + "ntrkshighpt"	).setBranchAlias(aliasprefix_ + "_ntrks_highpt"	);

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
  

  auto_ptr<vector<int> >	        puInfo_nPUvertices	( new vector<int> );
  auto_ptr<vector<int> >	        puInfo_bunchCrossing	( new vector<int> );
  auto_ptr<vector<vector<float> > >	puInfo_zpositions	( new vector<vector<float> > );
  auto_ptr<vector<vector<float> > >     puInfo_sumptlowpt	( new vector<vector<float> > );
  auto_ptr<vector<vector<float> > >	puInfo_sumpthighpt	( new vector<vector<float> > );
  auto_ptr<vector<vector<float> > >	puInfo_instLumi         ( new vector<vector<float> > );
  auto_ptr<vector<vector<int> > >	puInfo_ntrkslowpt	( new vector<vector<int> > );
  auto_ptr<vector<vector<int> > >	puInfo_ntrkshighpt	( new vector<vector<int> > );

  Handle<vector<PileupSummaryInfo> > puInfoH;
  bool bPuInfo=iEvent.getByLabel("addPileupInfo", puInfoH);

  if(bPuInfo) {
    for ( vector<PileupSummaryInfo>::const_iterator itr = puInfoH->begin(); 
	  itr != puInfoH->end(); ++itr )
      {
	puInfo_nPUvertices->push_back(   itr->getPU_NumInteractions() );
	puInfo_bunchCrossing->push_back( itr->getBunchCrossing()      );
	puInfo_zpositions->push_back(    itr->getPU_zpositions()      );
	puInfo_sumptlowpt->push_back(    itr->getPU_sumpT_lowpT()     );
	puInfo_sumpthighpt->push_back(   itr->getPU_sumpT_highpT()    );
	puInfo_ntrkslowpt->push_back(    itr->getPU_ntrks_lowpT()     );
	puInfo_ntrkshighpt->push_back(   itr->getPU_ntrks_highpT()    );
      }
  }
	
  iEvent.put(puInfo_nPUvertices,	"puInfonPUvertices"	);
  iEvent.put(puInfo_bunchCrossing,	"puInfobunchCrossing"	);
  iEvent.put(puInfo_zpositions,		"puInfozpositions"	);
  iEvent.put(puInfo_sumptlowpt,		"puInfosumptlowpt"	);
  iEvent.put(puInfo_sumpthighpt,	"puInfosumpthighpt"	);
  iEvent.put(puInfo_instLumi,		"puInfoinstLumi"	);
  iEvent.put(puInfo_ntrkslowpt,		"puInfontrkslowpt"	);
  iEvent.put(puInfo_ntrkshighpt,	"puInfontrkshighpt"	);

  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PUSummaryInfoMaker);
