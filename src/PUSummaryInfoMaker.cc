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
  
  produces<int>			(aliasprefix_ + "nPUvertices"	).setBranchAlias(aliasprefix_ + "_nPUvertices"	);
  produces<vector<float> >	(aliasprefix_ + "zpositions"	).setBranchAlias(aliasprefix_ + "_zpositions"	);
  produces<vector<float> >	(aliasprefix_ + "sumptlowpt"	).setBranchAlias(aliasprefix_ + "_sumpt_lowpt"	);
  produces<vector<float> >	(aliasprefix_ + "sumpthighpt"	).setBranchAlias(aliasprefix_+ "_sump_highpt"	);
  produces<vector<float> >	(aliasprefix_ + "instLumi"	).setBranchAlias(aliasprefix_+ "_instLumi"      );
  produces<vector<int>   >	(aliasprefix_ + "ntrkslowpt"	).setBranchAlias(aliasprefix_+ "_ntrks_lowpt"	);
  produces<vector<int>   >	(aliasprefix_ + "ntrkshighpt"	).setBranchAlias(aliasprefix_+ "_ntrks_highpt"	);

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
  

  auto_ptr<int>			puInfo_nPUvertices	(new int		);
  auto_ptr<vector<float> >	puInfo_zpositions	(new vector<float>	);
  auto_ptr<vector<float> >	puInfo_sumptlowpt	(new vector<float>	);
  auto_ptr<vector<float> >	puInfo_sumpthighpt	(new vector<float>	);
  auto_ptr<vector<float> >	puInfo_instLumi         (new vector<float>	);
  auto_ptr<vector<int>   >	puInfo_ntrkslowpt	(new vector<int>	);
  auto_ptr<vector<int>   >	puInfo_ntrkshighpt	(new vector<int>	);
  

  Handle<PileupSummaryInfo> puInfoH;
  bool bPuInfo=iEvent.getByLabel("addPileupInfo", puInfoH);

  if(bPuInfo) {
    *puInfo_nPUvertices = puInfoH->getPU_NumInteractions();
    *puInfo_zpositions  = puInfoH->getPU_zpositions();
    *puInfo_sumptlowpt  = puInfoH->getPU_sumpT_lowpT();
    *puInfo_sumpthighpt = puInfoH->getPU_sumpT_highpT();
    *puInfo_ntrkslowpt  = puInfoH->getPU_ntrks_lowpT();
    *puInfo_ntrkshighpt = puInfoH->getPU_ntrks_highpT();
    
  } else {
    *puInfo_nPUvertices = -9999;
    *puInfo_zpositions  = vector<float>();
    *puInfo_sumptlowpt  = vector<float>();
    *puInfo_sumpthighpt = vector<float>();
    *puInfo_ntrkslowpt  = vector<int>();
    *puInfo_ntrkshighpt = vector<int>();
  }

	
  iEvent.put(puInfo_nPUvertices,	"puInfonPUvertices"	);
  iEvent.put(puInfo_zpositions,		"puInfozpositions"	);
  iEvent.put(puInfo_sumptlowpt,		"puInfosumptlowpt"	);
  iEvent.put(puInfo_sumpthighpt,	"puInfosumpthighpt"	);
  iEvent.put(puInfo_instLumi,		"puInfoinstLumi"	);
  iEvent.put(puInfo_ntrkslowpt,		"puInfontrkslowpt"	);
  iEvent.put(puInfo_ntrkshighpt,	"puInfontrkshighpt"	);

  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PUSummaryInfoMaker);
