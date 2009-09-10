//-*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/HypGenMaker.h" 
//#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

HypGenMaker::HypGenMaker(const edm::ParameterSet& iConfig) {

  //muonsInputTag            = iConfig.getParameter<InputTag>("muonsInputTag"             );
  //electronsInputTag        = iConfig.getParameter<InputTag>("electronsInputTag"         );
  candToGenAssTag          = iConfig.getParameter<InputTag>("candToGenAssTag");  
  hypInputTag              = iConfig.getParameter<edm::InputTag>("hypInputTag");
  
  produces<vector<int> >           ("hypltmcid"                ).setBranchAlias("hyp_lt_mc_id"                 );
  produces<vector<int> >           ("hypltmcmotherid"          ).setBranchAlias("hyp_lt_mc_motherid"           );
  produces<vector<LorentzVector > >("hypltmcp4"                ).setBranchAlias("hyp_lt_mc_p4"                 );
  
  produces<vector<int> >           ("hypllmcid"                ).setBranchAlias("hyp_ll_mc_id"                 );
  produces<vector<int> >           ("hypllmcmotherid"          ).setBranchAlias("hyp_ll_mc_motherid"           );
  produces<vector<LorentzVector > >("hypllmcp4"                ).setBranchAlias("hyp_ll_mc_p4"                 );
  

}

HypGenMaker::~HypGenMaker()
{
}

void  HypGenMaker::beginJob(const edm::EventSetup&)
{
}

void HypGenMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void HypGenMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<int> >   hyp_lt_mc_id                 (new vector<int>);
  auto_ptr<vector<int> >   hyp_lt_mc_motherid           (new vector<int>);
  auto_ptr<vector<LorentzVector> > hyp_lt_mc_p4         (new vector<LorentzVector>);
  
  auto_ptr<vector<int> >   hyp_ll_mc_id                 (new vector<int>);
  auto_ptr<vector<int> >   hyp_ll_mc_motherid           (new vector<int>);
  auto_ptr<vector<LorentzVector> > hyp_ll_mc_p4         (new vector<LorentzVector>);


  // get MC particle collection
  //edm::Handle<reco::GenParticleCollection> genpsHandle;
  //iEvent.getByLabel(genParticlesInputTag, genpsHandle);

  //if( !genpsHandle.isValid() ) {
  //  edm::LogInfo("OutputInfo") << " failed to retrieve gen particle collection";
  //  edm::LogInfo("OutputInfo") << " HypGenMaker cannot continue...!";
  //  return;
  //}

  //PDG id of matched MC particle
  InputTag mus_mc_id_tag(candToGenAssTag.label(),"musmcid");
  Handle<vector<int> > mus_mc_id_h;
  iEvent.getByLabel(mus_mc_id_tag, mus_mc_id_h);
  const vector<int> *mus_mc_id = mus_mc_id_h.product();
  
  //PDG id of MC matched mother 
  InputTag mus_mc_motherid_tag(candToGenAssTag.label(),"musmcmotherid");
  Handle<vector<int> > mus_mc_motherid_h;
  iEvent.getByLabel(mus_mc_motherid_tag, mus_mc_motherid_h);
  const vector<int> *mus_mc_motherid = mus_mc_motherid_h.product();
  
  //muon mc P4
  InputTag mus_mc_p4_tag(candToGenAssTag.label(),"musmcp4");
  Handle<vector<LorentzVector> > mus_mc_p4_h;
  iEvent.getByLabel(mus_mc_p4_tag, mus_mc_p4_h);
  const vector<LorentzVector> *mus_mc_p4 = mus_mc_p4_h.product();
  
  //PDG id of matched MC particle
  InputTag els_mc_id_tag(candToGenAssTag.label(),"elsmcid");
  Handle<vector<int> > els_mc_id_h;
  iEvent.getByLabel(els_mc_id_tag, els_mc_id_h);
  const vector<int> *els_mc_id = els_mc_id_h.product();
  
  //PDG id of MC matched mother 
  InputTag els_mc_motherid_tag(candToGenAssTag.label(),"elsmcmotherid");
  Handle<vector<int> > els_mc_motherid_h;
  iEvent.getByLabel(els_mc_motherid_tag, els_mc_motherid_h);
  const vector<int> *els_mc_motherid = els_mc_motherid_h.product();
  
  //electron mc P4
  InputTag els_mc_p4_tag(candToGenAssTag.label(),"elsmcp4");
  Handle<vector<LorentzVector> > els_mc_p4_h;
  iEvent.getByLabel(els_mc_p4_tag, els_mc_p4_h);
  const vector<LorentzVector> *els_mc_p4 = els_mc_p4_h.product();

  // hyp_type
  InputTag hyp_type_tag(hypInputTag.label(),"hyptype");
  Handle<vector<int> > hyp_type_h;
  iEvent.getByLabel(hyp_type_tag, hyp_type_h);
  const vector<int> *hyp_type = hyp_type_h.product();

  //hyp_lt_index
  edm::InputTag hyp_lt_index_tag(hypInputTag.label(),"hypltindex");
  edm::Handle<std::vector<int> > hyp_lt_index_h;
  iEvent.getByLabel(hyp_lt_index_tag, hyp_lt_index_h);
  const vector<int> *hyp_lt_index = hyp_lt_index_h.product();

  //hyp_ll_index
  edm::InputTag hyp_ll_index_tag(hypInputTag.label(),"hypllindex");
  edm::Handle<std::vector<int> > hyp_ll_index_h;
  iEvent.getByLabel(hyp_ll_index_tag, hyp_ll_index_h);
  const vector<int> *hyp_ll_index = hyp_ll_index_h.product();


  for( unsigned int i=0; i<hyp_type->size(); i++ ) {
	if( hyp_type->at(i) == 0 ) { //mumu
	  hyp_lt_mc_id        ->push_back(mus_mc_id        ->at(hyp_lt_index->at(i))  );
      hyp_lt_mc_motherid  ->push_back(mus_mc_motherid  ->at(hyp_lt_index->at(i))  );
      hyp_lt_mc_p4        ->push_back(mus_mc_p4        ->at(hyp_lt_index->at(i))  );
	  
      hyp_ll_mc_id        ->push_back(mus_mc_id        ->at(hyp_ll_index->at(i))  );
      hyp_ll_mc_motherid  ->push_back(mus_mc_motherid  ->at(hyp_ll_index->at(i))  );
      hyp_ll_mc_p4        ->push_back(mus_mc_p4        ->at(hyp_ll_index->at(i))  );
	}														
	else if( hyp_type->at(i) == 3 ) { //ee						
      hyp_lt_mc_id        ->push_back(els_mc_id        ->at(hyp_lt_index->at(i))  );
      hyp_lt_mc_motherid  ->push_back(els_mc_motherid  ->at(hyp_lt_index->at(i))  );
      hyp_lt_mc_p4        ->push_back(els_mc_p4        ->at(hyp_lt_index->at(i))  );
	  														
      hyp_ll_mc_id        ->push_back(els_mc_id        ->at(hyp_ll_index->at(i))  );
      hyp_ll_mc_motherid  ->push_back(els_mc_motherid  ->at(hyp_ll_index->at(i))  );
      hyp_ll_mc_p4        ->push_back(els_mc_p4        ->at(hyp_ll_index->at(i))  );
	}														
	else if( hyp_type->at(i) == 1 ) { //mu tight, e loose		
	  hyp_lt_mc_id        ->push_back(mus_mc_id        ->at(hyp_lt_index->at(i))  );
	  hyp_lt_mc_motherid  ->push_back(mus_mc_motherid  ->at(hyp_lt_index->at(i))  );
	  hyp_lt_mc_p4        ->push_back(mus_mc_p4        ->at(hyp_lt_index->at(i))  );
	  														
	  hyp_ll_mc_id        ->push_back(els_mc_id        ->at(hyp_ll_index->at(i))  );
	  hyp_ll_mc_motherid  ->push_back(els_mc_motherid  ->at(hyp_ll_index->at(i))  );
	  hyp_ll_mc_p4        ->push_back(els_mc_p4        ->at(hyp_ll_index->at(i))  );
	}														
	else if( hyp_type->at(i) == 2 ) { //e tight, mu loose		
	  hyp_lt_mc_id        ->push_back(els_mc_id        ->at(hyp_lt_index->at(i))  );
	  hyp_lt_mc_motherid  ->push_back(els_mc_motherid  ->at(hyp_lt_index->at(i))  );
	  hyp_lt_mc_p4        ->push_back(els_mc_p4        ->at(hyp_lt_index->at(i))  );
	  														
	  hyp_ll_mc_id        ->push_back(mus_mc_id        ->at(hyp_ll_index->at(i))  );
	  hyp_ll_mc_motherid  ->push_back(mus_mc_motherid  ->at(hyp_ll_index->at(i))  );
	  hyp_ll_mc_p4        ->push_back(mus_mc_p4        ->at(hyp_ll_index->at(i))  );
	}
	else
	  cout << "\nHypGenMaker :: Bad hyp_type \n\n";
  }
	
  iEvent.put(hyp_lt_mc_id                 ,"hypltmcid"                   );
  iEvent.put(hyp_lt_mc_motherid           ,"hypltmcmotherid"             );
  iEvent.put(hyp_lt_mc_p4                 ,"hypltmcp4"                   );
  
  iEvent.put(hyp_ll_mc_id                 ,"hypllmcid"                   );
  iEvent.put(hyp_ll_mc_motherid           ,"hypllmcmotherid"             );
  iEvent.put(hyp_ll_mc_p4                 ,"hypllmcp4"                   );

}

//define this as a plug-in
DEFINE_FWK_MODULE(HypGenMaker);
