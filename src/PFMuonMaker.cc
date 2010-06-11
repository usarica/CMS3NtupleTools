//-*- C++ -*-
//
// Package:    PFMuonMaker
// Class:      PFMuonMaker
// 
/**\class PFMuonMaker PFMuonMaker.cc CMS2/PFMuonMaker/src/PFMuonMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMuonMaker.cc,v 1.1 2010/06/11 17:26:50 fgolf Exp $
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h" 

#include "CMS2/NtupleMaker/interface/PFMuonMaker.h"


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
PFMuonMaker::PFMuonMaker(const edm::ParameterSet& iConfig) {

     pfCandidatesTag_	= iConfig.getParameter<InputTag>("pfCandidatesTag");
     isoc_vm_tag_	= iConfig.getParameter<edm::InputTag>("isoc_vm_tag");
     ison_vm_tag_	= iConfig.getParameter<edm::InputTag>("ison_vm_tag");
     isop_vm_tag_	= iConfig.getParameter<edm::InputTag>("isop_vm_tag");
     pfAllMuons_tag_	= iConfig.getParameter<edm::InputTag>("pfAllMuons_tag");     

     produces<vector<LorentzVector>	> ("pfmusp4"                    ).setBranchAlias("pfmus_p4"			);
     produces<vector<LorentzVector>	> ("pfmusposAtEcalp4"		).setBranchAlias("pfmus_posAtEcal_p4"		);
     produces<vector<float>		> ("pfmusecalE"			).setBranchAlias("pfmus_ecalE"			);
     produces<vector<float>		> ("pfmushcalE"			).setBranchAlias("pfmus_hcalE"			);
     produces<vector<float>		> ("pfmusrawEcalE"		).setBranchAlias("pfmus_rawEcalE"		);
     produces<vector<float>		> ("pfmusrawHcalE"		).setBranchAlias("pfmus_rawHcalE"		);
     produces<vector<float>		> ("pfmuspS1E"			).setBranchAlias("pfmus_pS1E"			);
     produces<vector<float>		> ("pfmuspS2E"			).setBranchAlias("pfmus_pS2E"			);
     produces<vector<float>		> ("pfmusdeltaP"		).setBranchAlias("pfmus_deltaP"			);
     produces<vector<float>		> ("pfmusmvaepi"		).setBranchAlias("pfmus_mva_epi"		);
     produces<vector<float>		> ("pfmusmvaemu"		).setBranchAlias("pfmus_mva_emu"		);
     produces<vector<float>		> ("pfmusmvapimu"		).setBranchAlias("pfmus_mva_pimu"		);
     produces<vector<float>		> ("pfmusmvanothinggamma"	).setBranchAlias("pfmus_mva_nothing_gamma"	);
     produces<vector<float>		> ("pfmusmvanothingnh"		).setBranchAlias("pfmus_mva_nothing_nh"		);
     produces<vector<int>		> ("pfmuscharge"		).setBranchAlias("pfmus_charge"			);
     produces<vector<int>		> ("pfmusparticleId"		).setBranchAlias("pfmus_particleId"		);
     produces<vector<int>  		> ("pfmusflag"			).setBranchAlias("pfmus_flag"			);
     produces<vector<float>             > ("pfmusisoChargedHadrons"     ).setBranchAlias("pfmus_isoChargedHadrons"      );
     produces<vector<float>             > ("pfmusisoNeutralHadrons"     ).setBranchAlias("pfmus_isoNeutralHadrons"      );
     produces<vector<float>             > ("pfmusisoPhotons"            ).setBranchAlias("pfmus_isoPhotons"             );
}

PFMuonMaker::~PFMuonMaker() {

}


void  PFMuonMaker::beginRun(edm::Run&, const edm::EventSetup& es) {
     

}

void PFMuonMaker::beginJob() {


}

void PFMuonMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void PFMuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     auto_ptr<vector<LorentzVector> >	pfmus_p4		(new vector<LorentzVector>	);
     auto_ptr<vector<LorentzVector> >	pfmus_posAtEcal_p4	(new vector<LorentzVector>	);
     auto_ptr<vector<float> >		pfmus_ecalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_hcalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_rawEcalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_rawHcalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_pS1E		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_pS2E		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_deltaP		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_mva_epi		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_mva_emu		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_mva_pimu		(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_mva_nothing_gamma	(new vector<float>		);
     auto_ptr<vector<float> >		pfmus_mva_nothing_nh	(new vector<float>		);
     auto_ptr<vector<int> >		pfmus_charge		(new vector<int>		);
     auto_ptr<vector<int> >	        pfmus_particleId	(new vector<int>		);
     auto_ptr<vector<int> >	        pfmus_flag		(new vector<int>        	);
     auto_ptr<vector<float> >	        pfmus_isoChargedHadrons (new vector<float>		);
     auto_ptr<vector<float> >	        pfmus_isoNeutralHadrons (new vector<float>		);
     auto_ptr<vector<float> >	        pfmus_isoPhotons        (new vector<float>		);


     Handle<PFCandidateCollection> pfCandidatesHandle;
     iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
     const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

     edm::Handle<edm::ValueMap<double> > isocHandle;
     iEvent.getByLabel(isoc_vm_tag_, isocHandle);
     edm::ValueMap<double> isoc_data = *isocHandle;

     edm::Handle<edm::ValueMap<double> > isonHandle;
     iEvent.getByLabel(ison_vm_tag_, isonHandle);
     edm::ValueMap<double> ison_data = *isonHandle;

     edm::Handle<edm::ValueMap<double> > isopHandle;
     iEvent.getByLabel(isop_vm_tag_, isopHandle);
     edm::ValueMap<double> isop_data = *isopHandle;

     edm::Handle<reco::PFCandidateCollection> pfAllMuonsHandle;
     iEvent.getByLabel(pfAllMuons_tag_, pfAllMuonsHandle);

     for(PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++)
     {
	  if (pf_it->particleId() != PFCandidate::mu)
	       continue;

	  int pfflags = 0;
	  for(unsigned int i = 0; i < 17; i++) {
	       if(pf_it->flag((PFCandidate::Flags)i))
		    pfflags |= (1<<i);
	  }

	  pfmus_p4			->push_back(LorentzVector(pf_it->px(), pf_it->py(), pf_it->pz(), pf_it->p())	);
	  pfmus_posAtEcal_p4		->push_back(LorentzVector(pf_it->positionAtECALEntrance().x(),
							       pf_it->positionAtECALEntrance().y(),
							       pf_it->positionAtECALEntrance().z(),
							       0.0)							);
	  pfmus_ecalE			->push_back(pf_it->ecalEnergy()							);
	  pfmus_hcalE			->push_back(pf_it->hcalEnergy()							);
	  pfmus_rawEcalE		->push_back(pf_it->rawEcalEnergy()						);
	  pfmus_rawHcalE		->push_back(pf_it->rawHcalEnergy()						);
	  pfmus_pS1E			->push_back(pf_it->pS1Energy()							);
	  pfmus_pS2E			->push_back(pf_it->pS2Energy()							);
	  pfmus_deltaP			->push_back(pf_it->deltaP()							);
	  pfmus_mva_epi			->push_back(pf_it->mva_e_pi()							);
	  pfmus_mva_emu			->push_back(pf_it->mva_e_mu()							);
	  pfmus_mva_pimu		->push_back(pf_it->mva_pi_mu()							);
	  pfmus_mva_nothing_gamma	->push_back(pf_it->mva_nothing_gamma()						);
	  pfmus_mva_nothing_nh		->push_back(pf_it->mva_nothing_nh()						);
	  pfmus_charge			->push_back(pf_it->charge()							);
	  pfmus_particleId		->push_back(pf_it->translateTypeToPdgId(pf_it->particleId())			);
	  pfmus_flag			->push_back(pfflags                                                             );
    
     }//loop over candidate collection

     for(unsigned int  pfit = 0; pfit < pfAllMuonsHandle->size(); pfit++)
     {
	  reco::PFCandidateRef pfref(pfAllMuonsHandle, pfit);
	  float isoc = (isoc_data)[pfref];
	  float ison = (ison_data)[pfref];
	  float isop = (isop_data)[pfref];

	  pfmus_isoChargedHadrons	->push_back(isoc);
	  pfmus_isoNeutralHadrons	->push_back(ison);
	  pfmus_isoPhotons		->push_back(isop);
     }

     iEvent.put(pfmus_p4,		  	"pfmusp4"		);
     iEvent.put(pfmus_posAtEcal_p4,		"pfmusposAtEcalp4"	);
     iEvent.put(pfmus_ecalE,		 	"pfmusecalE"		);
     iEvent.put(pfmus_hcalE,		 	"pfmushcalE"		);
     iEvent.put(pfmus_rawEcalE,			"pfmusrawEcalE"		);
     iEvent.put(pfmus_rawHcalE,			"pfmusrawHcalE"		);
     iEvent.put(pfmus_pS1E,			"pfmuspS1E"		);
     iEvent.put(pfmus_pS2E,			"pfmuspS2E"		);
     iEvent.put(pfmus_deltaP,			"pfmusdeltaP"		);
     iEvent.put(pfmus_mva_epi,			"pfmusmvaepi"		);
     iEvent.put(pfmus_mva_emu,			"pfmusmvaemu"		);
     iEvent.put(pfmus_mva_pimu,			"pfmusmvapimu"		);
     iEvent.put(pfmus_mva_nothing_gamma,	"pfmusmvanothinggamma"	);
     iEvent.put(pfmus_mva_nothing_nh,		"pfmusmvanothingnh"	);
     iEvent.put(pfmus_charge,			"pfmuscharge"		);
     iEvent.put(pfmus_particleId,		"pfmusparticleId"	);
     iEvent.put(pfmus_flag,			"pfmusflag"		);
     iEvent.put(pfmus_isoChargedHadrons,        "pfmusisoChargedHadrons");
     iEvent.put(pfmus_isoNeutralHadrons,        "pfmusisoNeutralHadrons");
     iEvent.put(pfmus_isoPhotons,               "pfmusisoPhotons"       );
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMuonMaker);

