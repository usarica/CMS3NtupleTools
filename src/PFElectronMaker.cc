//-*- C++ -*-
//
// Package:    PFElectronMaker
// Class:      PFElectronMaker
// 
/**\class PFElectronMaker PFElectronMaker.cc CMS2/PFElectronMaker/src/PFElectronMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFElectronMaker.cc,v 1.2 2010/06/11 17:24:04 fgolf Exp $
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

#include "CMS2/NtupleMaker/interface/PFElectronMaker.h"


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
PFElectronMaker::PFElectronMaker(const edm::ParameterSet& iConfig) {

     pfCandidatesTag_		= iConfig.getParameter<InputTag>	("pfCandidatesTag"	);
     isoc_vm_tag_		= iConfig.getParameter<edm::InputTag>	("isoc_vm_tag"		);
     ison_vm_tag_		= iConfig.getParameter<edm::InputTag>	("ison_vm_tag"		);
     isop_vm_tag_		= iConfig.getParameter<edm::InputTag>	("isop_vm_tag"		);
     pfAllElectrons_tag_	= iConfig.getParameter<edm::InputTag>	("pfAllElectrons_tag"	);     

     produces<vector<LorentzVector>	> ("pfelsp4"                    ).setBranchAlias("pfels_p4"			);
     produces<vector<LorentzVector>	> ("pfelsposAtEcalp4"		).setBranchAlias("pfels_posAtEcal_p4"		);
     produces<vector<float>		> ("pfelsecalE"			).setBranchAlias("pfels_ecalE"			);
     produces<vector<float>		> ("pfelshcalE"			).setBranchAlias("pfels_hcalE"			);
     produces<vector<float>		> ("pfelsrawEcalE"		).setBranchAlias("pfels_rawEcalE"		);
     produces<vector<float>		> ("pfelsrawHcalE"		).setBranchAlias("pfels_rawHcalE"		);
     produces<vector<float>		> ("pfelspS1E"			).setBranchAlias("pfels_pS1E"			);
     produces<vector<float>		> ("pfelspS2E"			).setBranchAlias("pfels_pS2E"			);
     produces<vector<float>		> ("pfelsdeltaP"		).setBranchAlias("pfels_deltaP"			);
     produces<vector<float>		> ("pfelsmvaepi"		).setBranchAlias("pfels_mva_epi"		);
     produces<vector<float>		> ("pfelsmvaemu"		).setBranchAlias("pfels_mva_emu"		);
     produces<vector<float>		> ("pfelsmvapimu"		).setBranchAlias("pfels_mva_pimu"		);
     produces<vector<float>		> ("pfelsmvanothinggamma"	).setBranchAlias("pfels_mva_nothing_gamma"	);
     produces<vector<float>		> ("pfelsmvanothingnh"		).setBranchAlias("pfels_mva_nothing_nh"		);
     produces<vector<int>		> ("pfelscharge"		).setBranchAlias("pfels_charge"			);
     produces<vector<int>		> ("pfelsparticleId"		).setBranchAlias("pfels_particleId"		);
     produces<vector<int>  		> ("pfelsflag"			).setBranchAlias("pfels_flag"			);
     produces<vector<float>             > ("pfelsisoChargedHadrons"     ).setBranchAlias("pfels_isoChargedHadrons"      );
     produces<vector<float>             > ("pfelsisoNeutralHadrons"     ).setBranchAlias("pfels_isoNeutralHadrons"      );
     produces<vector<float>             > ("pfelsisoPhotons"            ).setBranchAlias("pfels_isoPhotons"             );
}

PFElectronMaker::~PFElectronMaker() {

}


void  PFElectronMaker::beginRun(edm::Run&, const edm::EventSetup& es) {
     

}

void PFElectronMaker::beginJob() {


}

void PFElectronMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void PFElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     auto_ptr<vector<LorentzVector> >	pfels_p4		(new vector<LorentzVector>	);
     auto_ptr<vector<LorentzVector> >	pfels_posAtEcal_p4	(new vector<LorentzVector>	);
     auto_ptr<vector<float> >		pfels_ecalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_hcalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_rawEcalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_rawHcalE		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_pS1E		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_pS2E		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_deltaP		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_mva_epi		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_mva_emu		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_mva_pimu		(new vector<float>		);
     auto_ptr<vector<float> >		pfels_mva_nothing_gamma	(new vector<float>		);
     auto_ptr<vector<float> >		pfels_mva_nothing_nh	(new vector<float>		);
     auto_ptr<vector<int> >		pfels_charge		(new vector<int>		);
     auto_ptr<vector<int> >		pfels_particleId	(new vector<int>		);
     auto_ptr<vector<int> >	        pfels_flag		(new vector<int>        	);
     auto_ptr<vector<float> >	        pfels_isoChargedHadrons (new vector<float>		);
     auto_ptr<vector<float> >	        pfels_isoNeutralHadrons (new vector<float>		);
     auto_ptr<vector<float> >	        pfels_isoPhotons        (new vector<float>		);  


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

     edm::Handle<reco::PFCandidateCollection> pfAllElectronsHandle;
     iEvent.getByLabel(pfAllElectrons_tag_, pfAllElectronsHandle);

     for(PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++)
     {
	  int pfflags = 0;
	  for(unsigned int i = 0; i < 17; i++)
	  {
	       if(pf_it->flag((PFCandidate::Flags)i))
		    pfflags |= (1<<i);
	  }

	  pfels_p4			->push_back(LorentzVector(pf_it->px(), pf_it->py(), pf_it->pz(), pf_it->p())	);
	  pfels_posAtEcal_p4		->push_back(LorentzVector(pf_it->positionAtECALEntrance().x(),
								  pf_it->positionAtECALEntrance().y(),
								  pf_it->positionAtECALEntrance().z(),
								  0.0)							);
	  pfels_ecalE			->push_back(pf_it->ecalEnergy()							);
	  pfels_hcalE			->push_back(pf_it->hcalEnergy()							);
	  pfels_rawEcalE		->push_back(pf_it->rawEcalEnergy()						);
	  pfels_rawHcalE		->push_back(pf_it->rawHcalEnergy()						);
	  pfels_pS1E			->push_back(pf_it->pS1Energy()							);
	  pfels_pS2E			->push_back(pf_it->pS2Energy()							);
	  pfels_deltaP			->push_back(pf_it->deltaP()							);
	  pfels_mva_epi			->push_back(pf_it->mva_e_pi()							);
	  pfels_mva_emu			->push_back(pf_it->mva_e_mu()							);
	  pfels_mva_pimu		->push_back(pf_it->mva_pi_mu()							);
	  pfels_mva_nothing_gamma	->push_back(pf_it->mva_nothing_gamma()						);
	  pfels_mva_nothing_nh		->push_back(pf_it->mva_nothing_nh()						);
	  pfels_charge			->push_back(pf_it->charge()							);
	  pfels_particleId		->push_back(pf_it->translateTypeToPdgId(pf_it->particleId())			);
	  pfels_flag			->push_back(pfflags                                                             );    
     }//loop over candidate collection

     for(unsigned int  pfit = 0; pfit < pfAllElectronsHandle->size(); pfit++)
     {
	  reco::PFCandidateRef pfref(pfAllElectronsHandle, pfit);
	  float isoc = (isoc_data)[pfref];
	  float ison = (ison_data)[pfref];
	  float isop = (isop_data)[pfref];

	  pfels_isoChargedHadrons	->push_back(isoc);
	  pfels_isoNeutralHadrons	->push_back(ison);
	  pfels_isoPhotons		->push_back(isop);
     }

     iEvent.put(pfels_p4,			"pfelsp4"		);
     iEvent.put(pfels_posAtEcal_p4,		"pfelsposAtEcalp4"	);
     iEvent.put(pfels_ecalE,			"pfelsecalE"		);
     iEvent.put(pfels_hcalE,			"pfelshcalE"		);
     iEvent.put(pfels_rawEcalE,			"pfelsrawEcalE"		);
     iEvent.put(pfels_rawHcalE,			"pfelsrawHcalE"		);
     iEvent.put(pfels_pS1E,			"pfelspS1E"		);
     iEvent.put(pfels_pS2E,			"pfelspS2E"		);
     iEvent.put(pfels_deltaP,			"pfelsdeltaP"		);
     iEvent.put(pfels_mva_epi,			"pfelsmvaepi"		);
     iEvent.put(pfels_mva_emu,			"pfelsmvaemu"		);
     iEvent.put(pfels_mva_pimu,			"pfelsmvapimu"		);
     iEvent.put(pfels_mva_nothing_gamma,	"pfelsmvanothinggamma"	);
     iEvent.put(pfels_mva_nothing_nh,		"pfelsmvanothingnh"	);
     iEvent.put(pfels_charge,			"pfelscharge"		);
     iEvent.put(pfels_particleId,		"pfelsparticleId"	);
     iEvent.put(pfels_flag,			"pfelsflag"		);
     iEvent.put(pfels_isoChargedHadrons,        "pfelsisoChargedHadrons");
     iEvent.put(pfels_isoNeutralHadrons,        "pfelsisoNeutralHadrons");
     iEvent.put(pfels_isoPhotons,               "pfelsisoPhotons"       );
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFElectronMaker);

