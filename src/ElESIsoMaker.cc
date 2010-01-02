// -*- C++ -*-
// $Id: ElESIsoMaker.cc,v 1.3 2010/01/02 02:47:59 kalavase Exp $

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CMS2/NtupleMaker/interface/ElESIsoMaker.h"
#include "CMS2/NtupleMaker/interface/ESCluster.h"

ElESIsoMaker::ElESIsoMaker(const edm::ParameterSet& iConfig)
{
	produces<std::vector<float> >  ("elsesJuraIso03").setBranchAlias("els_esJuraIso03");
	produces<std::vector<float> >  ("elsesJuraIso04").setBranchAlias("els_esJuraIso04");
        produces<std::vector<float> >  ("elsesJuraVeto").setBranchAlias("els_esJuraVeto");

	electronsInputTag_ 	= iConfig.getParameter<edm::InputTag>("electronsInputTag");
	esHitsInputTag_ 	= iConfig.getParameter<edm::InputTag>("esHitsInputTag");

        intRadius_ = iConfig.getParameter<double>("intRadius");
        etaSlice_ = iConfig.getParameter<double>("etaSlice");
	useNumCrystals_ = iConfig.getParameter<bool>("useNumCrystals");

        clusterAlgo_ = new ESClusterAlgo();
        topology_ = 0;
	
}


ElESIsoMaker::~ElESIsoMaker()
{
        delete clusterAlgo_;
        if (topology_) delete topology_;
}

void
ElESIsoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
	std::auto_ptr<std::vector<float> >  els_esJuraIso03( new std::vector<float> ) ;
        std::auto_ptr<std::vector<float> >  els_esJuraIso04( new std::vector<float> ) ;
        std::auto_ptr<std::vector<float> >  els_esJuraVeto( new std::vector<float> ) ;

	// get the electrons
        edm::Handle<reco::GsfElectronCollection> electronHandle;
        iEvent.getByLabel(electronsInputTag_, electronHandle);
        const reco::GsfElectronCollection *electrons = electronHandle.product();

	// get the ES hits
        edm::Handle<EcalRecHitCollection> ESHitsHandle;
        iEvent.getByLabel(esHitsInputTag_, ESHitsHandle);
        const EcalRecHitCollection *ESHits = ESHitsHandle.product();

        // get Calo Geometry
        edm::ESHandle<CaloGeometry> pG;
        iSetup.get<CaloGeometryRecord>().get(pG);
        const CaloGeometry* caloGeom = pG.product();

	// get the preshower topology (needed to make clusters)
        CaloSubdetectorTopology *topology_ = 0;
	  if (!topology_) topology_  = new EcalPreshowerTopology(pG); //offending line

	// make the ES clusters for this event
        std::vector<ESCluster> clusters;
        clusterAlgo_->cluster(ESHits, caloGeom, topology_, clusters);

        for (size_t j = 0; j < electrons->size(); ++j)
        {

		// get the electron calorimeter position
                float scEta = (*electrons)[j].superCluster()->eta();
                float scPhi = (*electrons)[j].superCluster()->phi();

		// only care about endcap electrons
		if (fabs(scEta) < 1.5) {
                	els_esJuraIso03->push_back(0.0);
        	        els_esJuraIso04->push_back(0.0);
	                els_esJuraVeto->push_back(0.0);
			continue;
		}
	
		// isolation variables
		float iso03 = 0.0;
		float iso04 = 0.0;
		float veto = 0.0;

		// identify which ES clusters are in the isolation cone
                for (size_t clus = 0; clus < clusters.size(); ++clus)
                {

                        float dEta = clusters[clus].eta() - scEta;
                        float dPhi = acos(cos(clusters[clus].phi() - scPhi));
                        float deltaR = sqrt(dEta*dEta + dPhi*dPhi);

			bool inVetoRegion = false;
                    	if(useNumCrystals_) {
				if ( fabs(dEta) < 0.00864*fabs(sinh(scEta))*etaSlice_ 
						&& deltaR < 0.00864*fabs(sinh(scEta))*intRadius_) 
					inVetoRegion = true;  
			} else {
	                        if (deltaR > intRadius_ && fabs(dEta) > etaSlice_)
					inVetoRegion = true;
			}

			if (!inVetoRegion) {
				if (deltaR < 0.4) iso04 += clusters[clus].energy() / cosh( clusters[clus].eta() );
				if (deltaR < 0.3) iso03 += clusters[clus].energy() / cosh( clusters[clus].eta() );
			}
			else {
				veto += clusters[clus].energy() / cosh( clusters[clus].eta() );
			}

                } // end loop on ES clusters

		els_esJuraIso03->push_back(iso03);
                els_esJuraIso04->push_back(iso04);
                els_esJuraVeto->push_back(veto);
		

	} // end loop on electrons
	iEvent.put(els_esJuraIso03, "elsesJuraIso03");
        iEvent.put(els_esJuraIso04, "elsesJuraIso04");
        iEvent.put(els_esJuraVeto, "elsesJuraVeto");
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElESIsoMaker);
