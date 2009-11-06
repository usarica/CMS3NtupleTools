
#include <vector>
#include "Math/VectorUtil.h"
#include "CMS2/NtupleMaker/interface/ElTkJuraIsoMaker.h"

// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"

ElTkJuraIsoMaker::ElTkJuraIsoMaker(const edm::ParameterSet& iConfig) 
{

	// cms2 inputs
	elsInputTag_            = iConfig.getParameter<edm::InputTag>("elsInputTag");
	trackInputTag_          = iConfig.getParameter<edm::InputTag>("trackInputTag");

	// track isolation configuration
	trackIsoExtRadius_    	= iConfig.getParameter<double>("trackIsoExtRadius");
	trackIsoInRadius_  	= iConfig.getParameter<double>("trackIsoInRadius");
	trackIsoJurassicWidth_ 	= iConfig.getParameter<double>("trackIsoJurassicWidth");
	trackIsoMinPt_        	= iConfig.getParameter<double>("trackIsoMinPt");
	trackIsoMind0_        	= iConfig.getParameter<double>("trackIsoMind0");
	trackIsoMinz0_        	= iConfig.getParameter<double>("trackIsoMinz0");

	//register your products
	produces<std::vector<float> >         ("elstkjuraiso").setBranchAlias("els_tkJuraIso");

}

ElTkJuraIsoMaker::~ElTkJuraIsoMaker(){}

// ------------ method called to produce the data  ------------
void ElTkJuraIsoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// electron p4
	edm::InputTag els_trk_p4_tag(elsInputTag_.label(),"elstrkp4");
	edm::Handle<std::vector<LorentzVector> > els_trk_p4_h;
	iEvent.getByLabel(els_trk_p4_tag, els_trk_p4_h);
	const std::vector<LorentzVector> *els_trk_p4 = els_trk_p4_h.product();

	// electron z0
	edm::InputTag els_z0_tag(elsInputTag_.label(),"elsz0");
	edm::Handle<std::vector<float> > els_z0_h;
	iEvent.getByLabel(els_z0_tag, els_z0_h);
	const std::vector<float> *els_z0 = els_z0_h.product();

	// track p4
	edm::InputTag trks_trk_p4_tag(trackInputTag_.label(),"trkstrkp4");
	edm::Handle<std::vector<LorentzVector> > trks_trk_p4_h;
	iEvent.getByLabel(trks_trk_p4_tag, trks_trk_p4_h);
	const std::vector<LorentzVector> *trks_trk_p4 = trks_trk_p4_h.product();

	// track d0
	edm::InputTag trks_d0_tag(trackInputTag_.label(),"trksd0corr");
	edm::Handle<std::vector<float> > trks_d0_h;
	iEvent.getByLabel(trks_d0_tag, trks_d0_h);
	const std::vector<float> *trks_d0 = trks_d0_h.product();

	// track z0
	edm::InputTag trks_z0_tag(trackInputTag_.label(),"trksz0");
	edm::Handle<std::vector<float> > trks_z0_h;
	iEvent.getByLabel(trks_z0_tag, trks_z0_h);
	const std::vector<float> *trks_z0 = trks_z0_h.product();

	std::auto_ptr<std::vector<float> > els_tkJuraIso        (new std::vector<float>);

	float isoSum;
	for (size_t j = 0; j < els_trk_p4->size(); ++j) 
	{

		isoSum = 0.0;
		for (size_t i = 0; i < trks_trk_p4->size(); ++i)
		{

			float dz0 = fabs((*trks_z0)[i] - (*els_z0)[j]);
			float dEta = (*trks_trk_p4)[i].Eta() - (*els_trk_p4)[j].Eta();
			float dR = ROOT::Math::VectorUtil::DeltaR((*els_trk_p4)[j], (*trks_trk_p4)[i]);
			float d0 = fabs((*trks_d0)[i]);
			float pT = (*trks_trk_p4)[i].Pt();

                        // is the track in the cone
                        if (dR < trackIsoInRadius_ || dR > trackIsoExtRadius_) continue;

                        // pt cut on tracks 
                        if (pT < trackIsoMinPt_) continue;

			// dz0 cut on tracks
			if (dz0 > trackIsoMinz0_) continue;

			// d0 cut
			if (d0 > trackIsoMind0_) continue;

			// is the track in the jurassic strip
			if (fabs(dEta) < trackIsoJurassicWidth_) continue;

			isoSum += pT;

		// end loop on tracks
		}

		// set the isolation value for this electron
		els_tkJuraIso->push_back(isoSum);

	// end loop on electrons
	}

	// put results into event
	iEvent.put(els_tkJuraIso, "elstkjuraiso");

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElTkJuraIsoMaker);

