#ifndef CMS2_ELTKJURAISOMAKER_H
#define CMS2_ELTKJURAISOMAKER_H

// -*- C++ -*-
//

#include <vector>
#include <memory>

// framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;

class ElTkJuraIsoMaker : public edm::EDProducer {

	public:

		explicit ElTkJuraIsoMaker(const edm::ParameterSet&);
		~ElTkJuraIsoMaker();
		virtual void produce(edm::Event&, const edm::EventSetup&);

	private:

		edm::InputTag elsInputTag_;
		edm::InputTag trackInputTag_;

		double trackIsoExtRadius_;   
		double trackIsoInRadius_;
		double trackIsoJurassicWidth_;
		double trackIsoMinPt_;
		double trackIsoMind0_;
		double trackIsoMinz0_;

};

#endif
