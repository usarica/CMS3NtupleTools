#ifndef NTUPLEMAKER_ElESIsoMaker_H
#define NTUPLEMAKER_ElESIsoMaker_H
// -*- C++ -*-
// $Id: ElESIsoMaker.h,v 1.2 2010/03/03 04:19:19 kalavase Exp $

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/ESClusterAlgo.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

class ElESIsoMaker : public edm::EDProducer {

	public:
   		explicit ElESIsoMaker(const edm::ParameterSet&);
   		~ElESIsoMaker();
 private:

   	void produce(edm::Event&, const edm::EventSetup&);

        edm::InputTag 	electronsInputTag_;
        edm::InputTag	esHitsInputTag_;

	double intRadius_;
	double etaSlice_;
	bool useNumCrystals_;

        ESClusterAlgo *clusterAlgo_;
        EcalPreshowerTopology *topology_;

	std::string aliasprefix_;
};

#endif
