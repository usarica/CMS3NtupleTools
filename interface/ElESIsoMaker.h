#ifndef NTUPLEMAKER_ElESIsoMaker_H
#define NTUPLEMAKER_ElESIsoMaker_H
// -*- C++ -*-
// $Id: ElESIsoMaker.h,v 1.1 2009/08/01 09:57:47 dlevans Exp $

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

};

#endif
