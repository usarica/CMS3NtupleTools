// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TriggerFilter
// 
/**\class TriggerFilter TriggerFilter.h CMS2/NtupleMaker/interface/TriggerFilter.h

Description: Filter on the OR of trigger bits

Implementation:
- get list of names of triggers as input
- event passes if any of these triggers passed

 */
//
//
#ifndef CMS2_TRIGGERFILTER_H
#define CMS2_TRIGGERFILTER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "TString.h"

//
// class decleration
//

class TriggerFilter : public edm::EDFilter {
	public:



		explicit TriggerFilter (const edm::ParameterSet&);
		~TriggerFilter();

	private:
		virtual void beginJob(const edm::EventSetup&) ;
		virtual bool filter(edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		bool passHLTTrigger(edm::Event& iEvent, TString trigName);
		bool passL1Trigger(edm::Event& iEvent, TString trigName);
		int getTriggerBit(edm::Event& iEvent, edm::InputTag bitInputTag);

		// ----------member data ---------------------------
		std::vector<std::string> filterExpressions_;
		TString menu_;
		bool hlt_;

};

#endif

