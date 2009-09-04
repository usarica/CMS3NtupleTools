//-*- C++ -*-
//
// Package:    TriggerFilter
// Class:      TriggerFilter
// 
/**\class TriggerFilter TriggerFilter.cc CMS2/src/TriggerFilter.cc

Description: filter for cms2

 */
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/TriggerFilter.h"

using namespace std;

//
// constructors and destructor
//

TriggerFilter::TriggerFilter(const edm::ParameterSet& iConfig) {

	filterExpressions_ = iConfig.getParameter<std::vector<std::string> >("filterExpressions");
	menu_ = TString(iConfig.getParameter<std::string>("menu"));
	menu_.ToLower();
	hlt_ = iConfig.getParameter<bool>("hlt");
}

TriggerFilter::~TriggerFilter() {}

void  TriggerFilter::beginJob(const edm::EventSetup&) {
}

void TriggerFilter::endJob() {
}

// ------------ method called to produce the data  ------------
bool TriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	//loop over input trigger names
	for ( vector<std::string>::const_iterator i = filterExpressions_.begin(), end = filterExpressions_.end();
			i != end; ++i) {

		// if the filter is the OR of HLT bits then
		// consider if the requested HLT passed
		if (hlt_) {
			if (passHLTTrigger(iEvent, *i)) return true;

		// otherwise, consider if the requested
		// L1 trigger passed
		} else {
			if (passL1Trigger(iEvent, *i)) return true;
		}

	}

	// if none of the requested triggers passed
	// then return false
	return false;

}

bool TriggerFilter::passHLTTrigger(edm::Event& iEvent, TString trigName) {

	edm::Handle<std::vector<TString> > trigNamesHandle;
	iEvent.getByLabel("hltMaker", Form("%strigNames", menu_.Data()), trigNamesHandle);
	const std::vector<TString> *evt_HLT_trigNames = trigNamesHandle.product();

	int trigIndx;
	vector<TString>::const_iterator begin_it = evt_HLT_trigNames->begin();
	vector<TString>::const_iterator end_it = evt_HLT_trigNames->end();
	vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
	if(found_it != end_it)
		trigIndx = found_it - begin_it;
	else {
		cout << "Cannot find Trigger " << trigName << endl;
		return 0;
	}

	if(trigIndx <= 31) {
		unsigned int bitmask = 1;
		bitmask <<= trigIndx;
		return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 32 && trigIndx <= 63) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 32);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 64 && trigIndx <= 95) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 64);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 96 && trigIndx <= 127) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 96);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 128 && trigIndx <= 159) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 128);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 160 && trigIndx <= 191) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 160);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 192 && trigIndx <= 223) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 192);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	if(trigIndx >= 224 && trigIndx <= 255) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 224);
                return getTriggerBit(iEvent, edm::InputTag("hltMaker", Form("%sbits", menu_.Data()))) & bitmask;
	}
	return 0;
}

bool TriggerFilter::passL1Trigger(edm::Event& iEvent, TString trigName) {

        edm::Handle<std::vector<TString> > trigNamesHandle;
        iEvent.getByLabel("l1Maker", "l1bitstrigNames", trigNamesHandle);
        const std::vector<TString> *evt_L1_trigNames = trigNamesHandle.product();

	int trigIndx;
	vector<TString>::const_iterator begin_it = evt_L1_trigNames->begin();
	vector<TString>::const_iterator end_it = evt_L1_trigNames->end();
	vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
	if(found_it != end_it)
		trigIndx = found_it - begin_it;
	else {
		cout << "Cannot find Trigger " << trigName << endl;
		return 0;
	}

	if(trigIndx <= 31) { 
		unsigned int bitmask = 1;
		bitmask <<= trigIndx;
                return getTriggerBit(iEvent, edm::InputTag("l1Maker", "l1bits1")) & bitmask;
	}
	if(trigIndx >= 32 && trigIndx <= 63) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 32);
                return getTriggerBit(iEvent, edm::InputTag("l1Maker", "l1bits2")) & bitmask;
	}
	if(trigIndx >= 64 && trigIndx <= 95) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 64);
                return getTriggerBit(iEvent, edm::InputTag("l1Maker", "l1bits3")) & bitmask;
	}
	if(trigIndx >= 96 && trigIndx <= 127) {
		unsigned int bitmask = 1;
		bitmask <<= (trigIndx - 96);
                return getTriggerBit(iEvent, edm::InputTag("l1Maker", "l1bits4")) & bitmask;
	}
	return 0;
}

int TriggerFilter::getTriggerBit(edm::Event& iEvent, edm::InputTag bitInputTag)
{
	edm::Handle<int> bitHandle;
	iEvent.getByLabel(bitInputTag, bitHandle);
	return *bitHandle;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerFilter);

