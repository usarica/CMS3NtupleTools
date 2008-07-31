// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      DuplicateFilter
// 
/**\class DuplicateFilter DuplicateFilter.h CMS2/NtupleMaker/interface/DuplicateFilter.h

Description: duplicate event remover for cms2

Implementation:
- removes a given list of events (identified by a bunch of floats in
  the absence of better event identifiers)  
- prints events processed to the message logger (so duplicate lists
  can be compiled from stdout

*/
//
// Original Author:  J. Mü.
//         Created:  Wed Jun 25 19:59:33 UTC 2008  
// $Id: DuplicateFilter.h,v 1.1 2008/07/31 04:38:01 jmuelmen Exp $
//
//
#ifndef CMS2_DUPLICATEFILTER_H
#define CMS2_DUPLICATEFILTER_H

// system include files
#include <memory>
#include <string>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//
// class decleration
//

class DuplicateFilter : public edm::EDFilter {
public:
     explicit DuplicateFilter (const edm::ParameterSet&);
  
private:
     virtual bool filter(edm::Event&, const edm::EventSetup&);
     
     struct DorkyEventIdentifier {
	  // this is a workaround for not having unique event id's in MC
	  unsigned long int run, event;
	  float trks_d0;
	  float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
	  bool operator < (const DorkyEventIdentifier &) const;
     };
     // ----------member data ---------------------------
     std::set<DorkyEventIdentifier> 	eventsToFilter;

     edm::InputTag	run_tag;
     edm::InputTag	event_tag;
     edm::InputTag	trks_d0_tag;
     edm::InputTag	hyp_lt_p4_tag;
};

#endif
