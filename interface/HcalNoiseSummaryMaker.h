// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/EventMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: HcalNoiseSummaryMaker.h,v 1.4 2012/04/22 21:31:32 dbarge Exp $
//
//

#ifndef NTUPLEMAKER_EVENTMAKER_H
#define NTUPLEMAKER_EVENTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TString.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

class HcalNoiseSummaryMaker : public edm::EDProducer {

 public:

  explicit HcalNoiseSummaryMaker (const edm::ParameterSet&);
  ~HcalNoiseSummaryMaker();

private:

  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::EDGetTokenT<HcalNoiseSummary> hcalNoiseSummaryToken_;
	std::string aliasprefix_;
  std::string branchprefix_;

};


#endif
