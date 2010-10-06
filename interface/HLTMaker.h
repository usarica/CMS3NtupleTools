// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HLTMaker
// 
/**\class HLTMaker CMS2/NtupleMaker/src/HLTMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: HLTMaker.h,v 1.8 2010/10/06 13:52:42 jribnik Exp $
//
//
#ifndef NTUPLEMAKER_HLTMAKER_H
#define NTUPLEMAKER_HLTMAKER_H

// system include files
#include <algorithm>
#include <string>
#include <vector>

// user include files
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TRegexp.h"
#include "TString.h"

class HLTMaker : public edm::EDProducer {
public:
  explicit HLTMaker(const edm::ParameterSet&);
  ~HLTMaker() {}
  
private:
  virtual void beginRun(edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}
  
  void fillTriggerObjectInfo(unsigned int,
			     std::vector<int>&,
			     std::vector<math::XYZTLorentzVectorF>&) const;
  bool doPruneTriggerName(const std::string&) const;
  
  edm::Handle<edm::TriggerResults> triggerResultsH_;
  edm::Handle<trigger::TriggerEvent> triggerEventH_;
  HLTConfigProvider hltConfig_;
  edm::ProcessHistory edmPH_;
  
  std::string processName_;
  bool fillTriggerObjects_;
  std::vector<std::string> prunedTriggerNames_;
  TString processNamePrefix_;
  std::string aliasprefix_;
};

#endif
