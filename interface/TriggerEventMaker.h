// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TriggerEventMaker
// 
/**\class TriggerEventMaker CMS2/NtupleMaker/src/TriggerEventMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TriggerEventMaker.h,v 1.1 2008/12/16 03:47:08 slava77 Exp $
//
//
#ifndef NTUPLEMAKER_TriggerEventMaker_H
#define NTUPLEMAKER_TriggerEventMaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//
// class decleration
//

class TriggerEventMaker : public edm::EDProducer {
public:
     explicit TriggerEventMaker (const edm::ParameterSet&);
      ~TriggerEventMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  void fillFilterInfo(const trigger::TriggerEvent* tevCP, const edm::InputTag& fName,
					 std::auto_ptr<std::vector<int> >& tidV, std::auto_ptr<std::vector<int> >& idV,
					 std::auto_ptr<std::vector<math::XYZTLorentzVector> >& p4V) const;
  // ----------member data ---------------------------
};


#endif
