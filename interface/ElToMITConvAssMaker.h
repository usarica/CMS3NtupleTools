// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElToMITConvAssMaker
// 
/**\class ElToMITConvAssMaker ElToMITConvAssMaker.h CMS2/NtupleMaker/interface/ElToMITConvAssMaker.h

Description: make associations between electrons and jets

*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Wed Jun 13 18:32:24 UTC 2010
// $Id: ElToMITConvAssMaker.h,v 1.2 2010/06/14 11:56:08 kalavase Exp $
//
//
#ifndef CMS2_ELTOMITCONVERSIONASSM_H
#define CMS2_ELTOJETASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class ElToMITConvAssMaker : public edm::EDProducer {
public:
  explicit ElToMITConvAssMaker (const edm::ParameterSet&);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag elsInputTag_;
  edm::InputTag mitConvMakerInputTag_;

};

#endif
