// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MITConversionMaker
// 
/**\class MITConversionMaker MITConversionMaker.cc CMS2/NtupleMaker/src/MITConversionMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Thu Jun  10 11:07:38 CDT 2010
// $Id: MITConversionMaker.h,v 1.3 2010/07/08 14:35:47 kalavase Exp $
//
//
#ifndef NTUPLEMAKER_FLAVORHISTORYPRODUCER_H
#define NTUPLEMAKER_FLAVORHISTORYPRODUCER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class MITConversionMaker : public edm::EDProducer {
public:
  explicit MITConversionMaker (const edm::ParameterSet&);
  ~MITConversionMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  // ----------member data ---------------------------
  edm::InputTag elsInputTag_;
  edm::InputTag mitConversionsTag_;
  edm::InputTag ctfTrksInputTag_;
  edm::InputTag beamSpotTag_;
  
};

#endif

