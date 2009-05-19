// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFMETMaker
// 
/**\class PFMETMaker.cc CMS2/NtupleMaker/src/PFMETMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMETMaker.h,v 1.1 2009/05/19 23:20:02 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_PFMETMAKER_H
#define NTUPLEMAKER_PFMETMAKER_H

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

class PFMETMaker : public edm::EDProducer {
public:
  explicit PFMETMaker (const edm::ParameterSet&);
  ~PFMETMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag pfMetInputTag;
};


#endif
