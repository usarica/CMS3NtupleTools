// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      JetToElAssMaker
// 
/**\class JetToElAssMaker JetToElAssMaker.h CMS2/NtupleMaker/interface/JetToElAssMaker.h

Description: make associations between jets and muons

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun 17 20:40:42 UTC 2008
// $Id: JetToElAssMaker.h,v 1.4 2010/03/02 19:24:11 fgolf Exp $
//
//
#ifndef CMS2_JETTOELASSMAKER_H
#define CMS2_JETTOELASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class JetToElAssMaker : public edm::EDProducer {
public:
  explicit JetToElAssMaker (const edm::ParameterSet&);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  double m_minDR;
  edm::InputTag jetInputTag;
  edm::InputTag elInputTag;
};

#endif
