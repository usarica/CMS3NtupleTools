// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventWeights
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventWeights.h,v 1.2 2008/08/05 01:24:16 jmuelmen Exp $
//
//
#ifndef CMS2_EVENTWEIGHTS_H
#define CMS2_EVENTWEIGHTS_H

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

class EventWeights : public edm::EDProducer {
public:
     explicit EventWeights (const edm::ParameterSet&);
      ~EventWeights();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag csa07InputTag;
  edm::InputTag eventInputTag;
  float topDilepton2Electron_nevts;
  float topDileptonMuonX_nevts;
  bool  IsTopDilepton2Electron;
  bool  IsTopDileptonMuonX;
  float Chowder_topDilepton2Electron_nevts;
  float Chowder_topDileptonMuonX_nevts;
     bool  IsNotSoup;
     float Xsec;
     int Nevts;

};


#endif
