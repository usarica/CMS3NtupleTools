// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElectronDuplicateRemover
// 
/**\class ElectronDuplicateRemover ElectronDuplicateRemover.cc CMS2/NtupleMaker/src/ElectronDuplicateRemover.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
#ifndef CMS2_SCMAKER_H
#define CMS2_SCMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

//
// class declaration
//

class ElectronDuplicateRemover : public edm::EDProducer {
public:
     explicit ElectronDuplicateRemover (const edm::ParameterSet&);
  
private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     // a comparitor function for closest electron to E/P = 1
     //bool betterElectron(const reco::GsfElectron &e1, const reco::GsfElectron &e2);

     // ----------member data ---------------------------

     // primary vertex collection
     edm::InputTag electronsInputTag_;

};

class BetterElectron : public std::binary_function<reco::GsfElectron, reco::GsfElectron, bool>
{
     public:
          bool operator() (const reco::GsfElectron &e1, const reco::GsfElectron &e2)
          {
               return (fabs(e1.eSuperClusterOverP() - 1) < fabs(e2.eSuperClusterOverP() - 1));
          }
};

#endif

