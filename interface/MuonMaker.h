// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuonMaker
// 
/**\class MuonMaker MuonMaker.cc CMS2/MuonMaker/src/MuonMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonMaker.h,v 1.13 2010/03/02 19:24:11 fgolf Exp $
//
//
#ifndef CMS2_MUONMAKER_H
#define CMS2_MUONMAKER_H

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"

//
// class declaration
//

class MuonMaker : public edm::EDProducer {
public:
     explicit MuonMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
  
      // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag beamSpotInputTag;
  std::string tevMuonsName;

};


#endif
