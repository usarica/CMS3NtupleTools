// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuToElsAssMaker
// 
/**\class MuToElsAssMaker MuToElsAssMaker.cc CMS2/MuToElsAssMaker/src/MuToElsAssMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToElsAssMaker.h,v 1.4 2010/03/03 04:20:03 kalavase Exp $
//
//
#ifndef CMS2_MUTOELASSMAKER_H
#define CMS2_MUTOELASSMAKER_H

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;

//
// class declaration
//

class MuToElsAssMaker : public edm::EDProducer {
public:
     explicit MuToElsAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
     double m_minDR;
     std::string aliasprefix_;
     edm::EDGetTokenT<std::vector<LorentzVector> > musToken_;
     edm::EDGetTokenT<std::vector<LorentzVector> > elsToken_;
};


#endif
