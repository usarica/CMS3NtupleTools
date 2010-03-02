// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HypDilepVertexMaker
//
/**\class HypDilepVertexMaker HypDilepVertexMaker.h CMS2/NtupleMaker/interface/HypDilepVertexMaker.h

 Description: fit the hyp leptons to a single vertex and add
              relevant branches to output
*/
//
// Original Author:  slava77
//         Created:  Thu Dec  17 11:07:38 CDT 2009
// $Id: HypDilepVertexMaker.h,v 1.2 2010/03/02 19:24:11 fgolf Exp $
//
//
#ifndef CMS2_HYPDILEPVERTEXMAKER_H
#define CMS2_HYPDILEPVERTEXMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// class declaration
//

class HypDilepVertexMaker : public edm::EDProducer {
 public:
  explicit HypDilepVertexMaker (const edm::ParameterSet&);

 private:
  //  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  //  virtual void endJob();

  // ----------member data --------------------------- 
  edm::InputTag recomuonsInputTag_;
  edm::InputTag cms2muonsInputTag_;

  edm::InputTag recoelectronsInputTag_;
  edm::InputTag cms2electronsInputTag_;

  edm::InputTag hypInputTag_;
};

#endif //CMS2_HYPDILEPVERTEXMAKER_H
