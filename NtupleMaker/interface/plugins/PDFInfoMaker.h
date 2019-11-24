// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PDFInfoMaker
// 
/**\class PDFInfoMaker PDFInfoMaker.cc CMS3/NtupleMaker/src/PDFInfoMaker.cc

   Description: <one line class summary>

   Implementation:
   see CMS2/src/PDFInfoMaker.cc

*/
//
// Original Author: Warren Andrews
//         Created: Thu Mar  12  2009
// 

#ifndef CMS2_PDFINFOMAKER_H
#define CMS2_PDFINFOMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//
// class decleration
//

class PDFInfoMaker : public edm::stream::EDProducer<> {
public:
  explicit PDFInfoMaker (const edm::ParameterSet&);
  ~PDFInfoMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  std::string genEventInfoInputTag_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;
  edm::EDGetTokenT<edm::HepMCProduct> hepmcToken;
  std::string aliasprefix_;
};

#endif
