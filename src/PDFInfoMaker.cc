//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PDFInfoMaker
// 
/**\class PDFInfoMaker PDFInfoMaker.cc CMS2/NtupleMaker/src/PDFInfoMaker.cc

Description: PDF Info

Implementation:
PDF Info variable definitions:

x1 = momentum fraction of first parton (index 4)
x2 = momentum fraction of second parton (index 5)
Q = total COM energy of proton-proton collision squared (10Tev)
id1 = pdg id of first parton (index 4)
id2 = pdg id of second parton (index 5)

*/
//
// Original Author:  Warren Andrews
//         Created:  
// 


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/PDFInfoMaker.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

PDFInfoMaker::PDFInfoMaker(const edm::ParameterSet& iConfig) {

  produces<float> ("pdfinfox1" ).setBranchAlias("pdfinfo_x1" );
  produces<float> ("pdfinfox2" ).setBranchAlias("pdfinfo_x2" );
  produces<float> ("pdfinfoscale" ).setBranchAlias("pdfinfo_scale" );
  produces<int> ("pdfinfoid1" ).setBranchAlias("pdfinfo_id1" );
  produces<int> ("pdfinfoid2" ).setBranchAlias("pdfinfo_id2" );
  
}


PDFInfoMaker::~PDFInfoMaker() {}

void  PDFInfoMaker::beginJob(const edm::EventSetup&) {
}

void PDFInfoMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PDFInfoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float> pdfinfo_x1( new float );
  auto_ptr<float> pdfinfo_x2( new float );
  auto_ptr<float> pdfinfo_scale(  new float );
  auto_ptr<int> pdfinfo_id1 ( new int );
  auto_ptr<int> pdfinfo_id2 ( new int );

   // get MC particle collection
  edm::Handle<edm::HepMCProduct> hepmcHandle;
  iEvent.getByType( hepmcHandle ); //not getByLabel
  const HepMC::GenEvent* evt = hepmcHandle->GetEvent();

  //assign
  *pdfinfo_x1 = evt->pdf_info()->x1();
  *pdfinfo_x2 = evt->pdf_info()->x2();
  *pdfinfo_scale = evt->pdf_info()->scalePDF();
  *pdfinfo_id1 = evt->pdf_info()->id1();
  *pdfinfo_id2 = evt->pdf_info()->id2();
  
  
  iEvent.put(pdfinfo_x1 ,"pdfinfox1" );
  iEvent.put(pdfinfo_x2 ,"pdfinfox2" );
  iEvent.put(pdfinfo_scale ,"pdfinfoscale" );
  iEvent.put(pdfinfo_id1 ,"pdfinfoid1" );
  iEvent.put(pdfinfo_id2 ,"pdfinfoid2" );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PDFInfoMaker);





  
