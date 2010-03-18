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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CMS2/NtupleMaker/interface/PDFInfoMaker.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

PDFInfoMaker::PDFInfoMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<float> (branchprefix+"x1"    ).setBranchAlias(aliasprefix_+"_x1"   );
  produces<float> (branchprefix+"x2"    ).setBranchAlias(aliasprefix_+"_x2"   );
  produces<float> (branchprefix+"scale" ).setBranchAlias(aliasprefix_+"_scale");
  produces<int>   (branchprefix+"id1"   ).setBranchAlias(aliasprefix_+"_id1"  );
  produces<int>   (branchprefix+"id2"   ).setBranchAlias(aliasprefix_+"_id2"  );  
}

PDFInfoMaker::~PDFInfoMaker()
{
}

void  PDFInfoMaker::beginJob()
{
}

void PDFInfoMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void PDFInfoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float> pdfinfo_x1   ( new float );
  auto_ptr<float> pdfinfo_x2   ( new float );
  auto_ptr<float> pdfinfo_scale( new float );
  auto_ptr<int>   pdfinfo_id1  ( new int   );
  auto_ptr<int>   pdfinfo_id2  ( new int   );

  // get MC particle collection
  edm::Handle<edm::HepMCProduct> hepmcHandle;
  iEvent.getByType( hepmcHandle ); //not getByLabel
  try {
    iEvent.getByType( hepmcHandle ); //not getByLabel
  }
  catch ( cms::Exception& ex ) {
    edm::LogError("PDFInfoMakerError") << "Error! can't get the gen product";
  }

  const HepMC::GenEvent* evt = hepmcHandle->GetEvent();
  const HepMC::PdfInfo* pdfInfo = 0;
  pdfInfo = evt->pdf_info();

  if (pdfInfo == 0) {
    edm::LogError("PDFInfoMakerError") << "Error! No PDF Info";
    return;
  }

  //assign
  *pdfinfo_x1    = pdfInfo->x1();
  *pdfinfo_x2    = pdfInfo->x2();
  *pdfinfo_scale = pdfInfo->scalePDF();
  *pdfinfo_id1   = pdfInfo->id1();
  *pdfinfo_id2   = pdfInfo->id2();
  
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(pdfinfo_x1   ,branchprefix+"x1"    );
  iEvent.put(pdfinfo_x2   ,branchprefix+"x2"    );
  iEvent.put(pdfinfo_scale,branchprefix+"scale" );
  iEvent.put(pdfinfo_id1  ,branchprefix+"id1"   );
  iEvent.put(pdfinfo_id2  ,branchprefix+"id2"   );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PDFInfoMaker);
