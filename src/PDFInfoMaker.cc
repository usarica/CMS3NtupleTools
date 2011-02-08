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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CMS2/NtupleMaker/interface/PDFInfoMaker.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

PDFInfoMaker::PDFInfoMaker(const edm::ParameterSet& iConfig) {

  genEventInfoInputTag_ = iConfig.getParameter<std::string>("genEventInfoInputTag");
  
  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<float> (branchprefix+"x1"    ).setBranchAlias(aliasprefix_+"_x1"   );
  produces<float> (branchprefix+"x2"    ).setBranchAlias(aliasprefix_+"_x2"   );
  produces<float> (branchprefix+"scale" ).setBranchAlias(aliasprefix_+"_scale");
  produces<float> (branchprefix+"pdf1"  ).setBranchAlias(aliasprefix_+"_pdf1" );
  produces<float> (branchprefix+"pdf2"  ).setBranchAlias(aliasprefix_+"_pdf2" );
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


  using namespace edm;
  
  auto_ptr<float> pdfinfo_x1   ( new float );
  auto_ptr<float> pdfinfo_x2   ( new float );
  auto_ptr<float> pdfinfo_scale( new float );
  auto_ptr<float> pdfinfo_pdf1 ( new float ); 
  auto_ptr<float> pdfinfo_pdf2 ( new float ); 
  auto_ptr<int>   pdfinfo_id1  ( new int   );
  auto_ptr<int>   pdfinfo_id2  ( new int   );


  //get the GenEventInfoProduct
  Handle<GenEventInfoProduct> hEvtInfo;
  iEvent.getByLabel(genEventInfoInputTag_, hEvtInfo);

  // get MC particle collection
  edm::Handle<edm::HepMCProduct> hepmcHandle;
  iEvent.getByType( hepmcHandle ); //not getByLabel
  const HepMC::GenEvent* evt = 0;
  const HepMC::PdfInfo* pdfinfo = 0;
  if(!hepmcHandle.failedToGet() ) {
    evt = hepmcHandle->GetEvent();
    pdfinfo = evt->pdf_info();
  }

  //try to get using the GenEventInfoProduct
  if(!hEvtInfo.failedToGet() && hEvtInfo->hasPDF()) {
    
    const gen::PdfInfo *pdf = hEvtInfo->pdf();
    
    *pdfinfo_x1		= pdf->x.first;
    *pdfinfo_x2		= pdf->x.second;
    *pdfinfo_scale	= pdf->scalePDF;
    *pdfinfo_pdf1       = pdf->xPDF.first;
    *pdfinfo_pdf2       = pdf->xPDF.second;
    *pdfinfo_id1	= pdf->id.first;
    *pdfinfo_id2	= pdf->id.second;
    
  } else if(pdfinfo != 0) {
    
    //assign
    *pdfinfo_x1    = pdfinfo->x1();
    *pdfinfo_x2    = pdfinfo->x2();
    *pdfinfo_scale = pdfinfo->scalePDF();
    *pdfinfo_pdf1  = pdfinfo->pdf1();
    *pdfinfo_pdf2  = pdfinfo->pdf2();
    *pdfinfo_id1   = pdfinfo->id1();
    *pdfinfo_id2   = pdfinfo->id2();

  } else {

    *pdfinfo_x1		= -9999;  
    *pdfinfo_x2		= -9999;  
    *pdfinfo_scale	= -9999;
    *pdfinfo_pdf1       = -9999;
    *pdfinfo_pdf2       = -9999;
    *pdfinfo_id1	= -9999; 
    *pdfinfo_id2	= -9999;
    
    return;
  }





  
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(pdfinfo_x1   ,branchprefix+"x1"    );
  iEvent.put(pdfinfo_x2   ,branchprefix+"x2"    );
  iEvent.put(pdfinfo_scale,branchprefix+"scale" );
  iEvent.put(pdfinfo_pdf1 ,branchprefix+"pdf1"  );
  iEvent.put(pdfinfo_pdf2 ,branchprefix+"pdf2"  );
  iEvent.put(pdfinfo_id1  ,branchprefix+"id1"   );
  iEvent.put(pdfinfo_id2  ,branchprefix+"id2"   );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PDFInfoMaker);
