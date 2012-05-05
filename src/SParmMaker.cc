// -*- C++ -*-
//
// Package:    SParmMaker
// Class:      SParmMaker
// 
/**\class SParmMaker SParmMaker.cc CMS2/NtupleMaker/src/SParmMaker.cc

   Description: copy SUSY mSUGRA parameters into the EDM event tree

   Implementation:
   - extract and fill variables
*/
//
// Original Ben Hooberman
// Created:  Wed Mar  24 12:23:38 CDT 2010
// 
//
//


// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "CMS2/NtupleMaker/interface/SParmMaker.h"

#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

//
// class declaration
//

//
// constructors and destructor
//

SParmMaker::SParmMaker(const edm::ParameterSet& iConfig) {
  
  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
	
  // parameters from configuration
  sparm_inputTag = iConfig.getParameter<edm::InputTag>("sparm_inputTag");

  // sparm names from configuration
  vsparms_ = iConfig.getUntrackedParameter<std::vector<std::string> >("vsparms");

  // product of this EDProducer
  produces<std::vector<TString> > (branchprefix+"comment"          ).setBranchAlias(aliasprefix_+"_comment"          ); //even though we need a single string, it must be stored as a vector
  produces<std::vector<TString> > (branchprefix+"names"            ).setBranchAlias(aliasprefix_+"_names"            );
  produces<std::vector<float> >   (branchprefix+"values"           ).setBranchAlias(aliasprefix_+"_values"           );
  produces<int>                   (branchprefix+"subProcessId"     ).setBranchAlias(aliasprefix_+"_subProcessId"     );
  produces<float>                 (branchprefix+"weight"           ).setBranchAlias(aliasprefix_+"_weight"           );
  produces<float>                 (branchprefix+"pdfWeight1"       ).setBranchAlias(aliasprefix_+"_pdfWeight1"       );
  produces<float>                 (branchprefix+"pdfWeight2"       ).setBranchAlias(aliasprefix_+"_pdfWeight2"       );
  produces<float>                 (branchprefix+"pdfScale"         ).setBranchAlias(aliasprefix_+"_pdfScale"         );
  produces<float>                 (branchprefix+"filterEfficiency" ).setBranchAlias(aliasprefix_+"_filterEfficiency" );
  produces<float>                 (branchprefix+"xsec"             ).setBranchAlias(aliasprefix_+"_xsec"             );
 
}

SParmMaker::~SParmMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SParmMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<std::vector<TString> > sparm_comment          (new std::vector<TString>(1) ); 
  std::auto_ptr<std::vector<TString> > sparm_names            (new std::vector<TString>    );	
  std::auto_ptr<std::vector<float> >   sparm_values           (new std::vector<float>      );
  std::auto_ptr<int>                   sparm_subProcessId     (new int(-9999)              );
  std::auto_ptr<float>                 sparm_weight           (new float(-9999.)           );
  std::auto_ptr<float>                 sparm_pdfWeight1       (new float(-9999.)           );
  std::auto_ptr<float>                 sparm_pdfWeight2       (new float(-9999.)           );  
  std::auto_ptr<float>                 sparm_pdfScale         (new float(-9999.)           );
  std::auto_ptr<float>                 sparm_filterEfficiency (new float(1.)               );
  std::auto_ptr<float>                 sparm_xsec             (new float(-9999.)           );

  
  // fill the user supplied parameter names
  for(size_t i=0; i<vsparms_.size(); i++){
	sparm_names->push_back( TString(vsparms_[i].c_str()) );
  }
  
  // first try cms sparm comments
  edm::Handle<LHEEventProduct> sparm_handle;  
  iEvent.getByLabel(sparm_inputTag, sparm_handle);
  if( sparm_handle.isValid() ){
	for (std::vector<std::string>::const_iterator it = sparm_handle->comments_begin(); it != sparm_handle->comments_end(); it++) {      
	  TString model_comment(*it);

	  // check if sparm comment is in expected format
	  TObjArray* space_tokens = model_comment.Tokenize(" ");
	  int space_tokens_length=space_tokens->GetEntries();
	  for(int i=space_tokens_length-1; i>=0; i--){
		if( ((TObjString*)space_tokens->At(i))->GetString() == "\n" ){
		  space_tokens_length--; //if they had designed TString sensibly, we wouldn't need to worry about entries in TObjstring that are new lines...
		}
	  }

	  // check to make sure the comment is in the right format
	  if( space_tokens_length < 3 ||
		  ((TObjString*)space_tokens->At(0))->GetString() != "#" ||
		  ((TObjString*)space_tokens->At(1))->GetString() != "model" ){
		throw cms::Exception(Form("SParmMaker: Sparm comment not in the form expected. Comment is \"%s\".",model_comment.Data()));
	  }
	  if(space_tokens_length >= 4){ // this isn't so general, but we don't know what else to expect in the string
			*sparm_xsec             = ((TObjString*)space_tokens->At(3))->GetString().Atof();
		if(space_tokens_length >= 5){
		  *sparm_filterEfficiency = ((TObjString*)space_tokens->At(4))->GetString().Atof();
		}
	  }
	
	  (*sparm_comment)[0]=TString(*it);
	  TString model_params = ((TObjString*)space_tokens->At(2))->GetString();
	  TObjArray* tokens = model_params.Tokenize("_");
	  for(int i = 1; i < tokens->GetEntries(); i++){ //the crap before the first "_" is not a parameter
		sparm_values->push_back( ((TObjString*)tokens->At(i))->GetString().Atof() );
	  }
	  delete space_tokens;
	  delete tokens;
	}
	  
	// now, get info about this event
	const lhef::HEPEUP lhe_info = sparm_handle->hepeup();
	*sparm_subProcessId = lhe_info.IDPRUP;
	*sparm_weight       = lhe_info.XWGTUP;
	*sparm_pdfWeight1   = lhe_info.XPDWUP.first;
	*sparm_pdfWeight2   = lhe_info.XPDWUP.second;
	  
	const gen::PdfInfo* pdf_info = sparm_handle->pdf();
	if (pdf_info != 0) {
	  *sparm_pdfScale     = pdf_info->scalePDF;
	}
	else {
	  *sparm_pdfScale = lhe_info.SCALUP;
	}
  }
  else{ // next try our custom susy branches if sparm comment is not present
	edm::Handle<double> susyScan_handles[6];
	for(int i=0; i<6; i++){
	  iEvent.getByLabel(Form("susyScanP%i",(i+1)), susyScan_handles[i]);
	  if( !susyScan_handles[i].isValid() ){ break; }
	  sparm_values->push_back(*(susyScan_handles[i]));
	}
  }



  if(sparm_values->size() != sparm_names->size()){
	// We want to make damn sure that the 2 vectors have a 1 to 1 mapping. If the user doesn't screw up, this exception should never be encountered.
	throw cms::Exception(Form("SParmMaker: Size of the vector containing sparm values (size=%i) does not match size of the vector containing sparm names (size=%i).",int(sparm_values->size()),int(sparm_names->size())) );
  }
  if(sparm_comment->size() != 1){
	throw cms::Exception(Form("SparmMaker: Some Jabroney tried to store too many values in the sparm_comment vector (size=%i). This should always be a size of 1.",int(sparm_comment->size())));
  }


  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
  iEvent.put(sparm_comment          ,branchprefix+"comment"          );
  iEvent.put(sparm_names            ,branchprefix+"names"            );
  iEvent.put(sparm_values           ,branchprefix+"values"           );
  iEvent.put(sparm_weight           ,branchprefix+"weight"           );
  iEvent.put(sparm_pdfWeight1       ,branchprefix+"pdfWeight1"       );
  iEvent.put(sparm_pdfWeight2       ,branchprefix+"pdfWeight2"       );
  iEvent.put(sparm_pdfScale         ,branchprefix+"pdfScale"         );
  iEvent.put(sparm_subProcessId     ,branchprefix+"subProcessId"     );
  iEvent.put(sparm_filterEfficiency ,branchprefix+"filterEfficiency" );
  iEvent.put(sparm_xsec             ,branchprefix+"xsec"             );
}


// ------------ method called once each job just before starting event loop  ------------
void SParmMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void SParmMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SParmMaker);
