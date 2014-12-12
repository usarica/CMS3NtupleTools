//-*- C++ -*-
//
// Package:    LepGenFilter
// Class:      LepGenFilter
// 
/**\class LepGenFilter LepGenFilter.cc CMS2/src/LepGenFilter.cc

Description: filter for cms2

Implementation:
see header file
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  1 11:07:38 PDT 2009
// $Id: LepGenFilter.cc,v 1.3 2013/06/23 13:18:28 dalfonso Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/LepGenFilter.h"
#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;

//
// constructors and destructor
//

LepGenFilter::LepGenFilter(const edm::ParameterSet& iConfig) {

   ntupleDaughters_  = iConfig.getParameter<bool>("ntupleDaughters");
   nGenLepsRequired_ = iConfig.getParameter<int>("nGenLepsRequired");

}


LepGenFilter::~LepGenFilter() {}

void  LepGenFilter::beginJob() {
}

void LepGenFilter::endJob() {
}


// ------------ method called to produce the data  ------------
bool LepGenFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  using namespace std;
  
  InputTag genid_tag("genMaker", "genpsid");
  Handle<vector<int> > genid_h;
  iEvent.getByLabel(genid_tag, genid_h);
  const vector<int> *v_genid = genid_h.product();

  InputTag genstatus_tag("genMaker", "genpsstatus");
  Handle<vector<int> > genstatus_h;
  iEvent.getByLabel(genstatus_tag, genstatus_h);
  const vector<int> *v_genstatus = genstatus_h.product();

  InputTag genlepdaughterid_tag("genMaker", "genpslepdaughterid");
  const vector<vector<int> > *v_genlepdaughterid;
  if(ntupleDaughters_) {
    Handle<vector<vector<int> > > genlepdaughterid_h;
    iEvent.getByLabel(genlepdaughterid_tag, genlepdaughterid_h);
    v_genlepdaughterid= genlepdaughterid_h.product();
  }
 

  int nGenLeps = 0;
  for(unsigned int idx = 0; idx < v_genid->size(); idx++) {
    
    int pid = abs(v_genid->at(idx));

    if(v_genstatus->at(idx) != 3) 
      continue;
    
    //if electron of muon in doc line, sum 
    if(pid == 11 || pid == 13)
      nGenLeps++;

    if(ntupleDaughters_) {
      //if a tau whose daughter is a e/mu, sum
      if(pid == 15) {
	for(unsigned int daughter = 0; daughter < v_genlepdaughterid->at(idx).size(); daughter++) {
	  int daughterid = abs(v_genlepdaughterid->at(idx).at(daughter));
	  
	  if(daughterid == 11 || daughterid == 13)
	    nGenLeps++;
	}//daughter loop
      }//if(pid == 15)
    }	    
    
  }//genid loop

  if(nGenLeps < nGenLepsRequired_)
    return false;

  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(LepGenFilter);
 
