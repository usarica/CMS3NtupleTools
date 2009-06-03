//-*- C++ -*-
//
// Package:    DilepGenFilter
// Class:      DilepGenFilter
// 
/**\class DilepGenFilter DilepGenFilter.cc CMS2/src/DilepGenFilter.cc

Description: filter for cms2

Implementation:
see header file
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  1 11:07:38 PDT 2009
// $Id: DilepGenFilter.cc,v 1.1 2009/06/03 00:21:30 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/DilepGenFilter.h"

#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;

//
// constructors and destructor
//

DilepGenFilter::DilepGenFilter(const edm::ParameterSet& iConfig) {
}


DilepGenFilter::~DilepGenFilter() {}

void  DilepGenFilter::beginJob(const edm::EventSetup&) {
}

void DilepGenFilter::endJob() {
}


// ------------ method called to produce the data  ------------
bool DilepGenFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
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
  Handle<vector<int> > genlepdaughterid_h;
  iEvent.getByLabel(genlepdaughterid_tag, genlepdaughterid_h);
  const vector<int> *v_genlepdaughterid = genlepdaughterid_h.product();
  

  int nGenLeps = 0;
  for(unsigned int idx = 0; idx < v_genid->size(); idx++) {
    
    int pid = abs(v_genid->at(idx));

    if(v_genstatus->at(idx) != 3) 
      continue;
    
    //if electron of muon in doc line, sum 
    if(pid == 11 || pid == 13)
      nGenLeps++;


    //if a tau whose daughter is a e/mu, sum
    if(pid == 15) {
      for(unsigned int daughter = 0; daughter < v_genlepdaughterid->size(); daughter++) {
	int daughterid = abs(v_genlepdaughterid->at(daughter));
	
	if(daughterid == 11 || daughterid == 13)
	  nGenLeps++;
      }//daughter loop
    }//if(pid == 15)
	    
    
  }//genid loop

  cout << nGenLeps << endl;
  if(nGenLeps < 2)
    return false;

  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DilepGenFilter);





  
