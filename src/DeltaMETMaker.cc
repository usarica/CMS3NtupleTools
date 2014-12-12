//-*- C++ -*-
//
// Package:    DeltaMETMaker
// Class:      DeltaMETMaker
// 
/**\class DeltaMETMaker DeltaMETMaker.cc CMS2/DeltaMETMaker/src/DeltaMETMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: DeltaMETMaker.cc,v 1.1 2010/05/19 21:58:59 fgolf Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/DeltaMETMaker.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

DeltaMETMaker::DeltaMETMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
       branchprefix.replace(branchprefix.find("_"),1,"");

  produces<float> (branchprefix+"dmetx" ).setBranchAlias(aliasprefix_+"_dmetx" );
  produces<float> (branchprefix+"dmety" ).setBranchAlias(aliasprefix_+"_dmety" );
  produces<float> (branchprefix+"dsumet").setBranchAlias(aliasprefix_+"_dsumet");
  
  cms2_metInputTag    = iConfig.getParameter<edm::InputTag>("cms2_metInputTag_"   );       
  cms2_metphiInputTag = iConfig.getParameter<edm::InputTag>("cms2_metphiInputTag_");     
  cms2_sumetInputTag  = iConfig.getParameter<edm::InputTag>("cms2_sumetInputTag_" );   
  metInputTag         = iConfig.getParameter<edm::InputTag>("metInputTag_"        ); 
}

DeltaMETMaker::~DeltaMETMaker()
{
}

void  DeltaMETMaker::beginJob()
{
}

void DeltaMETMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void DeltaMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float>   evt_dmetx     (new float);
  auto_ptr<float>   evt_dmety     (new float);
  auto_ptr<float>   evt_dsumet    (new float);

  edm::Handle<float> cms2_met_h;
  iEvent.getByLabel(cms2_metInputTag, cms2_met_h);

  edm::Handle<float> cms2_metphi_h;
  iEvent.getByLabel(cms2_metphiInputTag, cms2_metphi_h);

  edm::Handle<float> cms2_sumet_h;
  iEvent.getByLabel(cms2_sumetInputTag, cms2_sumet_h);

  edm::Handle<reco::CaloMETCollection> met_h;
  iEvent.getByLabel(metInputTag, met_h);

  float met      = (met_h->front()).et();
  float met_phi  = (met_h->front()).phi();
  float sumet    = (met_h->front()).sumEt();

  float new_met_x = met * cos(met_phi);
  float new_met_y = met * sin(met_phi);

  float old_met     = *cms2_met_h;
  float old_met_phi = *cms2_metphi_h;
  float old_sumet   = *cms2_sumet_h;
  float old_met_x = (old_met) * cos(old_met_phi);
  float old_met_y = (old_met) * sin(old_met_phi);

  *evt_dmetx   = new_met_x - old_met_x;
  *evt_dmety   = new_met_y - old_met_y;
  *evt_dsumet  = sumet - old_sumet;

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
       branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_dmetx , branchprefix + "dmetx" );
  iEvent.put(evt_dmety , branchprefix + "dmety" );
  iEvent.put(evt_dsumet, branchprefix + "dsumet");
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeltaMETMaker);




  
