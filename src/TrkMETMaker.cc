//-*- C++ -*-
//
// Package:    TrkMETMaker
// Class:      TrkMETMaker
// 
/**\class TrkMETMaker TrkMETMaker.cc CMS2/TrkMETMaker/src/TrkMETMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrkMETMaker.cc,v 1.2 2012/05/14 17:15:17 cerati Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CMS3/NtupleMaker/interface/TrkMETMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
TrkMETMaker::TrkMETMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
    
  pfcandInputTag_           = iConfig.getParameter<InputTag>("pfcandInputTag"                          );
  hypInputTag_              = iConfig.getParameter<InputTag>("hypInputTag"                             );
  vertexInputTag_           = iConfig.getParameter<InputTag>("vertexInputTag"                          );
  
  dzcut_                     = iConfig.getParameter<double>          ("dzcut");
  drcut_                     = iConfig.getParameter<double>          ("drcut");
  correctJets_               = iConfig.getParameter<bool>            ("correctJet");

  produces<vector<float> >           (branchprefix+"ltdPhiTrkMet"               ).setBranchAlias(aliasprefix_+"_lt_dPhi_TrkMet"                );
  produces<vector<float> >           (branchprefix+"lldPhiTrkMet"               ).setBranchAlias(aliasprefix_+"_ll_dPhi_TrkMet"                );

  produces<vector<float> >           (branchprefix+"met"             ).setBranchAlias(aliasprefix_+"_met"     );
  produces<vector<float> >           (branchprefix+"metPhi"          ).setBranchAlias(aliasprefix_+"_metPhi"  );
  produces<vector<float> >           (branchprefix+"sumet"           ).setBranchAlias(aliasprefix_+"_sumet"   );

}

TrkMETMaker::~TrkMETMaker() {

}


void  TrkMETMaker::beginRun(const edm::Run&, const edm::EventSetup& es) {
     

}

void TrkMETMaker::beginJob() {


}

void TrkMETMaker::endJob() {
}

double TrkMETMaker::dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

// ------------ method called to produce the data  ------------
void TrkMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<float> >		hyp_trkmet		   (new vector<float>		);
  auto_ptr<vector<float> >		hyp_trkmetPhi		   (new vector<float>		);
  auto_ptr<vector<float> >		hyp_trksumet		   (new vector<float>		);

  auto_ptr<vector<float> >		hyp_ltdPhiTrkMet	   (new vector<float>		);
  auto_ptr<vector<float> >		hyp_lldPhiTrkMet           (new vector<float>		);
    
  //pfcandidate p4
  InputTag pfcands_p4_tag(pfcandInputTag_.label(),"pfcandsp4");
  Handle<vector<LorentzVector> > pfcands_p4_h;
  iEvent.getByLabel(pfcands_p4_tag, pfcands_p4_h);
  const vector<LorentzVector> *pfcands_p4 = pfcands_p4_h.product();

  //pfcandidate charge
  InputTag pfcands_charge_tag(pfcandInputTag_.label(),"pfcandscharge");
  Handle<vector<int> > pfcands_charge_h;
  iEvent.getByLabel(pfcands_charge_tag, pfcands_charge_h);
  const vector<int> *pfcands_charge = pfcands_charge_h.product();

  //pfcandidate dz
  InputTag pfcands_dz_tag(pfcandInputTag_.label(),"pfcandsdz");
  Handle<vector<float> > pfcands_dz_h;
  iEvent.getByLabel(pfcands_dz_tag, pfcands_dz_h);
  const vector<float> *pfcands_dz = pfcands_dz_h.product();

  //hyp ll p4
  InputTag hyp_ll_p4_tag(hypInputTag_.label(),"hypllp4");
  Handle<vector<LorentzVector> > hyp_ll_p4_h;
  iEvent.getByLabel(hyp_ll_p4_tag, hyp_ll_p4_h);
  const vector<LorentzVector> *hyp_ll_p4 = hyp_ll_p4_h.product();

  //hyp lt p4
  InputTag hyp_lt_p4_tag(hypInputTag_.label(),"hypltp4");
  Handle<vector<LorentzVector> > hyp_lt_p4_h;
  iEvent.getByLabel(hyp_lt_p4_tag, hyp_lt_p4_h);
  const vector<LorentzVector> *hyp_lt_p4 = hyp_lt_p4_h.product();

  //hyp jets p4
  InputTag hyp_jets_p4_tag(hypInputTag_.label(),"hypjetsp4");
  Handle<vector<vector<LorentzVector> > > hyp_jets_p4_h;
  iEvent.getByLabel(hyp_jets_p4_tag, hyp_jets_p4_h);
  const vector<vector<LorentzVector> > *hyp_jets_p4 = hyp_jets_p4_h.product();

  const unsigned int npfcands = pfcands_p4->size();
  const unsigned int nhyps    = hyp_ll_p4->size();

  //-----------------------------------
  // loop over hypotheses
  //-----------------------------------

  for( unsigned int ihyp = 0 ; ihyp < nhyps ; ihyp++ ){

    float metx  = 0;
    float mety  = 0;
    float sumet = 0;
    
    //------------------------------
    // correct met for hyp leptons
    //------------------------------

    metx -= hyp_ll_p4->at(ihyp).Px();
    metx -= hyp_lt_p4->at(ihyp).Px();
    mety -= hyp_ll_p4->at(ihyp).Py();
    mety -= hyp_lt_p4->at(ihyp).Py();

    //------------------------------
    // correct met for hyp jets
    //------------------------------

    if( correctJets_ ){
      for( unsigned int ijet = 0 ; ijet < hyp_jets_p4->at(ihyp).size() ; ijet++ ){

	LorentzVector jet = (hyp_jets_p4->at(ihyp)).at(ijet);

	metx -= jet.Px();
	mety -= jet.Py();
      }
    }

    //-----------------------------------
    // loop over pfcandidates
    //-----------------------------------

    for( unsigned int ipf = 0 ; ipf < npfcands ; ipf++ ){

      // require charged pfcandidate
      if( pfcands_charge->at(ipf) == 0 ) continue;

      // don't correct for pfcandidates dr-matched to hyp leptons
      double dRll = deltaR( hyp_ll_p4->at(ihyp).eta() , hyp_ll_p4->at(ihyp).phi() , pfcands_p4->at(ipf).eta() , pfcands_p4->at(ipf).phi());
      double dRlt = deltaR( hyp_lt_p4->at(ihyp).eta() , hyp_lt_p4->at(ihyp).phi() , pfcands_p4->at(ipf).eta() , pfcands_p4->at(ipf).phi());
      if( dRll < drcut_ || dRlt < drcut_ ) continue;

      // now make dz requirement on pfcandidate
      if( pfcands_dz->at(ipf) > dzcut_ ) continue;

      // skip pfcandidates matched to jet
      if( correctJets_ ){

       	bool skip = false;

 	for( unsigned int ijet = 0 ; ijet < hyp_jets_p4->at(ihyp).size() ; ijet++ ){

	  LorentzVector jet = (hyp_jets_p4->at(ihyp)).at(ijet);
	  double dRjet = deltaR( jet.eta() , jet.phi() , pfcands_p4->at(ipf).eta() , pfcands_p4->at(ipf).phi());
       	  if( dRjet < 0.5 ) skip = true;
       	}

	if( skip ) continue;
      }

      // pfcandidate passes selection so correct the met
      metx  -= pfcands_p4->at(ipf).Px();
      mety  -= pfcands_p4->at(ipf).Py();
      sumet += pfcands_p4->at(ipf).Pt();

    }//pfcandidates

    hyp_lldPhiTrkMet->push_back ( deltaPhi(hyp_ll_p4->at(ihyp).phi(),atan2(mety,metx)) );
    hyp_ltdPhiTrkMet->push_back ( deltaPhi(hyp_lt_p4->at(ihyp).phi(),atan2(mety,metx)) );

    hyp_trkmet->push_back       ( sqrt(metx*metx+mety*mety)                            );
    hyp_trkmetPhi->push_back    ( atan2(mety,metx)                                     );
    hyp_trksumet->push_back     ( sumet                                                );
    
  }//hypotheses

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(hyp_trkmet			, branchprefix+"met"		    );
  iEvent.put(hyp_trkmetPhi	       	, branchprefix+"metPhi"		    );
  iEvent.put(hyp_trksumet	       	, branchprefix+"sumet"		    );
  iEvent.put(hyp_lldPhiTrkMet	       	, branchprefix+"lldPhiTrkMet"	    );
  iEvent.put(hyp_ltdPhiTrkMet	       	, branchprefix+"ltdPhiTrkMet"	    );
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkMETMaker);
