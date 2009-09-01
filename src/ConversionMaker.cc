// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ConversionMaker
// 
/**\class ConversionMaker ConversionMaker.cc CMS2/NtupleMaker/src/ConversionMaker.cc

Description: make associations between jets and muons

*/
//
// Original Author:  Puneeth Kalavase 
//         Created:  Wed Oct 15 18:32:24 UTC 2008
// $Id: ConversionMaker.cc,v 1.4 2009/09/01 07:58:43 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/ConversionMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "Math/VectorUtil.h"
#include "Math/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
using namespace std;
using namespace edm;

ConversionMaker::ConversionMaker(const ParameterSet& iConfig)
{ 
  produces<vector<int>   >   ("elsconvtkidx"   ).setBranchAlias("els_conv_tkidx"   );
  produces<vector<float> >   ("elsconvdcot"    ).setBranchAlias("els_conv_dcot"    );
  produces<vector<float> >   ("elsconvdist"    ).setBranchAlias("els_conv_dist" );

  produces<vector<int>   >   ("trksconvtkidx"   ).setBranchAlias("trks_conv_tkidx"   );
  produces<vector<float> >   ("trksconvdcot"    ).setBranchAlias("trks_conv_dcot"    );
  produces<vector<float> >   ("trksconvdist"    ).setBranchAlias("trks_conv_dist" );


  //electronMakerInputTag_ = iConfig.getParameter<InputTag>("electronMakerInputTag"   );
  //tracksMakerInputTag_   = iConfig.getParameter<InputTag>("tracksMakerInputTag"     );
  //bFieldInputTag_        = iConfig.getParameter<InputTag>("bFieldInputTag"          );
     
}

void ConversionMaker::produce(Event& iEvent, const EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  // make vectors to hold the information
  auto_ptr<vector<int>   > els_conv_tkidx        (new vector<int>    );
  auto_ptr<vector<float> > els_conv_dcot         (new vector<float>  );
  auto_ptr<vector<float> > els_conv_dist         (new vector<float>  );

  auto_ptr<vector<int>   > trks_conv_tkidx        (new vector<int>    );
  auto_ptr<vector<float> > trks_conv_dcot         (new vector<float>  );
  auto_ptr<vector<float> > trks_conv_dist         (new vector<float>  );


  // get electron p4s
  Handle<vector<LorentzVector> > els_trk_p4_h;
  iEvent.getByLabel("electronMaker", "elstrkp4", els_trk_p4_h);  
          
  //Handle<vector<float> > els_d0corr_h;
  //iEvent.getByLabel("electronMaker", "elsd0corr", els_d0corr_h);  
     
  Handle<vector<float> > els_d0_h;
  iEvent.getByLabel("electronMaker", "elsd0", els_d0_h);  
     
  Handle<vector<int> > els_charge_h;
  iEvent.getByLabel("electronMaker", "elscharge", els_charge_h);
     
  //get the electron to track association 
  Handle<vector<int> > els_trkidx_h;
  //  iEvent.getByLabel("elToTrackAssMaker", "elstrkidx", els_trkidx_h);
  iEvent.getByLabel("electronMaker", "elstrkidx", els_trkidx_h);
     
  //now get the Track quantities
  Handle<vector<LorentzVector> > trks_p4_h;
  iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);  
          
  Handle<vector<float> > trks_d0corr_h;
  iEvent.getByLabel("trackMaker", "trksd0corr", trks_d0corr_h);  

  Handle<vector<float> > trks_d0_h;
  iEvent.getByLabel("trackMaker", "trksd0", trks_d0_h);  
  
  Handle<vector<int> > trks_charge_h;
  iEvent.getByLabel("trackMaker", "trkscharge", trks_charge_h);

  Handle<float> evt_bField_h;
  iEvent.getByLabel("eventMaker", "evtbField", evt_bField_h);
  
  for(unsigned int elIndex = 0; elIndex < els_trk_p4_h->size();
      elIndex++) { //electron loop
       
    int closestTkIdx = -999;
    float mindEta = 999;
    //pair.first = dist, pair.second = dcot
    pair <float, float> p_elConvInfo = make_pair(-999,-999);  
    for(unsigned int tkIndex = 0; tkIndex < trks_p4_h->size();
	tkIndex++) {
      
      if(els_charge_h->at(elIndex) + trks_charge_h->at(tkIndex) != 0) continue;
      //don't include the electron's track
      
      if(els_trkidx_h->at(elIndex) == (int)tkIndex) continue;
      //look only in cone of 0.5, which is probably a little big, but whatever
      if(ROOT::Math::VectorUtil::DeltaR(els_trk_p4_h->at(elIndex), trks_p4_h->at(tkIndex)) > 0.5) 
	continue;
      if(fabs(els_trk_p4_h->at(elIndex).Eta() - trks_p4_h->at(tkIndex).Eta()) < mindEta) {
	mindEta = fabs(els_trk_p4_h->at(elIndex).Eta() - trks_p4_h->at(tkIndex).Eta());
	closestTkIdx = tkIndex;
	p_elConvInfo = getConversionInfo(els_trk_p4_h->at(elIndex), els_charge_h->at(elIndex),
					 els_d0_h->at(elIndex), trks_p4_h->at(tkIndex),
					 trks_charge_h->at(tkIndex), trks_d0_h->at(tkIndex),
					 *evt_bField_h.product());
      }//if dEta < minDeta
    }//for(unsigned int tkIndex = 0.....
    els_conv_tkidx    ->push_back(closestTkIdx            );
    els_conv_dist     ->push_back(p_elConvInfo.first          );
    els_conv_dcot     ->push_back(p_elConvInfo.second         );
    
  }//for(unsigned int elIndex = 0........


  //now do the same for the tracks
  for(unsigned int tk1Index = 0; tk1Index < trks_p4_h->size();
      tk1Index++) {
       
    int closestTkIdx = -999;
    float mindEta = 999;
    //pair.first=dist, pair.second = dcot
    pair <float, float> p_tkConvInfo = make_pair(-999,-999);
    for(unsigned int tk2Index = 0; tk2Index < trks_p4_h->size();
	tk2Index++) {
	 
      if(trks_charge_h->at(tk1Index) + trks_charge_h->at(tk2Index) != 0 ) continue;
      if(tk1Index == tk2Index) continue; //don't want to use the same track!
      //look only in cone of 0.5, which is probably a little big, but whatever
      if(ROOT::Math::VectorUtil::DeltaR(trks_p4_h->at(tk1Index), trks_p4_h->at(tk2Index)) > 0.5) 
	continue;
      
      if(fabs(trks_p4_h->at(tk1Index).Eta() - trks_p4_h->at(tk2Index).Eta()) < mindEta) {
	mindEta = fabs(trks_p4_h->at(tk1Index).Eta() - trks_p4_h->at(tk2Index).Eta());
	closestTkIdx = tk2Index;
	p_tkConvInfo = getConversionInfo(trks_p4_h->at(tk1Index), trks_charge_h->at(tk1Index),
					 trks_d0_h->at(tk1Index), trks_p4_h->at(tk2Index),
					 trks_charge_h->at(tk2Index), trks_d0_h->at(tk2Index),
					 *evt_bField_h.product());
      }//if dEta < mindEta
    }//tk2 loop
    trks_conv_tkidx    ->push_back(closestTkIdx            );
    trks_conv_dist     ->push_back(p_tkConvInfo.first          );
    trks_conv_dcot     ->push_back(p_tkConvInfo.second         );
  }//tk1 loop

  // store vectors
  iEvent.put(els_conv_tkidx,        "elsconvtkidx"   );
  iEvent.put(els_conv_dcot,         "elsconvdcot"    );
  iEvent.put(els_conv_dist,         "elsconvdist"    );
  iEvent.put(trks_conv_tkidx,       "trksconvtkidx"  );
  iEvent.put(trks_conv_dist,        "trksconvdist"   );
  iEvent.put(trks_conv_dcot,        "trksconvdcot"   );

}
		     


//-------------------------------------------------------------------------------------
//Function that does the work of finding the dist and dcot theta 
//--------------------------

pair<float, float> ConversionMaker::getConversionInfo(math::XYZTLorentzVector trk1_p4, 
						      int trk1_q, float trk1_d0, 
						      math::XYZTLorentzVector trk2_p4,
						      int trk2_q, float trk2_d0,
						      float bField) {
  
  
  double tk1Curvature = -0.3*bField*(trk1_q/trk1_p4.pt())/100.;
  double rTk1 = fabs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_p4.phi());
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_p4.phi());
  
  
  
  
  double tk2Curvature = -0.3*bField*(trk2_q/trk2_p4.pt())/100.;
  double rTk2 = fabs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_p4.phi());
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_p4.phi());
	 
  double dist = sqrt(pow(xTk1-xTk2, 2) + pow(yTk1-yTk2 , 2));
  dist = dist - (rTk1 + rTk2);

  double dcot = 1/tan(trk1_p4.theta()) - 1/tan(trk2_p4.theta());

  return make_pair(dist, dcot);
  
}


// ------------ method called once each job just before starting event loop  ------------
void ConversionMaker::beginJob(const EventSetup& es)
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ConversionMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ConversionMaker);


