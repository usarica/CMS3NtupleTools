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
// $Id: ConversionMaker.cc,v 1.7 2009/11/18 21:53:03 kalavase Exp $
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
#include "TMath.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace std;
using namespace edm;

ConversionMaker::ConversionMaker(const ParameterSet& iConfig)
{ 
  produces<vector<int>   >   ("trksconvtkidx"   ).setBranchAlias("trks_conv_tkidx" );
  produces<vector<float> >   ("trksconvdcot"    ).setBranchAlias("trks_conv_dcot"  );
  produces<vector<float> >   ("trksconvdist"    ).setBranchAlias("trks_conv_dist"  );

  minFracShHits_ = iConfig.getParameter<double>("minFracSharedHits");
     
}

void ConversionMaker::produce(Event& iEvent, const EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  auto_ptr<vector<int>   > trks_conv_tkidx       (new vector<int>    );
  auto_ptr<vector<float> > trks_conv_dcot        (new vector<float>  );
  auto_ptr<vector<float> > trks_conv_dist        (new vector<float>  );


    
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
  float bField = *evt_bField_h.product();
  

  //now do the same for the tracks
  for(unsigned int tk1Index = 0; tk1Index < trks_p4_h->size();
      tk1Index++) {
       
    int closestTkIdx = -9999;
    float mindR = 9999.;
    //pair.first=dist, pair.second = dcot
    pair <float, float> p_tkConvInfo = make_pair(-9999.,-9999.);
    for(unsigned int tk2Index = 0; tk2Index < trks_p4_h->size();
	tk2Index++) {
	 
      if(trks_charge_h->at(tk1Index) + trks_charge_h->at(tk2Index) != 0 ) continue;
      if(tk1Index == tk2Index) continue; //don't want to use the same track!
      //look only in cone of 0.3
      double dR = deltaR( trks_p4_h->at(tk1Index).eta(), trks_p4_h->at(tk1Index).phi(),
			  trks_p4_h->at(tk2Index).eta(), trks_p4_h->at(tk2Index).phi());
      
      if(dR  > 0.3)
	continue;
      
      if(dR < mindR) {
	mindR = dR;
	closestTkIdx = tk2Index;
	p_tkConvInfo = ConversionFinder::getConversionInfo(math::XYZTLorentzVector(trks_p4_h->at(tk1Index)), trks_charge_h->at(tk1Index),
							   trks_d0_h->at(tk1Index), math::XYZTLorentzVector(trks_p4_h->at(tk2Index)),
							   trks_charge_h->at(tk2Index), trks_d0_h->at(tk2Index),
							   bField);
      }//if(dR < mindR)
    }//tk2 loop

    trks_conv_tkidx    ->push_back(closestTkIdx                );
    trks_conv_dist     ->push_back(p_tkConvInfo.first          );
    trks_conv_dcot     ->push_back(p_tkConvInfo.second         );
  }//tk1 loop

  // store vectors
  iEvent.put(trks_conv_tkidx,       "trksconvtkidx"  );
  iEvent.put(trks_conv_dist,        "trksconvdist"   );
  iEvent.put(trks_conv_dcot,        "trksconvdcot"   );

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


