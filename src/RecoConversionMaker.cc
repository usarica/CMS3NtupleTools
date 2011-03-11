// -*- C++ -*-
//
// Package:    RecoConversionMaker
// Class:      RecoConversionMaker
// 
/**\class RecoConversionMaker RecoConversionMaker.cc CMS2/RecoConversionMaker/src/RecoConversionMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: RecoConversionMaker.cc,v 1.1 2011/03/11 02:08:47 kalavase Exp $
//
//

// system include files
#include <memory>
#include "Math/VectorUtil.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/RecoConversionMaker.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
// Propagator specific include files

using namespace std;
using namespace reco;
using namespace edm;
typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPointF Point;


//
// class decleration
//

//
// constructors and destructor
//
RecoConversionMaker::RecoConversionMaker(const edm::ParameterSet& iConfig) {
       
  recoConversionInputTag_ = iConfig.getParameter<edm::InputTag>("recoConversionInputTag");
  beamSpotInputTag_       = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");


  produces<vector<int> >		("convsisConverted"	).setBranchAlias("convs_isConverted"	);
  produces<vector<int> >		("convsalgo"		).setBranchAlias("convs_algo"		);
  produces<vector<vector<int> > >	("convstkidx"		).setBranchAlias("convs_tkidx"		);
  produces<vector<vector<int> > >	("convsnHitsBeforeVtx"	).setBranchAlias("convs_nHitsBeforeVtx"	);

  produces<vector<float> >              ("convsfitProb"         ).setBranchAlias("convs_fitProb"        );
  produces<vector<float> >              ("convsdl"              ).setBranchAlias("convs_dl"             );

  //this comes as a momentum vector from the Conversion object, not a p4 so the 
  produces<vector<LorentzVector> >      ("convsrefitPairMomp4"  ).setBranchAlias("convs_refitPairMom_p4");
  produces<vector<LorentzVector> >      ("convsvtxpos"          ).setBranchAlias("convs_vtxpos"         );

  //produces<vector<vector<float> > >     ("convsdlClosestHit"    ).setBranchAlias("convs_dlClosestHit"   );// signed decay length between nearest hit and conversion vertex
  //produces<vector<vector<float> > >     ("convsdlClosestHitErr" ).setBranchAlias("convs_dlClosestHitErr");// signed decay length error

}

void RecoConversionMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  auto_ptr<vector<int> >		convs_isConverted	(new vector<int>		);
  auto_ptr<vector<int> >		convs_algo		(new vector<int>		);
  auto_ptr<vector<vector<int> > >       convs_tkidx		(new vector<vector<int> >       );
  auto_ptr<vector<vector<int> > >       convs_nHitsBeforeVtx	(new vector<vector<int> >       );
  
  auto_ptr<vector<float> >              convs_fitProb           (new vector<float>              );          
  auto_ptr<vector<float> >              convs_dl                (new vector<float>              );
  auto_ptr<vector<LorentzVector> >      convs_refitPairMom_p4   (new vector<LorentzVector>      );
  auto_ptr<vector<LorentzVector> >      convs_vtxpos            (new vector<LorentzVector>      );
  //auto_ptr<vector<vector<float> > >     convs_dlClosestHit      (new vector<vector<float> >     );
  //auto_ptr<vector<vector<float> > >     convs_dlClosestHitErr   (new vector<vector<float> >     );
	   

  // get reco Conversions
  Handle<View<Conversion> > convs_h;
  iEvent.getByLabel(recoConversionInputTag_, convs_h);
  
  Handle<BeamSpot> beamSpotH;
  iEvent.getByLabel(beamSpotInputTag_, beamSpotH);

  for(View<Conversion>::const_iterator it = convs_h->begin();
      it != convs_h->end(); it++) { 

    convs_isConverted->push_back(it->isConverted());
    convs_algo->push_back(it->algo());
    vector<edm::RefToBase<reco::Track> > v_temp_trks = it->tracks();
    vector<int> v_temp_out;
    for(unsigned int i = 0; i < v_temp_trks.size(); i++) 
      v_temp_out.push_back(v_temp_trks.at(i).key());
    convs_tkidx->push_back(v_temp_out);

    v_temp_out.clear();
    vector<uint8_t> v_temp_nhits = it->nHitsBeforeVtx();
    for(unsigned int i = 0; i < v_temp_nhits.size(); i++) 
      v_temp_out.push_back(v_temp_nhits.at(i));
    convs_nHitsBeforeVtx->push_back(v_temp_out);

    convs_fitProb->push_back(ChiSquaredProbability( it->conversionVertex().chi2(),  it->conversionVertex().ndof() ));
    convs_dl->push_back(lxy(beamSpotH->position(), *it));

    convs_refitPairMom_p4->push_back(LorentzVector(it->refittedPair4Momentum()));
    convs_vtxpos->push_back(LorentzVector(it->conversionVertex().x(), 
					  it->conversionVertex().y(),
					  it->conversionVertex().z(),
					  0));
    
    /*
      std::vector<Measurement1DFloat> v_temp_dlClosestHit = it->dlClosestHitToVtx();
      vector<float> v_fltemp, v_fltemp_err;
      for(unsigned int i = 0; i < v_temp_dlClosestHit.size(); i++) {
      v_fltemp.push_back(v_temp_dlClosestHit.at(i).value());
      v_fltemp_err.push_back(v_temp_dlClosestHit.at(i).error());
      }
      convs_dlClosestHit->push_back(v_fltemp);
      convs_dlClosestHitErr->push_back(v_fltemp_err);
    */
      

      }//reco conversion loop



  iEvent.put(convs_isConverted		, "convsisConverted"	 );
  iEvent.put(convs_algo			, "convsalgo"		 );
  iEvent.put(convs_tkidx		, "convstkidx"		 );
  iEvent.put(convs_nHitsBeforeVtx	, "convsnHitsBeforeVtx"	 );

  iEvent.put(convs_fitProb              , "convsfitProb"         );
  iEvent.put(convs_dl                   , "convsdl"              );
  iEvent.put(convs_refitPairMom_p4      , "convsrefitPairMomp4"  );
  iEvent.put(convs_vtxpos               , "convsvtxpos"          );
  //iEvent.put(convs_dlClosestHit       , "convsdlClosestHit"    );
  //iEvent.put(convs_dlClosestHitErr	, "convsdlClosestHitErr" );
  
}


double RecoConversionMaker::lxy(const math::XYZPoint& myBeamSpot, const Conversion& conv) const {

  const reco::Vertex &vtx = conv.conversionVertex();
  if (!vtx.isValid()) return -9999.;

  math::XYZVector mom = conv.refittedPairMomentum();
  
  double dbsx = vtx.x() - myBeamSpot.x();
  double dbsy = vtx.y() - myBeamSpot.y();
  double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
  return lxy;  
  
}

// ------------ method called once each job just before starting event loop  ------------
void RecoConversionMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void RecoConversionMaker::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoConversionMaker);
