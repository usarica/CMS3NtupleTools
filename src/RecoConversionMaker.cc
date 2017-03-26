// -*- C++ -*-

#include <memory>

#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "CMS3/NtupleMaker/interface/RecoConversionMaker.h"

using namespace std;
using namespace reco;
using namespace edm;

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPointF Point;

RecoConversionMaker::RecoConversionMaker(const edm::ParameterSet& iConfig) {
       
  recoConversionToken = consumes<edm::View<Conversion> >(iConfig.getParameter<edm::InputTag>("recoConversionInputTag"));
  beamSpotToken = consumes<BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotInputTag"));


  produces<vector<int> >		("convsisConverted"	).setBranchAlias("convs_isConverted"	);
  produces<vector<int> >                ("convsquality"         ).setBranchAlias("convs_quality"        );
  produces<vector<int> >		("convsalgo"		).setBranchAlias("convs_algo"		);
  produces<vector<vector<int> > >	("convstkidx"		).setBranchAlias("convs_tkidx"		);
  produces<vector<vector<int> > >	("convstkalgo"		).setBranchAlias("convs_tkalgo"		);
  produces<vector<vector<int> > >	("convsnHitsBeforeVtx"	).setBranchAlias("convs_nHitsBeforeVtx"	);

  produces<vector<float> >              ("convsndof"            ).setBranchAlias("convs_ndof"           );
  produces<vector<float> >              ("convschi2"            ).setBranchAlias("convs_chi2"           );
  produces<vector<float> >              ("convsdl"              ).setBranchAlias("convs_dl"             );

}

void RecoConversionMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  auto_ptr<vector<int> >		convs_isConverted	(new vector<int>		);
  auto_ptr<vector<int> >		convs_algo		(new vector<int>		);
  auto_ptr<vector<int> >                convs_quality           (new vector<int>                );
  
  auto_ptr<vector<vector<int> > >       convs_tkidx		(new vector<vector<int> >       );
  auto_ptr<vector<vector<int> > >       convs_tkalgo            (new vector<vector<int> >       );
  auto_ptr<vector<vector<int> > >       convs_nHitsBeforeVtx	(new vector<vector<int> >       );
  
  auto_ptr<vector<float> >              convs_ndof              (new vector<float>              );          
  auto_ptr<vector<float> >              convs_chi2              (new vector<float>              );            
  auto_ptr<vector<float> >              convs_dl                (new vector<float>              );

  // get reco Conversions
  Handle<View<Conversion> > convs_h;
  iEvent.getByToken(recoConversionToken, convs_h);
  
  Handle<BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken, beamSpotH);

  for(View<Conversion>::const_iterator it = convs_h->begin(); it != convs_h->end(); it++){ 

    convs_isConverted->push_back(it->isConverted());
    convs_algo->push_back(it->algo());
    //quality 
    int qualityMask = 0;
    for(int iM = 0; iM < 32; ++iM){
      if(it->quality((Conversion::ConversionQuality)iM)) qualityMask |= 1 << iM;
    }
    convs_quality->push_back(qualityMask);
    
    vector<int> v_temp_out;
    vector<uint8_t> v_temp_nhits = it->nHitsBeforeVtx();
    for(unsigned int i = 0; i < v_temp_nhits.size(); i++) v_temp_out.push_back(v_temp_nhits.at(i));
    convs_nHitsBeforeVtx->push_back(v_temp_out);

    convs_ndof->push_back(it->conversionVertex().ndof() );
    convs_chi2->push_back(it->conversionVertex().chi2() );
    convs_dl->push_back(lxy(beamSpotH->position(), *it));

  }//reco conversion loop

  iEvent.put(convs_isConverted		, "convsisConverted"	 );
  iEvent.put(convs_quality          , "convsquality"         );
  iEvent.put(convs_algo		      	, "convsalgo"		     );
  iEvent.put(convs_tkidx	      	, "convstkidx"		     );
  iEvent.put(convs_tkalgo	      	, "convstkalgo"		     );
  iEvent.put(convs_nHitsBeforeVtx	, "convsnHitsBeforeVtx"	 );
  iEvent.put(convs_ndof             , "convsndof"            );
  iEvent.put(convs_chi2             , "convschi2"            );
  iEvent.put(convs_dl               , "convsdl"              );
  
}

double RecoConversionMaker::lxy(const math::XYZPoint& myBeamSpot, const Conversion& conv) const {

  const reco::Vertex &vtx = conv.conversionVertex();
  if (!vtx.isValid()) return -9999.;

  math::XYZVectorF mom = conv.refittedPairMomentum();
  
  double dbsx = vtx.x() - myBeamSpot.x();
  double dbsy = vtx.y() - myBeamSpot.y();
  double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
  return lxy;  
  
}

// ------------ method called once each job just before starting event loop  ------------
void RecoConversionMaker::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void RecoConversionMaker::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoConversionMaker);
