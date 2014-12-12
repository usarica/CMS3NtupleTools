//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      BeamSpotMaker
// 
/**\class BeamSpotMaker BeamSpotMaker.cc CMS2/NtupleMakerMaker/src/BeamSpotMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: BeamSpotMaker.cc,v 1.10 2010/03/18 02:11:48 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/BeamSpotMaker.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

BeamSpotMaker::BeamSpotMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  //p4 because we're not able to (yet) read XYZPointDs in bare root for some reason 
  //the 4th co-ordinate is 0
  produces<LorentzVector>   (branchprefix+"p4"         ).setBranchAlias(aliasprefix_+"p4"            );
  produces<int>             (branchprefix+"Type"       ).setBranchAlias(aliasprefix_+"Type"          );
  produces<float>           (branchprefix+"xErr"       ).setBranchAlias(aliasprefix_+"_xErr"         );
  produces<float>           (branchprefix+"yErr"       ).setBranchAlias(aliasprefix_+"_yErr"         );
  produces<float>           (branchprefix+"zErr"       ).setBranchAlias(aliasprefix_+"_zErr"         );
  produces<float>           (branchprefix+"sigmaZ"     ).setBranchAlias(aliasprefix_+"_sigmaZ"       );
  produces<float>           (branchprefix+"sigmaZErr"  ).setBranchAlias(aliasprefix_+"_sigmaZErr"    );
  produces<float>           (branchprefix+"dxdz"       ).setBranchAlias(aliasprefix_+"_dxdz"         );
  produces<float>           (branchprefix+"dxdzErr"    ).setBranchAlias(aliasprefix_+"_dxdzErr"      );
  produces<float>           (branchprefix+"dydz"       ).setBranchAlias(aliasprefix_+"_dydz"         );
  produces<float>           (branchprefix+"dydzErr"    ).setBranchAlias(aliasprefix_+"_dydzErr"      );
  produces<float>           (branchprefix+"Xwidth"     ).setBranchAlias(aliasprefix_+"_Xwidth"       );
  produces<float>           (branchprefix+"Ywidth"     ).setBranchAlias(aliasprefix_+"_Ywidth"       );
  produces<float>           (branchprefix+"XwidthErr"  ).setBranchAlias(aliasprefix_+"_XwidthErr"    );
  produces<float>           (branchprefix+"YwidthErr"  ).setBranchAlias(aliasprefix_+"_YwidthErr"    );
  produces<vector <float> > (branchprefix+"covMatrix"  ).setBranchAlias(aliasprefix_+"_covMatrix"    ); 
  
  beamSpotInputTag = iConfig.getParameter<InputTag>("beamSpotInputTag");
  
}


BeamSpotMaker::~BeamSpotMaker() {}

void  BeamSpotMaker::beginJob() {
}

void BeamSpotMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BeamSpotMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<LorentzVector>  evt_bs_p4          (new LorentzVector);
  auto_ptr<int>            evt_bsType         (new int          );
  auto_ptr<float>          evt_bs_xErr        (new float        );
  auto_ptr<float>          evt_bs_yErr        (new float        );
  auto_ptr<float>          evt_bs_zErr        (new float        );
  auto_ptr<float>          evt_bs_sigmaZ      (new float        );
  auto_ptr<float>          evt_bs_sigmaZErr   (new float        );
  auto_ptr<float>          evt_bs_dxdz        (new float        );
  auto_ptr<float>          evt_bs_dxdzErr     (new float        );
  auto_ptr<float>          evt_bs_dydz        (new float        );
  auto_ptr<float>          evt_bs_dydzErr     (new float        );
  auto_ptr<float>          evt_bs_Xwidth      (new float        );
  auto_ptr<float>          evt_bs_Ywidth      (new float        );
  auto_ptr<float>          evt_bs_XwidthErr   (new float        );
  auto_ptr<float>          evt_bs_YwidthErr   (new float        );
  auto_ptr<vector<float> > evt_bs_covMatrix   (new vector<float>);
  
  Handle<BeamSpot> beamSpotH;
  iEvent.getByLabel(beamSpotInputTag, beamSpotH);

  bool haveBeamSpot = true;
  if(!beamSpotH.isValid() )
    haveBeamSpot = false;
  
  *evt_bs_p4         = haveBeamSpot ? LorentzVector(beamSpotH->position().x(),
						    beamSpotH->position().y(),
						    beamSpotH->position().z(),
						    0) : LorentzVector(0,0,0,0);

  *evt_bsType        = haveBeamSpot ? beamSpotH->type()           : -9999;
  *evt_bs_xErr       = haveBeamSpot ? beamSpotH->x0Error()        : 0.0;
  *evt_bs_yErr       = haveBeamSpot ? beamSpotH->y0Error()        : 0.0;
  *evt_bs_zErr       = haveBeamSpot ? beamSpotH->z0Error()        : 0.0;
  *evt_bs_sigmaZ     = haveBeamSpot ? beamSpotH->sigmaZ()         : 0.0;
  *evt_bs_sigmaZErr  = haveBeamSpot ? beamSpotH->sigmaZ0Error()   : 0.0;
  *evt_bs_dxdz       = haveBeamSpot ? beamSpotH->dxdz()           : 0.0;
  *evt_bs_dxdzErr    = haveBeamSpot ? beamSpotH->dxdzError()      : 0.0;
  *evt_bs_dydz       = haveBeamSpot ? beamSpotH->dydz()           : 0.0;	
  *evt_bs_dydzErr    = haveBeamSpot ? beamSpotH->dydzError()      : 0.0;
  *evt_bs_Xwidth     = haveBeamSpot ? beamSpotH->BeamWidthX()     : 0.0;
  *evt_bs_Ywidth     = haveBeamSpot ? beamSpotH->BeamWidthY()     : 0.0;
  *evt_bs_XwidthErr  = haveBeamSpot ? beamSpotH->BeamWidthXError(): 0.0;
  *evt_bs_YwidthErr  = haveBeamSpot ? beamSpotH->BeamWidthYError(): 0.0;

  const unsigned int covMatrix_dim = 7;

  for( unsigned int i = 0; i < covMatrix_dim; i++ ) {
    for( unsigned int j = 0; j < covMatrix_dim; j++ ) {
      evt_bs_covMatrix->push_back( beamSpotH->covariance(i, j) );
    }
  }  


  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_bs_p4           , branchprefix+"p4"        );
  iEvent.put(evt_bsType          , branchprefix+"Type"      );

  if(haveBeamSpot) {
    iEvent.put(evt_bs_xErr       , branchprefix+"xErr"      );
    iEvent.put(evt_bs_yErr       , branchprefix+"yErr"      );
    iEvent.put(evt_bs_zErr       , branchprefix+"zErr"      );
    iEvent.put(evt_bs_sigmaZ     , branchprefix+"sigmaZ"    );
    iEvent.put(evt_bs_sigmaZErr  , branchprefix+"sigmaZErr" );
    iEvent.put(evt_bs_dxdz       , branchprefix+"dxdz"      );
    iEvent.put(evt_bs_dxdzErr    , branchprefix+"dxdzErr"   );
    iEvent.put(evt_bs_dydz       , branchprefix+"dydz"      );
    iEvent.put(evt_bs_dydzErr    , branchprefix+"dydzErr"   );
    iEvent.put(evt_bs_Xwidth     , branchprefix+"Xwidth"    );
    iEvent.put(evt_bs_Ywidth     , branchprefix+"Ywidth"    );
    iEvent.put(evt_bs_XwidthErr  , branchprefix+"XwidthErr" );
    iEvent.put(evt_bs_YwidthErr  , branchprefix+"YwidthErr" );
    iEvent.put(evt_bs_covMatrix  , branchprefix+"covMatrix" );
  }
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(BeamSpotMaker);
