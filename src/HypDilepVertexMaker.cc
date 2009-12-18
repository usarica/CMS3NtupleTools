#include "CMS2/NtupleMaker/interface/HypDilepVertexMaker.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "Math/VectorUtil.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;


HypDilepVertexMaker::HypDilepVertexMaker (const edm::ParameterSet& cfg)
{
  recomuonsInputTag_      = cfg.getParameter<edm::InputTag>("recomuonsInputTag"    );
  cms2muonsInputTag_      = cfg.getParameter<edm::InputTag>("cms2muonsInputTag"    );
  recoelectronsInputTag_  = cfg.getParameter<edm::InputTag>("recoelectronsInputTag");
  cms2electronsInputTag_  = cfg.getParameter<edm::InputTag>("cms2electronsInputTag");
  hypInputTag_            = cfg.getParameter<edm::InputTag>("hypInputTag"          );

  produces<vector<int> >          ("hypFVFitstatus"        ).setBranchAlias("hyp_FVFit_status"        );
  produces<vector<int> >          ("hypFVFitndf"           ).setBranchAlias("hyp_FVFit_ndf"           );
  produces<vector<float> >        ("hypFVFitchi2ndf"       ).setBranchAlias("hyp_FVFit_chi2ndf"       );
  produces<vector<float> >        ("hypFVFitprob"          ).setBranchAlias("hyp_FVFit_prob"          );
  produces<vector<float> >        ("hypFVFitv4cxx"         ).setBranchAlias("hyp_FVFit_v4cxx"         );
  produces<vector<float> >        ("hypFVFitv4cxy"         ).setBranchAlias("hyp_FVFit_v4cxy"         );
  produces<vector<float> >        ("hypFVFitv4cyy"         ).setBranchAlias("hyp_FVFit_v4cyy"         );
  produces<vector<float> >        ("hypFVFitv4czz"         ).setBranchAlias("hyp_FVFit_v4czz"         );
  produces<vector<float> >        ("hypFVFitv4czx"         ).setBranchAlias("hyp_FVFit_v4czx"         );
  produces<vector<float> >        ("hypFVFitv4czy"         ).setBranchAlias("hyp_FVFit_v4czy"         );
  produces<vector<LorentzVector> >("hypFVFitp4"            ).setBranchAlias("hyp_FVFit_p4"            );
  produces<vector<LorentzVector> >("hypFVFitv4"            ).setBranchAlias("hyp_FVFit_v4"            );
  
}

void HypDilepVertexMaker::produce(edm::Event& ev, const edm::EventSetup& es){
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder);

  edm::InputTag mus_p4_IT(cms2muonsInputTag_.label(), "musp4");
  Handle<vector<LorentzVector> > mus_p4_h;
  ev.getByLabel(mus_p4_IT, mus_p4_h);
  const vector<LorentzVector>* mus_p4 = mus_p4_h.product(); 

  Handle<vector<reco::Muon> > mus_h;
  ev.getByLabel(recomuonsInputTag_, mus_h);
  vector<reco::Muon> mus = *(mus_h.product());
  MatchUtilities::alignToLVector(*mus_p4, mus);

  edm::InputTag els_p4_IT(cms2electronsInputTag_.label(), "elsp4");
  Handle<vector<LorentzVector> > els_p4_h;
  ev.getByLabel(els_p4_IT, els_p4_h);
  const vector<LorentzVector>* els_p4 = els_p4_h.product(); 

  Handle<vector<reco::GsfElectron> > els_h;
  ev.getByLabel(recoelectronsInputTag_, els_h);
  vector<reco::GsfElectron> els = *(els_h.product());
  MatchUtilities::alignToLVector(*els_p4, els);

  //now read hyp-related stuff
  edm::InputTag hyp_ll_id_IT(hypInputTag_.label(), "hypllid");
  Handle<vector<int> > hyp_ll_id_h;
  ev.getByLabel(hyp_ll_id_IT, hyp_ll_id_h);
  const vector<int>* hyp_ll_id = hyp_ll_id_h.product();

  edm::InputTag hyp_ll_index_IT(hypInputTag_.label(), "hypllindex");
  Handle<vector<int> > hyp_ll_index_h;
  ev.getByLabel(hyp_ll_index_IT, hyp_ll_index_h);
  const vector<int>* hyp_ll_index = hyp_ll_index_h.product();
  
  edm::InputTag hyp_lt_id_IT(hypInputTag_.label(), "hypltid");
  Handle<vector<int> > hyp_lt_id_h;
  ev.getByLabel(hyp_lt_id_IT, hyp_lt_id_h);
  const vector<int>* hyp_lt_id = hyp_lt_id_h.product();

  edm::InputTag hyp_lt_index_IT(hypInputTag_.label(), "hypltindex");
  Handle<vector<int> > hyp_lt_index_h;
  ev.getByLabel(hyp_lt_index_IT, hyp_lt_index_h);
  const vector<int>* hyp_lt_index = hyp_lt_index_h.product();
  
  unsigned int iH = 0;
  unsigned int nH = hyp_ll_id->size();

  auto_ptr<vector<int> >           hyp_FVFit_status (new vector<int>          (nH, -9999)                                 );
  auto_ptr<vector<int> >           hyp_FVFit_ndf    (new vector<int>          (nH, -9999)                                 );
  auto_ptr<vector<float> >         hyp_FVFit_chi2ndf(new vector<float>        (nH, -9999)                                 );
  auto_ptr<vector<float> >         hyp_FVFit_prob   (new vector<float>        (nH, -9999)                                 );
  auto_ptr<vector<float> >         hyp_FVFit_v4cxx  (new vector<float>        (nH, -9999)                                 ); 
  auto_ptr<vector<float> >         hyp_FVFit_v4cxy  (new vector<float>        (nH, -9999)                                 ); 
  auto_ptr<vector<float> >         hyp_FVFit_v4cyy  (new vector<float>        (nH, -9999)                                 ); 
  auto_ptr<vector<float> >         hyp_FVFit_v4czz  (new vector<float>        (nH, -9999)                                 ); 
  auto_ptr<vector<float> >         hyp_FVFit_v4czx  (new vector<float>        (nH, -9999)                                 ); 
  auto_ptr<vector<float> >         hyp_FVFit_v4czy  (new vector<float>        (nH, -9999)                                 ); 
  auto_ptr<vector<LorentzVector> > hyp_FVFit_p4     (new vector<LorentzVector>(nH, LorentzVector(-9999,-9999,-9999,9999)));
  auto_ptr<vector<LorentzVector> > hyp_FVFit_v4     (new vector<LorentzVector>(nH, LorentzVector(-9999,-9999,-9999,9999)));


  for(; iH<nH; ++iH){
    const reco::Track* llTrack = 0;
    const reco::Track* ltTrack = 0;
    float llMass = 0.13957;
    float ltMass = 0.13957;

    int ihyp_ll_id = hyp_ll_id->at(iH);
    int ihyp_lt_id = hyp_lt_id->at(iH);
    int ihyp_ll_index = hyp_ll_index->at(iH);
    int ihyp_lt_index = hyp_lt_index->at(iH);


    if(abs(ihyp_ll_id)==11){
      if(els[ihyp_ll_index].track().isAvailable() && els[ihyp_ll_index].track().isNonnull()){
	llTrack = els[ihyp_ll_index].track().get();
      } else if(els[ihyp_ll_index].gsfTrack().isAvailable() && els[ihyp_ll_index].gsfTrack().isNonnull()){
	llTrack = els[ihyp_ll_index].gsfTrack().get();
      }
      llMass = 0.0005109989;
    }
    if(abs(ihyp_ll_id)==13){
      if(mus[ihyp_ll_index].track().isAvailable() && mus[ihyp_ll_index].track().isNonnull()){
	llTrack = mus[ihyp_ll_index].track().get();
      } //do nothing for muons without a track [might add a case for global mus with inner track coll dropped, not usefull here]
      llMass = 0.105658369;
    }

    //same for lt
    if(abs(ihyp_lt_id)==11){
      if(els[ihyp_lt_index].track().isAvailable() && els[ihyp_lt_index].track().isNonnull()){
	ltTrack = els[ihyp_lt_index].track().get();
      } else if(els[ihyp_lt_index].gsfTrack().isAvailable() && els[ihyp_lt_index].gsfTrack().isNonnull()){
	ltTrack = els[ihyp_lt_index].gsfTrack().get();
      }
      ltMass = 0.0005109989;
    }
    if(abs(ihyp_lt_id)==13){
      if(mus[ihyp_lt_index].track().isAvailable() && mus[ihyp_lt_index].track().isNonnull()){
	ltTrack = mus[ihyp_lt_index].track().get();
      } //do nothing for muons without a track [might add a case for global mus with inner track coll dropped, not usefull here]
      ltMass = 0.105658369;
    }

    if (llTrack == 0 || ltTrack == 0){
      LogDebug("HypDilepVertexMaker")<<"failed to get tracks to vertex";
      (*hyp_FVFit_status)[iH] = 2*(llTrack==0) + 4*(ltTrack==0);
      continue;
    }


    TransientTrack llTT = (*ttrackBuilder).build(llTrack);
    TransientTrack ltTT = (*ttrackBuilder).build(ltTrack);

    //this is a copy-paste from Onia2MuMu.cc
    vector<TransientTrack> t_tks;
    t_tks.push_back(llTT);
    t_tks.push_back(ltTT);
    KalmanVertexFitter kvf;
    TransientVertex tv = kvf.vertex(t_tks);
    
    if (!tv.isValid()){
      LogDebug("HypDilepVertexMaker")<<"failed to fit vertex";
      (*hyp_FVFit_status)[iH] = 8;
      continue;
    }

    (*hyp_FVFit_ndf)[iH]     = (int)tv.degreesOfFreedom();
    (*hyp_FVFit_chi2ndf)[iH] = tv.normalisedChiSquared();
    (*hyp_FVFit_prob)[iH]    = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
    (*hyp_FVFit_v4cxx)[iH]   = tv.positionError().cxx();
    (*hyp_FVFit_v4cxy)[iH]   = tv.positionError().cyx();
    (*hyp_FVFit_v4cyy)[iH]   = tv.positionError().cyy();
    (*hyp_FVFit_v4czz)[iH]   = tv.positionError().czz();
    (*hyp_FVFit_v4czx)[iH]   = tv.positionError().czx();
    (*hyp_FVFit_v4czy)[iH]   = tv.positionError().czy();
    if (tv.hasRefittedTracks()){
      GlobalVector llp3Ref = tv.refittedTracks()[0].trajectoryStateClosestToPoint(tv.position()).momentum();
      GlobalVector ltp3Ref = tv.refittedTracks()[1].trajectoryStateClosestToPoint(tv.position()).momentum();
      float llE = sqrt(llp3Ref.mag2()+llMass*llMass);
      float ltE = sqrt(ltp3Ref.mag2()+ltMass*ltMass);
      //this is not particulalrly numerically stable :(
      LorentzVector llp4Ref(llp3Ref.x(), llp3Ref.y(), llp3Ref.z(), llE);
      LorentzVector ltp4Ref(ltp3Ref.x(), ltp3Ref.y(), ltp3Ref.z(), ltE);
      
      (*hyp_FVFit_p4)[iH]      = llp4Ref+ltp4Ref; 
      (*hyp_FVFit_status)[iH] = 1;
    } else {
      (*hyp_FVFit_status)[iH] = 0;
    }
    (*hyp_FVFit_v4)[iH] = LorentzVector(tv.position().x(), tv.position().y(), tv.position().z(), 0);
  }
	
  ev.put(hyp_FVFit_status , "hypFVFitstatus" );
  ev.put(hyp_FVFit_ndf    , "hypFVFitndf"    );
  ev.put(hyp_FVFit_chi2ndf, "hypFVFitchi2ndf");
  ev.put(hyp_FVFit_prob   , "hypFVFitprob"   );
  ev.put(hyp_FVFit_v4cxx  , "hypFVFitv4cxx"  );
  ev.put(hyp_FVFit_v4cxy  , "hypFVFitv4cxy"  );
  ev.put(hyp_FVFit_v4cyy  , "hypFVFitv4cyy"  );
  ev.put(hyp_FVFit_v4czz  , "hypFVFitv4czz"  );
  ev.put(hyp_FVFit_v4czx  , "hypFVFitv4czx"  );
  ev.put(hyp_FVFit_v4czy  , "hypFVFitv4czy"  );
  ev.put(hyp_FVFit_p4     , "hypFVFitp4"     );
  ev.put(hyp_FVFit_v4     , "hypFVFitv4"     );
}

//define this as a plug-in 
DEFINE_FWK_MODULE(HypDilepVertexMaker);
