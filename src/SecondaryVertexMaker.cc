// -*- C++ -*-
//
// Package:    SecondaryVertexMaker
// Class:      SecondaryVertexMaker
// 
/**\class SecondaryVertexMaker SecondaryVertexMaker.cc CMS2/SecondaryVertexMaker/src/SecondaryVertexMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TMath.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class decleration
//

class SecondaryVertexMaker : public edm::EDProducer {
public:
  explicit SecondaryVertexMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  GlobalVector direction(const reco::Vertex& final, const reco::Vertex& original);
  double angle(const reco::Vertex& final, const reco::Vertex& original, LorentzVector p4);

  // ----------member data ---------------------------
  edm::InputTag primaryVertexInputTag_;
  edm::InputTag inclusiveVertexInputTag_;
  std::string aliasprefix_;
  const reco::Vertex* primaryVertex_;
};

//
// constructors and destructor
//
SecondaryVertexMaker::SecondaryVertexMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<std::vector<LorentzVector> >       (branchprefix+"position"  ).setBranchAlias(aliasprefix_+"_position"   );
  produces<std::vector<float> >               (branchprefix+"xError"    ).setBranchAlias(aliasprefix_+"_xError"     );
  produces<std::vector<float> >               (branchprefix+"yError"    ).setBranchAlias(aliasprefix_+"_yError"     );
  produces<std::vector<float> >               (branchprefix+"zError"    ).setBranchAlias(aliasprefix_+"_zError"     );
  produces<std::vector<float> >               (branchprefix+"chi2"      ).setBranchAlias(aliasprefix_+"_chi2"       ); 
  produces<std::vector<float> >               (branchprefix+"ndof"      ).setBranchAlias(aliasprefix_+"_ndof"       );
  produces<std::vector<float> >               (branchprefix+"prob"      ).setBranchAlias(aliasprefix_+"_prob"       );   
  produces<std::vector<int>   >               (branchprefix+"isKs"      ).setBranchAlias(aliasprefix_+"_isKs"       );
  produces<std::vector<int>   >               (branchprefix+"isLambda"  ).setBranchAlias(aliasprefix_+"_isLambda"   );
  produces<std::vector<int>   >               (branchprefix+"nTrks"     ).setBranchAlias(aliasprefix_+"_nTrks"      ); // number of tracks used in the fit
  produces<std::vector<LorentzVector> >       (branchprefix+"p4"        ).setBranchAlias(aliasprefix_+"_p4"         ); // sum of the vertex track p4
  produces<std::vector<LorentzVector> >       (branchprefix+"flight"    ).setBranchAlias(aliasprefix_+"_flight"     ); // flight direction wrt primary vertex
  produces<std::vector<LorentzVector> >       (branchprefix+"refitp4"   ).setBranchAlias(aliasprefix_+"_refitp4"    ); // vertex constraint p4
  produces<std::vector<float> >               (branchprefix+"dist3Dval" ).setBranchAlias(aliasprefix_+"_dist3Dval"  ); // distance from primary vertex (3D, value)
  produces<std::vector<float> >               (branchprefix+"dist3Dsig" ).setBranchAlias(aliasprefix_+"_dist3Dsig"  ); // distance from primary vertex (3D, significance)
  produces<std::vector<float> >               (branchprefix+"distXYval" ).setBranchAlias(aliasprefix_+"_distXYval"  ); // distance from primary vertex (XY, value)
  produces<std::vector<float> >               (branchprefix+"distXYsig" ).setBranchAlias(aliasprefix_+"_distXYsig"  ); // distance from primary vertex (XY, significance)
  produces<std::vector<float> >               (branchprefix+"anglePV"   ).setBranchAlias(aliasprefix_+"_anglePV"    ); // angle between vertex momentum and direction from the primary vertex   
  produces<std::vector<int> >                 (branchprefix+"mc3id"     ).setBranchAlias(aliasprefix_+"_mc3_id"     ); // direction matching for hard scatter gen particles
  produces<std::vector<LorentzVector> >       (branchprefix+"mc3p4"     ).setBranchAlias(aliasprefix_+"_mc3_p4"     ); //   --- // ---

  primaryVertexInputTag_   = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
  inclusiveVertexInputTag_ = iConfig.getParameter<edm::InputTag>("inclusiveVertexInputTag");
}

void SecondaryVertexMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get the primary vertices
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(primaryVertexInputTag_, primaryVertices);
  if (! primaryVertices.isValid() ) {
    edm::LogError("SecondaryVertexMakerError") << "Error! can't get the primary vertices";
  }
  edm::Handle<reco::VertexCollection> secVertices;
  iEvent.getByLabel(inclusiveVertexInputTag_, secVertices);
  if (! secVertices.isValid() ) {
    edm::LogError("SecondaryVertexMakerError") << "Error! can't get the secondary vertices";
  }
  edm::Handle<reco::GenParticleCollection> genps;
  iEvent.getByLabel("genParticles",genps);

  primaryVertex_ = 0;
  double primaryVertexSumPt(0);
  for ( reco::VertexCollection::const_iterator vtx = primaryVertices->begin();
        vtx != primaryVertices->end(); ++vtx){
    if ( !vtx->isValid() || vtx->isFake() ) continue;
    double sumpt(0);
    for (reco::Vertex::trackRef_iterator trk = vtx->tracks_begin();
         trk != vtx->tracks_end(); ++trk)
      sumpt += (*trk)->pt();
    if ( sumpt > primaryVertexSumPt ){
      primaryVertexSumPt = sumpt;
      primaryVertex_ = &*vtx;
    }
  }

  std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_position          (new std::vector<LorentzVector>       );
  std::auto_ptr<std::vector<float> >               vector_vtxs_xError            (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_yError            (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_zError            (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_chi2              (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_ndof              (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_prob              (new std::vector<float>               );
  std::auto_ptr<std::vector<int>   >               vector_vtxs_isKs              (new std::vector<int>                 );
  std::auto_ptr<std::vector<int>   >               vector_vtxs_isLambda          (new std::vector<int>                 );
  std::auto_ptr<std::vector<int>   >               vector_vtxs_nTrk              (new std::vector<int>                 );
  std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_p4                (new std::vector<LorentzVector>       );
  std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_flight            (new std::vector<LorentzVector>       );
  std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_refitp4           (new std::vector<LorentzVector>       );
  std::auto_ptr<std::vector<float> >               vector_vtxs_dist3Dval         (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_dist3Dsig         (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_distXYval         (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_distXYsig         (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_anglePV           (new std::vector<float>               );
  std::auto_ptr<std::vector<int> >                 vector_vtxs_mc3id             (new std::vector<int>                 );
  std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_mc3p4             (new std::vector<LorentzVector>       );
     
  if ( primaryVertex_ ){
    for ( std::vector<reco::Vertex>::const_iterator vertex=secVertices->begin();
	  vertex!=secVertices->end(); ++vertex )
      {
	vector_vtxs_position->push_back( LorentzVector( vertex->position().x(), vertex->position().y(), vertex->position().z(), 0 ) );
	vector_vtxs_prob->push_back( TMath::Prob(vertex->chi2(), (int)vertex->ndof()) );

	std::vector<const reco::Track*> vtxtracks;
	LorentzVector p4;
	for(std::vector<reco::TrackBaseRef>::const_iterator trk=vertex->tracks_begin();
	    trk!=vertex->tracks_end(); ++trk){
	  vtxtracks.push_back(trk->get());
	  p4 += math::PtEtaPhiMLorentzVector((*trk)->pt(), (*trk)->eta(), (*trk)->phi(),0.139570);;
	}
	vector_vtxs_p4->push_back( p4 );

	LorentzVector p4Ref;
	if (vertex->hasRefittedTracks()){
	  const std::vector<reco::Track>& reftrks(vertex->refittedTracks());
	  for(std::vector<reco::Track>::const_iterator trk=reftrks.begin();
	      trk!=reftrks.end();++trk){
	    p4Ref += math::PtEtaPhiMLorentzVector(trk->pt(), trk->eta(), trk->phi(), 0.139570);
	  }
	}
	vector_vtxs_refitp4->push_back( p4Ref );

	// distance and errors
	VertexDistance3D vdist3D;
	Measurement1D dist3D = vdist3D.distance(*vertex,*primaryVertex_);
	vector_vtxs_dist3Dval->push_back(dist3D.value());
	vector_vtxs_dist3Dsig->push_back(dist3D.significance());
	VertexDistanceXY vdistXY;
	Measurement1D distXY = vdistXY.distance(*vertex,*primaryVertex_);
	vector_vtxs_distXYval->push_back(distXY.value());
	vector_vtxs_distXYsig->push_back(distXY.significance());

	vector_vtxs_anglePV->push_back(angle(*vertex, *primaryVertex_, p4));
	vector_vtxs_nTrk->push_back(vtxtracks.size());

	vector_vtxs_xError->push_back( vertex->xError() );
	vector_vtxs_yError->push_back( vertex->yError() );
	vector_vtxs_zError->push_back( vertex->zError() );
	vector_vtxs_chi2->push_back( vertex->chi2() );
	vector_vtxs_ndof->push_back( vertex->ndof() );
	
	bool kshort(false);
	bool lambda(false);
	if ( vtxtracks.size() == 2 ){
	  if (vertex->hasRefittedTracks()){
	    // good resolution (~7MeV for Ks mass)
	    const std::vector<reco::Track>& reftrks(vertex->refittedTracks());
	    math::PtEtaPhiMLorentzVector pion1(  reftrks.at(0).pt(), reftrks.at(0).eta(), reftrks.at(0).phi(), 0.139570);
	    math::PtEtaPhiMLorentzVector pion2(  reftrks.at(1).pt(), reftrks.at(1).eta(), reftrks.at(1).phi(), 0.139570);
	    math::PtEtaPhiMLorentzVector proton1(reftrks.at(0).pt(), reftrks.at(0).eta(), reftrks.at(0).phi(), 0.938272);
	    math::PtEtaPhiMLorentzVector proton2(reftrks.at(1).pt(), reftrks.at(1).eta(), reftrks.at(1).phi(), 0.938272);
	    if ( fabs((pion1+pion2).mass()-0.498) < 0.03 )   kshort = true;
	    if ( fabs((proton1+pion2).mass()-1.116) < 0.01 ) lambda = true;
	    if ( fabs((proton2+pion1).mass()-1.116) < 0.01 ) lambda = true;
	  } else {
	    // resolution is much worse (~50MeV for Ks mass)
	    math::PtEtaPhiMLorentzVector pion1(  vtxtracks.at(0)->pt(), vtxtracks.at(0)->eta(), vtxtracks.at(0)->phi(), 0.139570);
	    math::PtEtaPhiMLorentzVector pion2(  vtxtracks.at(1)->pt(), vtxtracks.at(1)->eta(), vtxtracks.at(1)->phi(), 0.139570);
	    math::PtEtaPhiMLorentzVector proton1(vtxtracks.at(0)->pt(), vtxtracks.at(0)->eta(), vtxtracks.at(0)->phi(), 0.938272);
	    math::PtEtaPhiMLorentzVector proton2(vtxtracks.at(1)->pt(), vtxtracks.at(1)->eta(), vtxtracks.at(1)->phi(), 0.938272);
	    if ( fabs((pion1+pion2).mass()-0.498) < 0.10 )   kshort = true;
	    if ( fabs((proton1+pion2).mass()-1.116) < 0.03 ) lambda = true;
	    if ( fabs((proton2+pion1).mass()-1.116) < 0.03 ) lambda = true;
	  }
	}
	vector_vtxs_isKs->push_back(kshort);
	vector_vtxs_isLambda->push_back(lambda);

	int genid(0);
	LorentzVector genp4;
	LorentzVector flightDir;
	if ( genps.isValid() ){
	  double drMin = 999;
	  for (reco::GenParticleCollection::const_iterator ps = genps->begin();
	       ps != genps->end(); ++ps)
	    {
	      if ( ps->status()!=3 ) continue;
	      GlobalVector dir = direction(*vertex,*primaryVertex_);
	      flightDir = LorentzVector(dir.x(),dir.y(),dir.z(),0);
	      double dr = reco::deltaR(dir,*ps);
	      if ( dr > 0.5 ) continue;
	      if ( dr > drMin ) continue;
	      genid = ps->pdgId();
	      genp4 = ps->p4();
	      drMin = dr;
	    }
	}
	vector_vtxs_mc3id->push_back(genid);
	vector_vtxs_mc3p4->push_back(genp4);
	vector_vtxs_flight->push_back(flightDir);
      }
  }
  // store into the event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(vector_vtxs_position,          branchprefix+"position"          );
  iEvent.put(vector_vtxs_xError,            branchprefix+"xError"            );
  iEvent.put(vector_vtxs_yError,            branchprefix+"yError"            );
  iEvent.put(vector_vtxs_zError,            branchprefix+"zError"            );
  iEvent.put(vector_vtxs_chi2,              branchprefix+"chi2"              );
  iEvent.put(vector_vtxs_ndof,              branchprefix+"ndof"              );
  iEvent.put(vector_vtxs_prob,              branchprefix+"prob"              );
  iEvent.put(vector_vtxs_isKs,              branchprefix+"isKs"              );
  iEvent.put(vector_vtxs_isLambda,          branchprefix+"isLambda"          );
  iEvent.put(vector_vtxs_nTrk,              branchprefix+"nTrks"             );
  iEvent.put(vector_vtxs_p4,                branchprefix+"p4"                );
  iEvent.put(vector_vtxs_flight,            branchprefix+"flight"            );
  iEvent.put(vector_vtxs_refitp4,           branchprefix+"refitp4"           );
  iEvent.put(vector_vtxs_dist3Dval,         branchprefix+"dist3Dval"         );
  iEvent.put(vector_vtxs_dist3Dsig,         branchprefix+"dist3Dsig"         );
  iEvent.put(vector_vtxs_distXYval,         branchprefix+"distXYval"         );
  iEvent.put(vector_vtxs_distXYsig,         branchprefix+"distXYsig"         );
  iEvent.put(vector_vtxs_anglePV,           branchprefix+"anglePV"           );
  iEvent.put(vector_vtxs_mc3id,             branchprefix+"mc3id"             );
  iEvent.put(vector_vtxs_mc3p4,             branchprefix+"mc3p4"     );
}

// ------------ method called once each job just before starting event loop  ------------
void SecondaryVertexMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void SecondaryVertexMaker::endJob() 
{
}

GlobalVector
SecondaryVertexMaker::direction(const reco::Vertex& final, const reco::Vertex& original)
{
  return GlobalPoint(Basic3DVector<float>(final.position())) - 
    GlobalPoint(Basic3DVector<float>(original.position()));
}
double
SecondaryVertexMaker::angle(const reco::Vertex& final, const reco::Vertex& original, 
				  LorentzVector p4){
  GlobalVector dir = direction(final,original);
  return acos((dir.x()*p4.x()+dir.y()*p4.y()+dir.z()*p4.z())/dir.mag()/p4.P());
}

//define this as a plug-in
DEFINE_FWK_MODULE(SecondaryVertexMaker);

