#include <vector>
#include <functional>
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "CMS2/NtupleMaker/interface/HypIsoMaker.h"
// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace edm;

HypIsoMaker::HypIsoMaker(const edm::ParameterSet& iConfig) 
{
  //inputs: edm collections
  emObjectProducer_               = iConfig.getParameter<edm::InputTag>("emObjectProducer"); //'uniqueElectrons'
  muonsInputTag_                  = iConfig.getParameter<edm::InputTag>("muonsInputTag");
  //cms2 branches
  cms2elsInputTag_                = iConfig.getParameter<edm::InputTag>("cms2elsInputTag");
  cms2musInputTag_                = iConfig.getParameter<edm::InputTag>("cms2musInputTag");
  trackInputTag_                  = iConfig.getParameter<edm::InputTag>("trackInputTag");
  hypInputTag_                    = iConfig.getParameter<edm::InputTag>("hypInputTag");
  ecalBarrelRecHitProducer_       = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHitProducer");
  ecalBarrelRecHitCollection_     = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHitCollection");
  ecalEndcapRecHitProducer_       = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHitProducer");
  ecalEndcapRecHitCollection_     = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHitCollection");
  caloTowersInputTag_             = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");
   
  //recompute variables--whether or not to use els/mus_iso vars
  recomputeEcalIso_ = iConfig.getParameter<bool>("recomputeEcalIso");
  recomputeTrckIso_ = iConfig.getParameter<bool>("recomputeTrckIso");
	
  //vetos
  elsEcalVetoRadBarrel_  = iConfig.getParameter<double>("elsEcalVetoRadBarrel"  );
  elsEcalVetoRadEndcap_  = iConfig.getParameter<double>("elsEcalVetoRadEndcap"	);
  elsEcalExtCone_        = iConfig.getParameter<double>("elsEcalExtCone"		);
  musEcalVetoRadBarrel_  = iConfig.getParameter<double>("musEcalVetoRadBarrel"	);
  musEcalVetoRadEndcap_  = iConfig.getParameter<double>("musEcalVetoRadEndcap"	);
  musEcalExtCone_        = iConfig.getParameter<double>("musEcalExtCone"        );
  
  IsoJurassicWidth_   	= iConfig.getParameter<double>("jurassicWidth");
  
  elsEtMinBarrel_       = iConfig.getParameter<double>("elsetMinBarrel");
  elsEMinBarrel_        = iConfig.getParameter<double>("elseMinBarrel");
  elsEtMinEndcap_       = iConfig.getParameter<double>("elsetMinEndcap");
  elsEMinEndcap_        = iConfig.getParameter<double>("elseMinEndcap");

  musEtMinBarrel_       = iConfig.getParameter<double>("musetMinBarrel");
  musEMinBarrel_        = iConfig.getParameter<double>("museMinBarrel");
  musEtMinEndcap_       = iConfig.getParameter<double>("musetMinEndcap");
  musEMinEndcap_        = iConfig.getParameter<double>("museMinEndcap");

  trackIsoExtRadius_    = iConfig.getParameter<double>("trackIsoExtRadius");
  trackIsoElsInRadius_  = iConfig.getParameter<double>("trackIsoElsInRadius");
  trackIsoMusInRadius_  = iConfig.getParameter<double>("trackIsoMusInRadius");
  trackIsoMinPt_        = iConfig.getParameter<double>("trackIsoMinPt");
  trackIsoMind0_        = iConfig.getParameter<double>("trackIsoMind0");
  trackIsoMinz0_        = iConfig.getParameter<double>("trackIsoMinz0");

  // options
  useIsolEt_ = iConfig.getParameter<bool>("useIsolEt");

  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  muonparameters_.loadParameters( parameters );

  //register your products
  produces<vector<float> >         ("hypltecaliso"             ).setBranchAlias("hyp_lt_ecaliso");
  produces<vector<float> >         ("hypllecaliso"             ).setBranchAlias("hyp_ll_ecaliso");
  produces<vector<float> >         ("hyplttrkiso"              ).setBranchAlias("hyp_lt_trkiso");
  produces<vector<float> >         ("hyplltrkiso"              ).setBranchAlias("hyp_ll_trkiso");
}


HypIsoMaker::~HypIsoMaker(){}


// ------------ method called to produce the data  ------------
void
HypIsoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get the filtered electrons
  edm::Handle< edm::View<reco::Candidate> > emObjectHandle;
  iEvent.getByLabel(emObjectProducer_,emObjectHandle);

  edm::Handle< edm::View<reco::Muon> > muonHandle;
  iEvent.getByLabel(muonsInputTag_, muonHandle);

  //get CMS2 branches

  // hyp_lt p4
  InputTag hyp_lt_p4_tag(hypInputTag_.label(),"hypltp4");
  Handle<vector<LorentzVector> > hyp_lt_p4_h;
  iEvent.getByLabel(hyp_lt_p4_tag, hyp_lt_p4_h);
  const vector<LorentzVector> *hyp_lt_p4 = hyp_lt_p4_h.product();

  // hyp_ll p4
  InputTag hyp_ll_p4_tag(hypInputTag_.label(),"hypllp4");
  Handle<vector<LorentzVector> > hyp_ll_p4_h;
  iEvent.getByLabel(hyp_ll_p4_tag, hyp_ll_p4_h);
  const vector<LorentzVector> *hyp_ll_p4 = hyp_ll_p4_h.product();

  // hyp_lt d0
  InputTag hyp_lt_d0_tag(hypInputTag_.label(),"hypltd0");
  Handle<vector<float> > hyp_lt_d0_h;
  iEvent.getByLabel(hyp_lt_d0_tag, hyp_lt_d0_h);
  const vector<float> *hyp_lt_d0 = hyp_lt_d0_h.product();

  // hyp_ll d0
  InputTag hyp_ll_d0_tag(hypInputTag_.label(),"hyplld0");
  Handle<vector<float> > hyp_ll_d0_h;
  iEvent.getByLabel(hyp_ll_d0_tag, hyp_ll_d0_h);
  const vector<float> *hyp_ll_d0 = hyp_ll_d0_h.product();

  // hyp_lt z0
  InputTag hyp_lt_z0_tag(hypInputTag_.label(),"hypltz0");
  Handle<vector<float> > hyp_lt_z0_h;
  iEvent.getByLabel(hyp_lt_z0_tag, hyp_lt_z0_h);
  const vector<float> *hyp_lt_z0 = hyp_lt_z0_h.product();

  // hyp_ll z0
  InputTag hyp_ll_z0_tag(hypInputTag_.label(),"hypllz0");
  Handle<vector<float> > hyp_ll_z0_h;
  iEvent.getByLabel(hyp_ll_z0_tag, hyp_ll_z0_h);
  const vector<float> *hyp_ll_z0 = hyp_ll_z0_h.product();

  // hyp_lt_id -- +/-13 or +/-11 (no mc match)
  edm::InputTag hyp_lt_id_tag(hypInputTag_.label(),"hypltid");
  edm::Handle<std::vector<int> > hyp_lt_id_h;
  iEvent.getByLabel(hyp_lt_id_tag, hyp_lt_id_h);
  const vector<int> *hyp_lt_id = hyp_lt_id_h.product();

  // hyp_ll_id -- +/-13 or +/-11 (no mc match)
  edm::InputTag hyp_ll_id_tag(hypInputTag_.label(),"hypllid");
  edm::Handle<std::vector<int> > hyp_ll_id_h;
  iEvent.getByLabel(hyp_ll_id_tag, hyp_ll_id_h);
  const vector<int> *hyp_ll_id = hyp_ll_id_h.product();

  //hyp_lt_index
  edm::InputTag hyp_lt_index_tag(hypInputTag_.label(),"hypltindex");
  edm::Handle<std::vector<int> > hyp_lt_index_h;
  iEvent.getByLabel(hyp_lt_index_tag, hyp_lt_index_h);
  const vector<int> *hyp_lt_index = hyp_lt_index_h.product();

  //hyp_ll_index
  edm::InputTag hyp_ll_index_tag(hypInputTag_.label(),"hypllindex");
  edm::Handle<std::vector<int> > hyp_ll_index_h;
  iEvent.getByLabel(hyp_ll_index_tag, hyp_ll_index_h);
  const vector<int> *hyp_ll_index = hyp_ll_index_h.product();

  // electrons
  edm::InputTag els_ecalIso_tag(cms2elsInputTag_.label(), "elsecalIso");
  Handle<vector<float> > els_ecalIso_h;
  iEvent.getByLabel(els_ecalIso_tag, els_ecalIso_h);
  const vector<float> *els_ecalIso = els_ecalIso_h.product();

  edm::InputTag els_tkIso_tag(cms2elsInputTag_.label(), "elstkIso");
  Handle<vector<float> > els_tkIso_h;
  iEvent.getByLabel(els_tkIso_tag, els_tkIso_h);
  const vector<float> *els_tkIso = els_tkIso_h.product();

  // muons
  edm::InputTag mus_iso03_sumPt_tag(cms2musInputTag_.label(), "musiso03sumPt");
  Handle<vector<float> > mus_iso03_sumPt_h;
  iEvent.getByLabel(mus_iso03_sumPt_tag, mus_iso03_sumPt_h);
  const vector<float> *mus_iso03_sumPt = mus_iso03_sumPt_h.product();

  edm::InputTag mus_iso03_emEt_tag(cms2musInputTag_.label(), "musiso03emEt");
  Handle<vector<float> > mus_iso03_emEt_h;
  iEvent.getByLabel(mus_iso03_emEt_tag, mus_iso03_emEt_h);
  const vector<float> *mus_iso03_emEt = mus_iso03_emEt_h.product();

  //cms2 track branches
  //track p4
  InputTag trks_trk_p4_tag(trackInputTag_.label(),"trkstrkp4");
  Handle<vector<LorentzVector> > trks_trk_p4_h;
  iEvent.getByLabel(trks_trk_p4_tag, trks_trk_p4_h);
  const vector<LorentzVector> *trks_trk_p4 = trks_trk_p4_h.product();

  //track d0
  InputTag trks_d0_tag(trackInputTag_.label(),"trksd0");
  Handle<vector<float> > trks_d0_h;
  iEvent.getByLabel(trks_d0_tag, trks_d0_h);
  const vector<float> *trks_d0 = trks_d0_h.product();

  //track z0
  InputTag trks_z0_tag(trackInputTag_.label(),"trksz0");
  Handle<vector<float> > trks_z0_h;
  iEvent.getByLabel(trks_z0_tag, trks_z0_h);
  const vector<float> *trks_z0 = trks_z0_h.product();

  auto_ptr<vector<float> > els_newiso                   (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_ecaliso               (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_ecaliso               (new vector<float>);
  auto_ptr<vector<float> > hyp_lt_trkiso                (new vector<float>);
  auto_ptr<vector<float> > hyp_ll_trkiso                (new vector<float>);

  //track isolation

  double trackltiso = 0;

  LorentzVector prip4;
  vector<LorentzVector> excp4;
  vector<int> excid;
  float pri_d0;
  float pri_z0;
  int priid;
  //debugging vars
  //bool printevt = false;
  //int prltidx = 0;
  //int prllidx = 0;
  
  //Dilepton hypotheses : lt
  for( unsigned int i=0; i < hyp_lt_p4->size(); i++ ) {
	prip4 = hyp_lt_p4->at(i);
	pri_d0 = hyp_lt_d0->at(i);
	pri_z0 = hyp_lt_z0->at(i);
	priid = hyp_lt_id->at(i);
	excp4.push_back(hyp_ll_p4->at(i));
	excid.push_back(hyp_ll_id->at(i));
	trackltiso = track_iso( prip4, pri_d0, pri_z0, priid, excp4, excid, trks_trk_p4, trks_d0, trks_z0);

	if( recomputeTrckIso_ )
	  hyp_lt_trkiso->push_back( trackltiso );
	else { // use existing variable
	  if( TMath::Abs(priid) == 11 )
		hyp_lt_trkiso->push_back( els_tkIso->at( hyp_lt_index->at(i) ) - trackltiso );
	  else if( TMath::Abs(priid) == 13 )
		hyp_lt_trkiso->push_back( mus_iso03_sumPt->at( hyp_lt_index->at(i) ) - trackltiso );
	}
  
	excp4.clear(); //don't accumulate for next hypothesis
	excid.clear();
	
	//if( abs( trackltiso - hyp_lt_tkIso->at(i) ) > 0.1 ) { //debugging
	  //printevt = true;
	  //cout << "my lt iso  egm lt iso\n";
	  //cout << trackltiso << "  " << hyp_lt_tkIso->at(i) << endl;
	  //prltidx = i;
	  //}
  }

  double tracklliso = 0;

  //Dilepton hypotheses : ll
  for( unsigned int i=0; i < hyp_ll_p4->size(); i++ ) {
	prip4 = hyp_ll_p4->at(i);
	pri_d0 = hyp_ll_d0->at(i);
	pri_z0 = hyp_ll_z0->at(i);
	priid = hyp_ll_id->at(i);
	excp4.push_back(hyp_lt_p4->at(i));
	excid.push_back(hyp_lt_id->at(i));
	tracklliso = track_iso( prip4, pri_d0, pri_z0, priid, excp4, excid, trks_trk_p4, trks_d0, trks_z0);

	if( recomputeTrckIso_ )
	  hyp_ll_trkiso->push_back( tracklliso );
	else { //use existing variable
	  if( TMath::Abs(priid) == 11 )
		hyp_ll_trkiso->push_back( els_tkIso->at( hyp_ll_index->at(i) ) - tracklliso );
	  else if( TMath::Abs(priid) == 13 )
		hyp_ll_trkiso->push_back( mus_iso03_sumPt->at( hyp_ll_index->at(i) ) - tracklliso );
	}

	excp4.clear();
	excid.clear();

	//if( abs( tracklliso - hyp_ll_tkIso->at(i) ) > 0.1 ) {
	  //printevt = true;
	  //cout << "my ll iso  egm ll iso\n";
	  //cout << tracklliso << "  " << hyp_ll_tkIso->at(i) << endl;
	  //prllidx = i;
	//}
  }

  //debugging
  //if( printevt ) {
  //	cout << "hyp_lt pt eta phi: " << hyp_lt_p4->at(prltidx).pt() << "  " << hyp_lt_p4->at(prltidx).eta() << "  " << hyp_lt_p4->at(prltidx).phi() << endl;
  //	cout << "hyp_ll pt eta phi: " << hyp_ll_p4->at(prllidx).pt() << "  " << hyp_ll_p4->at(prllidx).eta() << "  " << hyp_ll_p4->at(prllidx).phi() << endl;
  //	cout << "\nidx  trk_pt  d lt d0   d lt z0    r2lt      d ll d0    d ll z0    r2ll \n";
  //	for(unsigned int i=0;i<trks_trk_p4->size();i++) {
  //	  if( trks_trk_p4->at(i).pt() >= 1.0 ) {
  //		double r2lt = ROOT::Math::VectorUtil::DeltaR(hyp_lt_p4->at(prltidx), trks_trk_p4->at(i));
  //		double d0lt = TMath::Abs( hyp_lt_d0->at(prltidx) - trks_d0->at(i) );
  //		double z0lt = TMath::Abs( hyp_lt_z0->at(prltidx) - trks_z0->at(i) );
  //
  //		double r2ll = ROOT::Math::VectorUtil::DeltaR(hyp_ll_p4->at(prllidx), trks_trk_p4->at(i));
  //		double d0ll = TMath::Abs( hyp_ll_d0->at(prllidx) - trks_d0->at(i) );
  //		double z0ll = TMath::Abs( hyp_ll_z0->at(prllidx) - trks_z0->at(i) );
  //		cout << i << "  " << trks_trk_p4->at(i).pt() << "  " << d0lt << "  " << z0lt << "  " << r2lt << "  " << d0ll << "  " << z0ll << "  " << r2ll << endl;
  //	  }
  //	}
  //	cout << endl;
  //}


  //new ecal iso from hyp
  
  //cout << "\nmake sure els block is aligned with uniqueElectrons collection\n";
  //for( unsigned int i=0;i<els_p4->size();i++) {
  //	cout << "gsfel: " << emObjectHandle->at(i).pt() << "  " << emObjectHandle->at(i).eta() << "  " << emObjectHandle->at(i).phi() << endl;
  //	cout << "ntpel: " << els_p4->at(i).pt() << "  " << els_p4->at(i).eta() << "  " << els_p4->at(i).phi() << endl;
  //}
  //if( hyp_lt_p4->size() > 7 ) { //another debug if
  //cout << "\nhypidx  ltpt   llpt   ltid   llid    dr\n";
  //}
  
  vector<int> excidx;
  //Dilepton hypotheses : lt
  for( unsigned int i=0; i < hyp_lt_p4->size(); i++ ) {
	excid.push_back(hyp_ll_id->at(i));
	excidx.push_back(hyp_ll_index->at(i));
	double isoValue = 0;

	//double dr = ROOT::Math::VectorUtil::DeltaR( hyp_lt_p4->at(i), hyp_ll_p4->at(i) ); //debugging
	//if( hyp_lt_p4->size() > 7 || dr < 0.5 ) { //my event is the only one with so many hyps of the 100 i loop on
	//  cout << i << "  " << hyp_lt_p4->at(i).pt() << "  " << hyp_ll_p4->at(i).pt() << "  " << hyp_lt_id->at(i) << "  " << hyp_ll_id->at(i) << "  " << dr << endl;
	//}

	if(useIsolEt_) isoValue = getHypSum( hyp_lt_id->at(i), hyp_lt_index->at(i), excid, excidx, true, emObjectHandle, muonHandle, iEvent, iSetup );
	
 	else           isoValue = getHypSum( hyp_lt_id->at(i), hyp_lt_index->at(i), excid, excidx, false, emObjectHandle, muonHandle, iEvent, iSetup );
	

	excid.clear(); //must clear before relooping
	excidx.clear();

	if( recomputeTrckIso_ )
	  hyp_lt_ecaliso->push_back( isoValue );
	else { //use existing variable
	  if( TMath::Abs(hyp_lt_id->at(i)) == 11 )
		hyp_lt_ecaliso->push_back( els_ecalIso->at( hyp_lt_index->at(i) ) - isoValue );
	  else if( TMath::Abs(hyp_lt_id->at(i)) == 13 )
		hyp_lt_ecaliso->push_back( mus_iso03_emEt->at( hyp_lt_index->at(i) ) - isoValue );
	}
  }


  //Dilepton hypotheses : ll
  for( unsigned int i=0; i < hyp_ll_p4->size(); i++ ) {
	excid.push_back(hyp_lt_id->at(i));
	excidx.push_back(hyp_lt_index->at(i));
	double isoValue = 0;

	//double dr = ROOT::Math::VectorUtil::DeltaR( hyp_lt_p4->at(i), hyp_ll_p4->at(i) ); //debugging
	//if( hyp_lt_p4->size() > 7 || dr < 0.5 ) { //my event is the only one with so many hyps of the 100 i loop on
	//  cout << i << "  " << hyp_lt_p4->at(i).pt() << "  " << hyp_ll_p4->at(i).pt() << "  " << hyp_lt_id->at(i) << "  " << hyp_ll_id->at(i) << "  " << dr << endl;
	//}
	
	if(useIsolEt_) isoValue = getHypSum( hyp_ll_id->at(i), hyp_ll_index->at(i), excid, excidx, true, emObjectHandle, muonHandle, iEvent, iSetup );
	
 	else           isoValue = getHypSum( hyp_ll_id->at(i), hyp_ll_index->at(i), excid, excidx, false, emObjectHandle, muonHandle, iEvent, iSetup );
	
	excid.clear(); //must clear before relooping
	excidx.clear();

	if( recomputeTrckIso_ )
	  hyp_ll_ecaliso->push_back( isoValue );
	else { //use existing variable
	  if( TMath::Abs(hyp_ll_id->at(i)) == 11 )
		hyp_ll_ecaliso->push_back( els_ecalIso->at( hyp_ll_index->at(i) ) - isoValue );
	  else if( TMath::Abs(hyp_lt_id->at(i)) == 13 )
		hyp_ll_ecaliso->push_back( mus_iso03_emEt->at( hyp_ll_index->at(i) ) - isoValue );
	}
  }


  iEvent.put(hyp_lt_ecaliso, "hypltecaliso");
  iEvent.put(hyp_ll_ecaliso, "hypllecaliso");
  iEvent.put(hyp_lt_trkiso, "hyplttrkiso");
  iEvent.put(hyp_ll_trkiso, "hyplltrkiso");

}

//Ecal isolation for electrons and muons which excludes second hyp
double HypIsoMaker::getHypSum( int objid, int objidx, vector<int> excid, vector<int> excidx, bool returnEt, edm::Handle< edm::View<reco::Candidate> > emObjectHandle, edm::Handle< edm::View<reco::Muon> > muonHandle, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{

  double energySum = 0.;
  double vetoSum = 0.;
  double extRadius = 0.;
  double intRadius = 0.;
  GlobalPoint pclu(0.,0.,0.);
  double etaclus = 0;
  double phiclus = 0;
  //vectors of hits (in case of electrons) or towers (muons) to loop over
  vector<GlobalPoint> energypos;
  vector<double> energyval;

  //set up above vars based on el/mu
  if( abs( objid ) == 11 ) { //electron
	extRadius = elsEcalExtCone_;
	//Take the SC position
	const reco::Candidate* emObject = &(emObjectHandle->at(objidx)); //get gsf electron from uniqueElectrons which corresponds to hyp
	reco::SuperClusterRef sc = emObject->get<reco::SuperClusterRef>();
	math::XYZPoint theCaloPosition = sc.get()->position();
	pclu = GlobalPoint( theCaloPosition.x(), theCaloPosition.y(), theCaloPosition.z() );
	etaclus = pclu.eta(); 
	phiclus = pclu.phi();

	if( TMath::Abs( etaclus ) <= 1.479 )
	  intRadius = elsEcalVetoRadBarrel_;
	else
	  intRadius = elsEcalVetoRadEndcap_;

	// Get Ecal hits barrel
	edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle;
	iEvent.getByLabel(ecalBarrelRecHitProducer_.label(),ecalBarrelRecHitCollection_.label(), ecalBarrelRecHitHandle);

	// Get Ecal hits endcap
	edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
	iEvent.getByLabel(ecalEndcapRecHitProducer_.label(), ecalEndcapRecHitCollection_.label(),ecalEndcapRecHitHandle);
  
	//create the meta hit collections inorder that we can pass them into the isolation objects
	EcalRecHitMetaCollection ecalBarrelHits(*ecalBarrelRecHitHandle);
	EcalRecHitMetaCollection ecalEndcapHits(*ecalEndcapRecHitHandle);

	if( !&ecalBarrelHits || !&ecalEndcapHits ) {
	  cout << "HypIsoMaker : ERROR bad ecal hits collections\n\n";
	  return 0; //check hits collection
	}

	//Get Calo Geometry
	edm::ESHandle<CaloGeometry> pG;
	iSetup.get<CaloGeometryRecord>().get(pG);
	const CaloGeometry* caloGeom = pG.product();
	
	//set up the geometry and selector
	const CaloSubdetectorGeometry* subdet[2]; // barrel+endcap
	subdet[0] = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
	subdet[1] = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);


	//set position vector--use hits for electrons
	for(int subdetnr=0; subdetnr<=1 ; subdetnr++){  // look in barrel and endcap
	  double etLow = 0;
	  double eLow = 0;
	  if( subdetnr == 0 ) { //0 is Barrel
		etLow = elsEtMinBarrel_;
		eLow = elsEMinBarrel_;
		//if( abs( objid ) == 11 )      intRadius = elsEcalVetoRadBarrel_;
		//else if( abs( objid ) == 13 ) intRadius = musEcalVetoRadBarrel_;
	  }
	  else if( subdetnr == 1 ) { //1 is Endcap
		etLow = elsEtMinEndcap_;
		eLow = elsEMinEndcap_;
		//if( abs( objid ) == 11 )      intRadius = elsEcalVetoRadEndcap_;
		//else if( abs( objid ) == 13 ) intRadius = musEcalVetoRadEndcap_;
	  }
	  else {
		cout << "HypIsoMaker : ERROR bad subdet(ector) loop index--should be 0 or 1 \n\n";
		return 0;
	  }

	  CaloSubdetectorGeometry::DetIdSet chosen = subdet[subdetnr]->getCells(pclu, extRadius);// select cells around cluster
	  CaloRecHitMetaCollectionV::const_iterator bhits = ecalBarrelHits.end();
	  CaloRecHitMetaCollectionV::const_iterator ehits = ecalEndcapHits.end();
	  for (CaloSubdetectorGeometry::DetIdSet::const_iterator i = chosen.begin ();i!= chosen.end ();++i){//loop selected cells

		bhits = ecalBarrelHits.find(*i); // find selected cell among rechits
		ehits = ecalEndcapHits.find(*i); // find selected cell among rechits
		
		if( bhits != ecalBarrelHits.end() ) { //check if the detid (cell) in question has a hit
			  
		  const GlobalPoint & position = caloGeom->getPosition(*i);
		  double energy = bhits->energy();
		  double et = energy*position.perp()/position.mag();
		  //if( et == 0. && energy == 0. ) cout << "zero energy in cell " << endl; //debug
		  if( et < etLow || energy < eLow ) continue; //thresholds on e/et of cells

		  energypos.push_back( position );
		  energyval.push_back( energy );
		}
		else if( ehits != ecalEndcapHits.end() ) { // add rechit only if available
		  const GlobalPoint & position = caloGeom->getPosition(*i);
		  double energy = ehits->energy();
		  double et = energy*position.perp()/position.mag();
		  //if( et == 0. && energy == 0. ) cout << "zero energy in cell " << endl; //debug
		  if( et < etLow || energy < eLow ) continue; //thresholds on e/et of cells

		  energypos.push_back( position );
		  energyval.push_back( energy );
		}
	  }
	}
	//cout << "newfn nhits " << energypos.size() << "  " << energyval.size() << endl;
  }
  else if( abs( objid ) == 13 ) { //muon
	//const reco::Muon* muon = &(muonHandle->at(objidx)); //get the muon which corresponds to obj
	//const reco::TrackRef& track = muon->innerTrack(); //get its innerTrack
	const reco::TrackRef& track = muonHandle->at(objidx).innerTrack(); //get innerTrack of muon obj
	TrackDetectorAssociator trackAssociator_;
	trackAssociator_.useDefaultPropagator();
	TrackDetMatchInfo info = trackAssociator_.associate(iEvent, iSetup, *(track.get()), muonparameters_);
	math::XYZPoint point = info.trkGlobPosAtEcal;
	pclu = GlobalPoint( point.x(), point.y(), point.z() ); //and finally, we have the point
	etaclus = pclu.eta();
	phiclus = pclu.phi();

	extRadius = musEcalExtCone_ ;
	if( TMath::Abs( etaclus ) <= 1.479 )
	  intRadius = musEcalVetoRadBarrel_;
	else
	  intRadius = musEcalVetoRadEndcap_;
	//debugging
	//cout << "muon " << muonHandle->at(objidx).p4().pt() << "  eta " << muonHandle->at(objidx).p4().eta() << endl;

	//set position vector--use towers for muons
	edm::Handle<CaloTowerCollection> caloTowers;
	iEvent.getByLabel(caloTowersInputTag_, caloTowers);
	for(CaloTowerCollection::const_iterator tower = caloTowers->begin();
		tower != caloTowers->end(); ++tower) {

	  const GlobalPoint & position = tower->emPosition();
	  //apply thresholds based on tower position (not mu)
	  if( TMath::Abs( position.eta() ) < 1.479 &&
		  (tower->emEt() < musEtMinBarrel_ || tower->emEnergy() < musEMinBarrel_) )
		continue;
	  else if( tower->emEt() < musEtMinEndcap_ || tower->emEnergy() < musEMinEndcap_ ) 
		continue;

	  double etadiff = position.eta() - etaclus;
	  double phidiff = deltaPhi( (double)position.phi(), phiclus );
	  //only keep if in outer cone
	  if( (etadiff*etadiff + phidiff*phidiff) > extRadius*extRadius ) continue;
	  
	  energypos.push_back( position );
	  energyval.push_back( tower->emEnergy() );
	  //debugging
	  //double energy = tower->emEnergy();
	  //double et = energy*position.perp()/position.mag();
	  //cout << "tower " << et << "  dr " << sqrt(etadiff*etadiff + phidiff*phidiff) << endl;
	  //if( abs(et - tower->emEt()) > 0.01 )
	  //cout << "et disagreement: from pos " << et << "  from tower " << tower->emEt() << endl;
	}
  } 
  else {
	cout << "HypIsoMaker : ERROR bad hyp_lx_id--not abs = 11 || 13\n\n";
	return 0;
  }

  if( energypos.size() != energyval.size() ) {
	cout << "HypIsoMaker : ERROR bad energy vector size\n\n";
	return 0;
  }

  //loop on vector which was filled above, do exclusion
  for( unsigned int i=0; i < energypos.size(); i++) {

	double eta = energypos[i].eta();		  
	double phi = energypos[i].phi();		  
	double etaDiff = eta - etaclus;		  
	double phiDiff= deltaPhi(phi,phiclus);

	double energy = energyval[i]; 
	double et = energy*energypos[i].perp()/energypos[i].mag();
	  
	//loop on exclusion leptons
	bool exclude = false;
	for( unsigned int k=0; k < excid.size(); k++ ) {
	  double excintRadius = 0;
	  double excoutRadius = 0;
	  GlobalPoint pclux(0,0,0);
	  //Take the SC position for exclusion els
	  if( TMath::Abs(excid.at(k)) == 11 ) { //exclusion electrons
		//get gsf electron from uniqueElectrons which corresponds to exc el
		const reco::Candidate* emObjectx = &(emObjectHandle->at(excidx.at(k))); 
		reco::SuperClusterRef scx = emObjectx->get<reco::SuperClusterRef>();
		math::XYZPoint theCaloPositionx = scx.get()->position();
		pclux = GlobalPoint( theCaloPositionx.x(), theCaloPositionx.y(), theCaloPositionx.z() );
		excoutRadius = elsEcalExtCone_;
		if( TMath::Abs(pclux.eta()) <= 1.479 ) excintRadius = elsEcalVetoRadBarrel_; //barrel
		else excintRadius = elsEcalVetoRadEndcap_; //endcap
	  }
	  else if( TMath::Abs(excid.at(k)) == 13 ) { //exclusion muons
		//get innerTrack of muon obj, then get its position at ecal
		const reco::TrackRef& track = muonHandle->at(excidx.at(k)).innerTrack();
		TrackDetectorAssociator trackAssociator_;
		trackAssociator_.useDefaultPropagator();
		TrackDetMatchInfo info = trackAssociator_.associate(iEvent, iSetup, *(track.get()), muonparameters_);
		math::XYZPoint point = info.trkGlobPosAtEcal;
		pclux = GlobalPoint( point.x(), point.y(), point.z() ); //and finally, we have the point
		excoutRadius = musEcalExtCone_;
		if( TMath::Abs(pclux.eta()) <= 1.479 ) excintRadius = musEcalVetoRadBarrel_; //barrel
		else excintRadius = musEcalVetoRadEndcap_; //endcap
	  }
	  else {
		cout << "HypIsoMaker : ERROR bad hyp_lx_id--not abs = 11 || 13\n\n";
		continue;
	  }
			
	  double etadiffx = pclux.eta() - eta;
	  double phidiffx = pclux.phi() - phi;
	  double drx = sqrt( etadiffx*etadiffx + phidiffx*phidiffx );

	  if( drx < excoutRadius && ( ( TMath::Abs(excid.at(k)) == 11 && TMath::Abs(etadiffx) < IsoJurassicWidth_ )
								  || drx < excintRadius ) ) {
		//cout << "excluding " << et << "  dr " << drx << "  etadiff " << etadiffx << endl;
		exclude = true;
		break;
	  }
	}

	//original exclusion region
	if( TMath::Abs(objid) == 11 && TMath::Abs(etaDiff) < IsoJurassicWidth_ ) continue;  // jurassic strip cut for els only

	if( sqrt( etaDiff*etaDiff + phiDiff*phiDiff ) < intRadius ) continue; // exclusion cone cut

	//check second hyp veto--only add to vetoSum if wasn't vetoed by primary hyp
	if( exclude ) {
	  if(returnEt) vetoSum += et;
	  else vetoSum += energy;
	  continue;
	}
  
	if(returnEt) energySum += et;
	else energySum += energy;
	
  } //end loop on 'energy vectors'

  //cout << "\t" << energySum << "\t" << etaclusx << "\t" << phiclusx << endl;
  //cout << "\t" << energySum << endl;
  if( recomputeEcalIso_ ) //only return excluded sum to subtract off default
	return energySum;
  else
	return vetoSum;
}


//instead of passing indicies, pass p4s of e or mu to exclude, track vector, and the id's to set the cone sizes
double HypIsoMaker::track_iso(LorentzVector prip4, float pri_d0, float pri_z0, int priid, vector<LorentzVector> excp4, vector<int> excid,
							  const vector<LorentzVector> *trksp4, const vector<float> *trks_d0, const vector<float> *trks_z0) {
  //return is the sum of the pt of all tracks in the range specified in config around the lepton, and exclude cone around both

  double isolation = 0;
  double vetopt = 0;
  double pri_mindr = 0;
  double exc_mindr = 0;
  if( TMath::Abs(priid) == 11 ) pri_mindr = trackIsoElsInRadius_;
  else if( TMath::Abs(priid) == 13 ) pri_mindr = trackIsoMusInRadius_;
  else cout << "HypIsoMaker: bad hyp_lx_id\n\n";

  for( unsigned int i=0; i < trksp4->size(); i++ ) {
	//cuts on track quality
	if( TMath::Abs( trks_d0->at(i) - pri_d0 ) > trackIsoMind0_
		|| TMath::Abs( trks_z0->at(i) - pri_z0 ) > trackIsoMinz0_
		|| trksp4->at(i).pt() < trackIsoMinPt_ )
	  continue;

	double dR1 = ROOT::Math::VectorUtil::DeltaR( prip4, trksp4->at(i) );
	bool exclude = false;
	for( unsigned int j=0; j < excp4.size(); j++ ) {
	  double dR2 = ROOT::Math::VectorUtil::DeltaR( excp4[j], trksp4->at(i) );
	  if( abs(excid.at(j)) == 11 ) exc_mindr = trackIsoElsInRadius_;
	  else if( abs(excid.at(j)) == 13 ) exc_mindr = trackIsoMusInRadius_;
	  else cout << "HypIsoMaker: bad hyp_lx_id\n\n";
  
	  if( dR2 < exc_mindr && dR1 < trackIsoExtRadius_ ) { //exclude this track because of second hyp
		exclude = true;
		break;
	  }
	}
	if( exclude ) {
	  vetopt += trksp4->at(i).pt();
	  continue;
	}

	if( dR1 < pri_mindr || dR1 > trackIsoExtRadius_ ) //primary exclusion cone
	  continue; 
	
	isolation += trksp4->at(i).pt();
  }

  if( recomputeTrckIso_ ) //only return excluded pt to subtract off default
	return isolation;
  else
	return vetopt;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HypIsoMaker);
