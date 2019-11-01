//-*- C++ -*- 
//
// Package:    ElectronMaker
// Class:      ElectronMaker
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.cc,v 1.89 2012/08/16 00:00:27 slava77 Exp $


//System include files
#include <memory>
#include <math.h>

//User include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CMS3/NtupleMaker/interface/plugins/ElectronMaker.h"
#include "CMS3/NtupleMaker/interface/plugins/MatchUtilities.h"
#include "CMS3/NtupleMaker/interface/plugins/MCUtilities.h"
#include "CMS3/NtupleMaker/interface/EgammaFiduciality.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Math/VectorUtil.h"
#include "TVector2.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

using namespace reco;
using namespace edm;
using namespace std;

typedef math::XYZPoint Point;
typedef Ref<edmNew::DetSetVector<SiStripCluster>,SiStripCluster > ClusterRef;
typedef Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > pixel_ClusterRef;

// constructors and destructor
ElectronMaker::ElectronMaker(const ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasPrefix")),
  year_(iConfig.getParameter<int>("year")),

  //beamSpot_tag_                (iConfig.getParameter<edm::InputTag>("beamSpotTag")),

  trksInputTag_                (iConfig.getParameter<edm::InputTag>("trksInputTag")),
  gsftracksInputTag_           (iConfig.getParameter<edm::InputTag>("gsftracksInputTag")),

  //cms2scsseeddetidInputTag_  (iConfig.getParameter<edm::InputTag> ("cms2scsseeddetidInputTag"     )),
  eidLHTag_                    (iConfig.getParameter<edm::InputTag>("eidLHTag")),

  ebReducedRecHitCollectionTag (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollectionTag")),
  eeReducedRecHitCollectionTag (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollectionTag")),
  esReducedRecHitCollectionTag (iConfig.getParameter<edm::InputTag>("esReducedRecHitCollectionTag")),
  //pfIsoCharged03InputTag       (iConfig.getParameter<edm::InputTag> ("pfIsoCharged03InputTag"   )),
  //pfIsoGamma03InputTag         (iConfig.getParameter<edm::InputTag> ("pfIsoGamma03InputTag"     )),
  //pfIsoNeutral03InputTag       (iConfig.getParameter<edm::InputTag> ("pfIsoNeutral03InputTag"   )),
  //pfIsoCharged04InputTag       (iConfig.getParameter<edm::InputTag> ("pfIsoCharged04InputTag"   )),
  //pfIsoGamma04InputTag         (iConfig.getParameter<edm::InputTag> ("pfIsoGamma04InputTag"     )),
  //pfIsoNeutral04InputTag       (iConfig.getParameter<edm::InputTag> ("pfIsoNeutral04InputTag"   )),
  rhoInputTag_                 (iConfig.getParameter<edm::InputTag>("rhoInputTag")),

  minAbsDist_                  (iConfig.getParameter<double>("minAbsDist")),
  minAbsDcot_                  (iConfig.getParameter<double>("minAbsDcot")),
  minSharedFractionOfHits_     (iConfig.getParameter<double>("minSharedFractionOfHits")),

  pfCandidates(nullptr)
{
  //get setup parameters

  electronsToken  = consumes<edm::View<pat::Electron>  >(iConfig.getParameter<edm::InputTag>("electronsInputTag"));

  miniIsoChgValueMapToken_   = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("miniIsoChgValueMap"));
  miniIsoAllValueMapToken_   = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("miniIsoAllValueMap"));

  beamSpotToken  = consumes<LorentzVector>(iConfig.getParameter<edm::InputTag>("beamSpotInputTag"));

  pfJetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  //pfCandsToken  = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandsInputTag"));
  vtxToken  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));

  bFieldToken  = consumes<float>(iConfig.getParameter<edm::InputTag>("bFieldInputTag"));

  recoConversionToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("recoConversionInputTag"));

  ebReducedRecHitCollection = mayConsume<EcalRecHitCollection>(ebReducedRecHitCollectionTag);
  eeReducedRecHitCollection = mayConsume<EcalRecHitCollection>(eeReducedRecHitCollectionTag);
  esReducedRecHitCollection = mayConsume<EcalRecHitCollection>(esReducedRecHitCollectionTag);

  produces<pat::ElectronCollection>().setBranchAlias(aliasprefix_);
}

ElectronMaker::~ElectronMaker(){}

void  ElectronMaker::beginRun(const edm::Run&, const EventSetup& es){}

void ElectronMaker::beginJob(){}

void ElectronMaker::endJob(){}

// ------------ method called to produce the data  ------------
void ElectronMaker::produce(Event& iEvent, const EventSetup& iSetup){
  auto result = std::make_unique<pat::ElectronCollection>();

  // uncertainties and corrections
  // somewhat complicated: see 
  // http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_1_2/doc/html/d5/d4b/GsfElectron_8h-source.html
  // note that if ecalEnergy == eSC depends on if further ecal corrections have been applied to the electron
  // after its construction

  // --- Get Input Collections --- //

  ///////////////
  // Electrons //
  ///////////////
  Handle<View<pat::Electron> > els_h;
  iEvent.getByToken(electronsToken, els_h);
  if (!els_h.isValid()) throw cms::Exception("ElectronMaker::produce: error getting electron collection from Event!");

  // Jets
  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);
  //Handle<GsfElectronCollection> els_coll_h;
  //iEvent.getByLabel(electronsInputTag_, els_coll_h);    

  //////////////
  // PF Cands //
  //////////////
  // iEvent.getByToken(pfCandsToken, packPfCand_h);
  // if (!packPfCand_h.isValid()) throw cms::Exception("ElectronMaker::produce: error getting packed pfcands from Event!");
  // pfCandidates  = packPfCand_h.product();

  ////////////
  // Vertex //
  ////////////
  iEvent.getByToken(vtxToken, vertexHandle);
  if (!vertexHandle.isValid()) throw cms::Exception("ElectronMaker::produce: error getting vertex collection from Event!");

  /////////////////
  // Conversions //
  /////////////////
  iEvent.getByToken(recoConversionToken, convs_h);
  if (!convs_h.isValid()) throw cms::Exception("ElectronMaker::produce: error getting conversion collection");

  ///////////////////////////
  // TransientTrackBuilder //
  ///////////////////////////
  //ESHandle<TransientTrackBuilder> theTTBuilder;
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  ////////////////////////////////////////////////
  // Get tools to get cluster shape information //
  ////////////////////////////////////////////////

  //////////////
  // Beamspot //
  //////////////
  Handle<LorentzVector> beamSpotH;
  iEvent.getByToken(beamSpotToken, beamSpotH);
  const Point beamSpot = beamSpotH.isValid() ? Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0, 0, 0);
  //Handle<reco::BeamSpot> beamspot_h;
  //iEvent.getByLabel(beamSpot_tag_, beamspot_h);
  //const reco::BeamSpot &beamSpotreco = *(beamspot_h.product()); 
  //if ( beamSpotreco.x0() == 1234567 ) ; // Avoid "unused variable" error while the function using this variable is inactive

  ///////////////////////
  // rho for isolation //
  ///////////////////////
  //edm::Handle<float> rhoIso_h;
  //iEvent.getByLabel(rhoInputTag_, rhoIso_h);
  //float rhoIso = *(rhoIso_h.product());

  edm::Handle<edm::ValueMap<float> > miniIsoChg_values;
  edm::Handle<edm::ValueMap<float> > miniIsoAll_values;

  // Corrected Isolation using NanoAOD
  iEvent.getByToken(miniIsoChgValueMapToken_, miniIsoChg_values);
  iEvent.getByToken(miniIsoAllValueMapToken_, miniIsoAll_values);

  //////////////////////////
  // get cms2scsseeddetid //
  //////////////////////////
  // InputTag cms2scsseeddetid_tag(cms2scsseeddetidInputTag_.label(),"scsdetIdSeed");
  // Handle<vector<int> > cms2scsseeddetid_h;
  // iEvent.getByLabel(cms2scsseeddetid_tag, cms2scsseeddetid_h); 
  // const vector<int> *cms2scsseeddetid = cms2scsseeddetid_h.product();

  //////////////////////////////
  // Get the ele<->PFCand map //
  //////////////////////////////
  // edm::Handle<edm::ValueMap<std::vector<reco::PFCandidateRef > > > eleToParticleBasedIsoMapHandle;
  // InputTag particleBase(string("particleBasedIsolation"),string("gedGsfElectrons"));  
  // iEvent.getByLabel(particleBase, eleToParticleBasedIsoMapHandle);    
  // edm::ValueMap<std::vector<reco::PFCandidateRef > >   eleToParticleBasedIsoMap =  *(eleToParticleBasedIsoMapHandle.product());

  // --- Fill --- //

  /////////////////////////
  // Loop Over Electrons //
  /////////////////////////

  size_t evt_nels = els_h->size(); result->reserve(evt_nels);
  size_t elsIndex = 0;
  for (View<pat::Electron>::const_iterator el = els_h->begin(); el != els_h->end(); el++, elsIndex++){
    //pat::Electron el_result(*(el->get<reco::GsfElectron const*>())); // Clone the gsfElectron. This is the single electron to be put into the resultant collection
    pat::Electron el_result(*el); // Clone the gsfElectron. This is the single electron to be put into the resultant collection

    ////////////////
    // References //
    ////////////////
    const GsfTrackRef el_track = el->gsfTrack(); // Embedded GSF Track for miniAOD
    const RefToBase<pat::Electron> gsfElRef = els_h->refAt(elsIndex);
    const TrackRef ctfTkRef = el->closestCtfTrackRef(); // Embedded CTF Track for miniAOD 

    ////////////
    // Vertex //
    ////////////
    const VertexCollection* vertexCollection = vertexHandle.product();
    VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
    int firstGoodVertexIdx = 0;
    for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx){
      // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
      // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
      if (  /*!vtx->isFake() &&*/ !(vtx->chi2()==0 && vtx->ndof()==0) &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0){
        firstGoodVertex = vtx;
        break;
      }
    }

    //////////////////////
    // Fiduciality Mask //
    //////////////////////
    int fiducialityMask = 0;  // the enum is in interface/EgammaFiduciality.h
    if (el->isEB()) fiducialityMask |= 1 << ISEB;
    if (el->isEBEEGap()) fiducialityMask |= 1 << ISEBEEGAP;
    if (el->isEE()) fiducialityMask |= 1 << ISEE;
    if (el->isEEGap()) fiducialityMask |= 1 << ISEEGAP;
    if (el->isEBEtaGap()) fiducialityMask |= 1 << ISEBETAGAP;
    if (el->isEBPhiGap()) fiducialityMask |= 1 << ISEBPHIGAP;
    if (el->isEEDeeGap()) fiducialityMask |= 1 << ISEEDEEGAP;
    if (el->isEERingGap()) fiducialityMask |= 1 << ISEERINGGAP;
    if (el->isGap()) fiducialityMask |= 1 << ISGAP;
    el_result.addUserInt("fid_mask", fiducialityMask);

    ///////////////////////////
    // Corrections & Seeding //
    ///////////////////////////
    int electronTypeMask = 0;
    if (el->isEcalEnergyCorrected()) electronTypeMask |= 1 << ISECALENERGYCORRECTED;
    if (el->trackerDrivenSeed()) electronTypeMask |= 1 << ISTRACKERDRIVEN;
    if (el->ecalDrivenSeed()) electronTypeMask |= 1 << ISECALDRIVEN;
    if (el->passingCutBasedPreselection()) electronTypeMask |= 1 << ISCUTPRESELECTED;
    // el->ecalDriven() == el->ecalDrivenSeed() && el->passingCutBasedPreselection()
    if (el->passingMvaPreselection()) electronTypeMask |= 1 << ISMVAPRESELECTED;
    //if ( el->isMomentumCorrected() ) electronTypeMask |= 1 << ISMOMENTUMCORRECTED; // Depricated in CMSSW_4_2x ( DataFormats/EgammaCandidates/interface/GsfElectron.h )
    el_result.addUserInt("type_mask", electronTypeMask);

    ///////////////////
    // Predefined ID //
    ///////////////////
    // el_result.addUserInt("category", classify(gsfElRef)); // this is the sani classification
    setMVAIdUserVariables(el, el_result, "ElectronMVAEstimatorRun2Fall17IsoV2", "Fall17V2_Iso");
    setMVAIdUserVariables(el, el_result, "ElectronMVAEstimatorRun2Fall17NoIsoV2", "Fall17V2_NoIso");
    setCutBasedIdUserVariables(el, el_result, "cutBasedElectronID-Fall17-94X-V2-veto", "Fall17V2_CutBased_Veto");
    setCutBasedIdUserVariables(el, el_result, "cutBasedElectronID-Fall17-94X-V2-loose", "Fall17V2_CutBased_Loose");
    setCutBasedIdUserVariables(el, el_result, "cutBasedElectronID-Fall17-94X-V2-medium", "Fall17V2_CutBased_Medium");
    setCutBasedIdUserVariables(el, el_result, "cutBasedElectronID-Fall17-94X-V2-tight", "Fall17V2_CutBased_Tight");

    //////////////
    // Electron //
    //////////////
    //el_result.addUserFloat("ecalEnergy", el->ecalEnergy());  // energy corrections and uncertainties
    //el_result.addUserFloat("ecalEnergyError", el->ecalEnergyError());
    el_result.addUserFloat("ecalEnergy", el->correctedEcalEnergy());  // energy corrections and uncertainties
    el_result.addUserFloat("ecalEnergyError", el->correctedEcalEnergyError());
    el_result.addUserFloat("trackMomentumError", el->trackMomentumError());

    //-- Scale and smearing corrections are now stored in the miniAOD https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    float uncorrected_pt = el->pt();
    float uncorrected_energy = el->energy();
    el_result.addUserFloat("scale_smear_corr", el->userFloat("ecalTrkEnergyPostCorr") / uncorrected_energy); // get scale/smear correction factor directly from miniAOD

    // the p4 of the electron is the uncorrected one
    el_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, el->eta(), el->phi(), el->mass()));

    //get all scale uncertainties and their breakdown
    el_result.addUserFloat("scale_smear_corr_scale_totalUp", el->userFloat("energyScaleUp") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_scale_statUp", el->userFloat("energyScaleStatUp") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_scale_systUp", el->userFloat("energyScaleSystUp") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_scale_gainUp", el->userFloat("energyScaleGainUp") / uncorrected_energy);

    el_result.addUserFloat("scale_smear_corr_scale_totalDn", el->userFloat("energyScaleDown") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_scale_statDn", el->userFloat("energyScaleStatDown") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_scale_systDn", el->userFloat("energyScaleSystDown") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_scale_gainDn", el->userFloat("energyScaleGainDown") / uncorrected_energy);

    //get all smearing uncertainties and their breakdown
    el_result.addUserFloat("scale_smear_corr_smear_totalUp", el->userFloat("energySigmaUp") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_smear_rhoUp", el->userFloat("energySigmaRhoUp") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_smear_phiUp", el->userFloat("energySigmaPhiUp") / uncorrected_energy);

    el_result.addUserFloat("scale_smear_corr_smear_totalDn", el->userFloat("energySigmaDown") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_smear_rhoDn", el->userFloat("energySigmaRhoDown") / uncorrected_energy);
    el_result.addUserFloat("scale_smear_corr_smear_phiDn", el->userFloat("energySigmaPhiDown") / uncorrected_energy);

    /////////////
    // Vectors //
    /////////////
    el_result.addUserFloat("trk_pt", el_track->pt());
    el_result.addUserFloat("trk_eta", el_track->eta());
    el_result.addUserFloat("trk_phi", el_track->phi());

    math::XYZVectorF p3In  = el->trackMomentumAtVtx();
    el_result.addUserFloat("trk_atvtx_pt", p3In.R());
    el_result.addUserFloat("trk_atvtx_eta", p3In.Eta());
    el_result.addUserFloat("trk_atvtx_phi", p3In.Phi());

    math::XYZVectorF p3Out = el->trackMomentumOut();
    el_result.addUserFloat("trk_out_pt", p3Out.R());
    el_result.addUserFloat("trk_out_eta", p3Out.Eta());
    el_result.addUserFloat("trk_out_phi", p3Out.Phi());

    el_result.addUserFloat("vtx_x", el->vx());
    el_result.addUserFloat("vtx_y", el->vy());
    el_result.addUserFloat("vtx_z", el->vz());

    el_result.addUserFloat("trk_vtx_x", el_track->vx());
    el_result.addUserFloat("trk_vtx_y", el_track->vy());
    el_result.addUserFloat("trk_vtx_z", el_track->vz());

    ///////////////
    // Isolation //
    ///////////////
    el_result.addUserFloat("hcalDepth1TowerSumEt", el->dr03HcalDepth1TowerSumEt());
    el_result.addUserFloat("ecalIso", el->dr03EcalRecHitSumEt());
    el_result.addUserFloat("hcalIso", el->dr03HcalTowerSumEt());
    el_result.addUserFloat("tkIso", el->dr03TkSumPt());

    el_result.addUserFloat("ecalIso04", el->dr04EcalRecHitSumEt());
    el_result.addUserFloat("hcalIso04", el->dr04HcalTowerSumEt());
    el_result.addUserFloat("tkIso04", el->dr04TkSumPt());

    //////////////////
    // PF Isolation //
    //////////////////
    GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    el_result.addUserFloat("pfChargedHadronIso", pfIso.sumChargedHadronPt);
    el_result.addUserFloat("pfNeutralHadronIso", pfIso.sumNeutralHadronEt);
    el_result.addUserFloat("pfPhotonIso", pfIso.sumPhotonEt);
    el_result.addUserFloat("pfPUIso", pfIso.sumPUPt);

    //////////////////
    // Supercluster //
    //////////////////
    el_result.addUserFloat("etaSC", el->superCluster()->eta());
    el_result.addUserFloat("phiSC", el->superCluster()->phi());
    el_result.addUserFloat("eSC", el->superCluster()->energy());
    el_result.addUserFloat("eSCRaw", el->superCluster()->rawEnergy());
    el_result.addUserFloat("eSCPresh", el->superCluster()->preshowerEnergy());
    el_result.addUserFloat("sigmaIEtaIEta", el->sigmaIetaIeta());
    el_result.addUserFloat("etaSCwidth", el->superCluster()->etaWidth());
    el_result.addUserFloat("phiSCwidth", el->superCluster()->phiWidth());

    // We used to make these using the cluster tools, but now we can take them directly from RECO electron
    el_result.addUserFloat("sigmaIPhiIPhi", el->sigmaIphiIphi());

    // Take these directly from the PAT electron of the miniAOD
    el_result.addUserFloat("sigmaIPhiIPhi_full5x5", el->full5x5_sigmaIphiIphi());
    el_result.addUserFloat("sigmaEtaEta_full5x5", el->full5x5_sigmaEtaEta());
    el_result.addUserFloat("sigmaIEtaIEta_full5x5", el->full5x5_sigmaIetaIeta());
    el_result.addUserFloat("r9_full5x5", el->full5x5_r9());
    el_result.addUserFloat("r9", el->r9());
    el_result.addUserFloat("e1x5_full5x5", el->full5x5_e1x5());
    el_result.addUserFloat("e5x5_full5x5", el->full5x5_e5x5());
    el_result.addUserFloat("e2x5Max_full5x5", el->full5x5_e2x5Max());

    ///////////////////////////////////////////////////////
    // Get cluster info that is not stored in the object //
    ///////////////////////////////////////////////////////
    //
    // This is a fix for accessing SC information in reminiAOD_V2
    //
    size_t numberOfClusters =  el->superCluster()->clusters().size();
    bool missing_clusters = false;
    if (numberOfClusters > 0) missing_clusters = !el->superCluster()->clusters()[numberOfClusters-1].isAvailable();

    size_t numberOfPSClusters =  el->superCluster()->preshowerClusters().size();
    bool missing_PSclusters = false;
    if (numberOfPSClusters > 0) missing_PSclusters = !el->superCluster()->preshowerClusters()[numberOfPSClusters-1].isAvailable();

    if (!(missing_clusters || missing_PSclusters) && (el->p4().pt() > 5.)){
      // The commented ones are already available above! Keeping it here for reference

      //	els_scRawEnergy             = el->superCluster()->rawEnergy();
      //	els_scCalibratedEnergy      = el->superCluster()->energy();
      //	els_scPreshowerEnergy       = el->superCluster()->preshowerEnergy();
      //	els_scEta                   = el->superCluster()->position().Eta();
      //	els_scPhi                   = el->superCluster()->position().Phi();
      //	els_scPhiWidth              = el->superCluster()->phiWidth();
      //	els_scEtaWidth              = el->superCluster()->etaWidth();
      //	els_scSeedRawEnergy         = el->superCluster()->seed()->energy();
      //	els_scSeedCalibratedEnergy  = el->superCluster()->seed()->energy();
      //  el_result.addUserFloat("scSeedE5x5", clusterTools_->e5x5(*(el->superCluster()->seed())));
      //	el_result.addUserFloat("scSeedE2x5max", clusterTools_->e2x5Max(*(el->superCluster()->seed())));
      //  el_result.addUserFloat("scSeedSigmaIetaIeta", see);
      //  el_result.addUserFloat("scSeedSigmaIphiIphi", spp); 


      // The one below is kept for historical reasons
      el_result.addUserFloat("SC_seed_energy", el->superCluster()->seed()->energy());
      el_result.addUserFloat("SC_seed_eta", el->superCluster()->seed()->eta());
    }
    else{
      el_result.addUserFloat("SC_seed_energy", -999.);
      el_result.addUserFloat("SC_seed_eta", -999.);
    }
    //
    //            //
    //            const BasicCluster&  clRef              = *(el->superCluster()->seed());
    //            const vector<float>& covs               = clusterTools_->covariances(clRef);                         // get the covariances computed in 5x5 around the seed
    //            const vector<float>& lcovs              = clusterTools_->localCovariances(clRef);                    // get the local covariances computed in a 5x5 around the seed
    //            const vector<float>  localCovariancesSC = clusterTools_->scLocalCovariances(*(el->superCluster()));  // get the local covariances computed using all crystals in the SC
    //
    //            //
    //get from RECO            el_result.addUserFloat("sigmaIPhiIPhi   ",  isfinite(lcovs[2])              ? lcovs[2] > 0               ? sqrt(lcovs[2]) : -1 * sqrt(-1 * lcovs[2])                             : -9999. );
    //
    //            //


    ////////
    // ID //
    ////////
    //el_result.addUserFloat("hOverE                        ",  el->hadronicOverEm()                 );
    el_result.addUserFloat("hOverE", el->hcalOverEcal());
    el_result.addUserFloat("full5x5_hOverE", el->full5x5_hcalOverEcal());
    el_result.addUserFloat("eOverPIn", el->eSuperClusterOverP());
    el_result.addUserFloat("eOverPOut", el->eEleClusterOverPout());
    el_result.addUserFloat("fbrem", el->fbrem());

    el_result.addUserFloat("dEtaIn", el->deltaEtaSuperClusterTrackAtVtx());
    el_result.addUserFloat("dEtaOut", el->deltaEtaSeedClusterTrackAtCalo());
    el_result.addUserFloat("dPhiIn", el->deltaPhiSuperClusterTrackAtVtx());
    el_result.addUserFloat("dPhiOut", el->deltaPhiSeedClusterTrackAtCalo());
    el_result.addUserFloat("isGsfCtfScPixChargeConsistent", el->isGsfCtfScPixChargeConsistent());

    ////////////
    // Tracks //
    ////////////
    float pt       = el_track->pt();
    float p        = el_track->p();
    float q        = el_track->charge();
    float pz       = el_track->pz();
    float trkpterr = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0, 0))+2*pt*p/q*pz*(el_track->covariance(0, 1))+ pz*pz*(el_track->covariance(1, 1))) : -9999.;
    el_result.addUserInt("charge", el->charge());
    el_result.addUserInt("trk_charge", el_track->charge());
    el_result.addUserInt("SCcharge", el->scPixCharge());
    el_result.addUserFloat("chi2", el_track->chi2());
    el_result.addUserFloat("ndof", el_track->ndof());
    el_result.addUserFloat("d0Err", el_track->d0Error());
    el_result.addUserFloat("z0Err", el_track->dzError());
    el_result.addUserFloat("ptErr", trkpterr);
    el_result.addUserFloat("ptErrGsf", el_track->ptError());
    el_result.addUserFloat("validHits", el_track->numberOfValidHits());
    el_result.addUserFloat("lostHits", el_track->numberOfLostHits());
    if (firstGoodVertex!=vertexCollection->end()){
      el_result.addUserFloat("dxyPV", el_track->dxy(firstGoodVertex->position()));
      el_result.addUserFloat("dzPV", el_track->dz(firstGoodVertex->position()));
    }
    else{
      el_result.addUserFloat("dxyPV", -999.);
      el_result.addUserFloat("dzPV", -999.);
    }

    el_result.addUserFloat("dz_firstPV", el_track->dz((vertexCollection->begin())->position()));
    el_result.addUserFloat("dxy_firstPV", el_track->dxy((vertexCollection->begin())->position()));

    /////////
    // CTF //
    /////////

    if (ctfTkRef.isNonnull()){
      //el_result.addUserFloat("trkshFrac", static_cast<float>( el->shFracInnerHits() )                                  );
      el_result.addUserFloat("trkshFrac", static_cast<float>(el->ctfGsfOverlap()));
      el_result.addUserFloat("trkdr", deltaR(el_track->eta(), el_track->phi(), ctfTkRef->eta(), ctfTkRef->phi()));
      el_result.addUserFloat("ckf_chi2", ctfTkRef->chi2());
      el_result.addUserFloat("ckf_ndof", ctfTkRef->ndof());
      el_result.addUserInt("ckf_laywithmeas", ctfTkRef->hitPattern().trackerLayersWithMeasurement());
      el_result.addUserInt("ckf_charge", ctfTkRef->charge());
    }
    else{
      el_result.addUserFloat("trkshFrac", -999.);
      el_result.addUserFloat("trkdr", -999.);
      el_result.addUserFloat("ckf_chi2", -999.);
      el_result.addUserFloat("ckf_ndof", -999.);
      el_result.addUserFloat("ckf_laywithmeas", -999.);
      el_result.addUserFloat("ckf_charge", -999);
    }


    //        ////////////////////
    //        // Regular Vertex //
    //        ////////////////////        
    //        TransientTrack tt = theTTBuilder->build(el->gsfTrack());
    //    
    //        if ( firstGoodVertex!=vertexCollection->end() ) {
    //            Measurement1D ip3D_regular = IPTools::absoluteImpactParameter3D(tt, *firstGoodVertex).second;
    //            //
    //            el_result.addUserFloat("ip3d", ip3D_regular.value() );
    //            el_result.addUserFloat("ip3derr", ip3D_regular.error() );
    //        } else {
    //            //
    //            el_result.addUserFloat("ip3d", -999. );
    //            el_result.addUserFloat("ip3derr", -999. );
    //        }


    //Impact Parameters
    el_result.addUserFloat("ip3d", el->dB(pat::Electron::PV3D));
    el_result.addUserFloat("ip3derr", el->edB(pat::Electron::PV3D));
    el_result.addUserFloat("ip2d", el->dB(pat::Electron::PV2D));
    el_result.addUserFloat("ip2derr", el->edB(pat::Electron::PV2D));
    //el_result.addUserFloat("ip3d", el->ip3d()); // miniAOD


    /////////////////
    // Hit Pattern //
    /////////////////
    // Redesign according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingHitPatternRedesign
    const HitPattern& pattern = el_track->hitPattern();
    //const HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
    //const HitPattern& p_outer = el_track->trackerExpectedHitsOuter();
    el_result.addUserInt("n_missing_inner_hits", pattern.numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    el_result.addUserInt("n_missing_outer_hits", pattern.numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
    el_result.addUserInt("n_valid_pixel_hits", pattern.numberOfValidPixelHits());
    el_result.addUserInt("n_lost_pixel_hits", pattern.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)); // Not sure about this. Could be MISSING_INNER_HITS instead.
    el_result.addUserInt("n_tracker_layers", pattern.trackerLayersWithMeasurement());
    el_result.addUserInt("n_tracker_layers_3D", pattern.pixelLayersWithMeasurement() + pattern.numberOfValidStripLayersWithMonoAndStereo());
    el_result.addUserInt("n_tracker_layers_lost", pattern.trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));

    /////////////////
    // Conversions //
    /////////////////
    reco::ConversionRef conv_ref = ConversionTools::matchedConversion(*el, convs_h, beamSpot);
    float vertexFitProbability = -1;
    if (!conv_ref.isNull()) {
      const reco::Vertex &vtx = conv_ref.get()->conversionVertex();
      if (vtx.isValid()) vertexFitProbability = TMath::Prob(vtx.chi2(), vtx.ndof());
    }
    el_result.addUserFloat("conv_vtx_prob", vertexFitProbability);
    //cout<<"Found electron with pt eta phi "<<el->p4().pt() <<" "<< el->p4().eta() <<" "<< el->p4().phi()<<" and vertexFitProbability "<<vertexFitProbability<<endl;

    //////////////////////////////////////////////
    // Flag For Vertex Fit Conversion Rejection //
    //////////////////////////////////////////////
    el_result.addUserInt("conv_vtx_flag", static_cast<int>(!el->passConversionVeto())); // PAT variable: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc#467

    //////////////////////
    // genMatch miniAOD //
    //////////////////////
    LorentzVector mc_p4(0, 0, 0, 0);
    const reco::GenParticle* gen = el->genParticle();
    if (gen){
      mc_p4 = gen->p4();
      el_result.addUserFloat("mc_patMatch_id", gen->pdgId());
      el_result.addUserFloat("mc_patMatch_pt", mc_p4.pt());
      el_result.addUserFloat("mc_patMatch_eta", mc_p4.eta());
      el_result.addUserFloat("mc_patMatch_phi", mc_p4.phi());
      el_result.addUserFloat("mc_patMatch_mass", mc_p4.mass());
      el_result.addUserFloat("mc_patMatch_dr", ROOT::Math::VectorUtil::DeltaR(gen->p4(), el->p4()));
    }
    else{
      el_result.addUserFloat("mc_patMatch_id", -999);
      el_result.addUserFloat("mc_patMatch_pt", 0);
      el_result.addUserFloat("mc_patMatch_eta", 0);
      el_result.addUserFloat("mc_patMatch_phi", 0);
      el_result.addUserFloat("mc_patMatch_mass", 0);
      el_result.addUserFloat("mc_patMatch_dr", -999.);
    }

    //////////////////////
    // mini-isolation   //
    //////////////////////

    // float minichiso     = 0.;
    // float mininhiso     = 0.;
    // float miniemiso     = 0.;
    // float minidbiso     = 0.;
    // elMiniIso(el, true, 0.0, minichiso, mininhiso, miniemiso, minidbiso);
    // el_result.addUserFloat("miniIso_uncor   ",  minichiso + mininhiso + miniemiso );
    // el_result.addUserFloat("miniIso_ch      ",  minichiso );
    // el_result.addUserFloat("miniIso_nh      ",  mininhiso );
    // el_result.addUserFloat("miniIso_em      ",  miniemiso );
    // el_result.addUserFloat("miniIso_db      ",  minidbiso );
    auto el2 = el->clone();
    auto miniiso = el2->miniPFIsolation();
    el_result.addUserFloat("miniIso_uncor", miniiso.chargedHadronIso() + miniiso.neutralHadronIso() + miniiso.photonIso());
    el_result.addUserFloat("miniIso_ch", miniiso.chargedHadronIso());
    el_result.addUserFloat("miniIso_nh", miniiso.neutralHadronIso());
    el_result.addUserFloat("miniIso_em", miniiso.photonIso());
    el_result.addUserFloat("miniIso_db", miniiso.puChargedHadronIso());
    delete el2;

    ///////////////////////////
    // PFCluster isolation   //
    ///////////////////////////
    el_result.addUserFloat("ecalPFClusterIso", el->ecalPFClusterIso());
    el_result.addUserFloat("hcalPFClusterIso", el->hcalPFClusterIso());

    result->emplace_back(el_result);
  } // end Loop on Electrons

  // Put the electron collection into the event
  iEvent.put(std::move(result));
}

//----------------------------------------------------------------------------
// Electron Id classification function (a flag for the Sani type class)
//----------------------------------------------------------------------------
int ElectronMaker::classify(const RefToBase<pat::Electron> &electron) {
  double eOverP = electron->eSuperClusterOverP();
  double fbrem = electron->fbrem();

  int cat;
  if ((electron->isEB() && fbrem<0.06) || (electron->isEE() && fbrem<0.1)) cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) cat=0;
  else cat=2;
  return cat;
}

//little labour-saving function to get the reference to the ValueMap
template<typename T> const ValueMap<T>& ElectronMaker::getValueMap(const Event& iEvent, InputTag& inputTag){
  Handle<ValueMap<T> > handle;
  iEvent.getByLabel(inputTag, handle);
  return *(handle.product());
}

double ElectronMaker::electronIsoValuePF(const GsfElectron& el, const Vertex& vtx, float coner, float minptn, float dzcut, float footprintdr, float gammastripveto, float elestripveto, int filterId){

  float pfciso = 0.;
  float pfniso = 0.;
  float pffootprint = 0.;
  float pfjurveto = 0.;
  float pfjurvetoq = 0.;

  //TrackRef siTrack     = el.closestCtfTrackRef();
  TrackRef siTrack     = el.closestTrack();
  GsfTrackRef gsfTrack = el.gsfTrack();

  if (gsfTrack.isNull() && siTrack.isNull()) return -9999.;

  float eldz = gsfTrack.isNonnull() ? gsfTrack->dz(vtx.position()) : siTrack->dz(vtx.position());
  float eleta = el.eta();

  for (PFCandidateCollection::const_iterator pf=pfCand_h->begin(); pf<pfCand_h->end(); ++pf){

    float pfeta = pf->eta();
    float dR = deltaR(pfeta, pf->phi(), eleta, el.phi());
    if (dR>coner) continue;

    float deta = fabs(pfeta - eleta);
    int pfid = abs(pf->pdgId());
    float pfpt = pf->pt();

    if (filterId!=0 && filterId!=pfid) continue;

    if (pf->charge()==0) {
      //neutrals
      if (pfpt>minptn) {
        pfniso+=pfpt;
        if (dR<footprintdr && pfid==130) pffootprint+=pfpt;
        if (deta<gammastripveto && pfid==22)  pfjurveto+=pfpt;
      }
    }
    else {
      //charged  
      //avoid double counting of electron itself
      //if either the gsf or the ctf track are shared with the candidate, skip it
      const TrackRef pfTrack  = pf->trackRef();
      if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
      //below pfid==1 is commented out: in some cases the pfCand has a gsf even if it is not an electron... this is to improve the sync with MIT
      if (/*pfid==11 &&*/ pf->gsfTrackRef().isNonnull()) {
        if (gsfTrack.isNonnull() && gsfTrack.key()==pf->gsfTrackRef().key()) continue;
      }
      //check electrons with gsf track
      if (pfid==11 && pf->gsfTrackRef().isNonnull()) {
        if (fabs(pf->gsfTrackRef()->dz(vtx.position()) - eldz)<dzcut) {//dz cut
          pfciso+=pfpt;
          if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
        }
        continue;//and avoid double counting
      }
      //then check anything that has a ctf track
      if (pfTrack.isNonnull()) {//charged (with a ctf track)
        if (fabs(pfTrack->dz(vtx.position()) - eldz)<dzcut) {//dz cut
          pfciso+=pfpt;
          if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
        }
      }
    }
  }
  return pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq;
}

void ElectronMaker::elIsoCustomCone(edm::View<pat::Electron>::const_iterator& el, float dr, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso){
  chiso     = 0.;
  nhiso     = 0.;
  emiso     = 0.;
  dbiso     = 0.;
  float deadcone_ch = 0.;
  float deadcone_pu = 0.;
  float deadcone_ph = 0.;

  double phi = el->p4().Phi();
  double eta = el->p4().Eta();
  double pi = M_PI;

  // veto cones only in the endcap for electrons
  if (useVetoCones && fabs(el->superCluster()->eta()) > 1.479) {
    deadcone_ch = 0.015;
    deadcone_pu = 0.015;
    deadcone_ph = 0.08;
  }
  for (pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++) {
    float id = pf_it->pdgId();
    if (fabs(id) != 211 && fabs(id) != 130 && fabs(id) != 22) continue;

    double deltaPhi = phi-pf_it->p4().Phi();
    if (deltaPhi > pi) deltaPhi -= 2.0*pi;
    else if (deltaPhi <= -pi) deltaPhi += 2.0*pi;
    deltaPhi = fabs(deltaPhi);
    if (deltaPhi > dr) continue;
    double deltaEta = fabs(pf_it->p4().Eta()-eta);
    if (deltaEta > dr) continue;
    double thisDR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

    if (thisDR>dr) continue;
    float pt = pf_it->p4().pt();
    if (fabs(id)==211) {
      if (pf_it->fromPV() > 1 && (!useVetoCones || thisDR > deadcone_ch)) chiso+=pt;
      else if ((pf_it->fromPV() <= 1) && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) dbiso+=pt;
    }
    if (fabs(id)==130 && (pt > ptthresh)) nhiso+=pt;
    if (fabs(id)==22 && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_ph)) emiso+=pt;
  }

  return;
}

void ElectronMaker::elMiniIso(edm::View<pat::Electron>::const_iterator& el, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float &dbiso){
  float pt = el->p4().pt();
  float dr = 0.2;
  if (pt>50) dr = 10./pt;
  if (pt>200) dr = 0.05;
  elIsoCustomCone(el, dr, useVetoCones, ptthresh, chiso, nhiso, emiso, dbiso);
  return;
}

void ElectronMaker::setMVAIdUserVariables(edm::View<pat::Electron>::const_iterator& el, pat::Electron& el_result, std::string const& id_name, std::string const& id_identifier) const{
  if (el->hasUserFloat(id_name+"RawValues")){
    el_result.addUserFloat("id_"+id_identifier+"_RawVal", el->userFloat(id_name+"RawValues"));
    el_result.addUserInt("id_"+id_identifier+"_Cat", el->userFloat(id_name+"Categories"));
  }
  else throw cms::Exception("ElectronMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}
void ElectronMaker::setCutBasedIdUserVariables(edm::View<pat::Electron>::const_iterator& el, pat::Electron& el_result, std::string const& id_name, std::string const& id_identifier) const{
  if (el->isElectronIDAvailable(id_name)){
    el_result.addUserInt("id_"+id_identifier+"_Bits", el->userInt(id_name));
  }
  else throw cms::Exception("ElectronMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}


//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
