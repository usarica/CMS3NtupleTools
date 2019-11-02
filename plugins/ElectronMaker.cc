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
  size_t electronIndex = 0;
  for (View<pat::Electron>::const_iterator el = els_h->begin(); el != els_h->end(); el++, electronIndex++){
    pat::Electron electron_result(*el); // Clone the gsfElectron. This is the single electron to be put into the resultant collection

    ////////////////
    // References //
    ////////////////
    const GsfTrackRef gsfTrack = el->gsfTrack(); // Embedded GSF Track for miniAOD
    const RefToBase<pat::Electron> gsfElectron = els_h->refAt(electronIndex);
    const TrackRef ctfTrack = el->closestCtfTrackRef(); // Embedded CTF Track for miniAOD 

    ////////////
    // Vertex //
    ////////////
    const VertexCollection* vertexCollection = vertexHandle.product();
    VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
    int firstGoodVertexIdx = 0;
    for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx){
      // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
      // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
      if (/*!vtx->isFake() &&*/ !(vtx->chi2()==0 && vtx->ndof()==0) &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0){
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
    electron_result.addUserInt("fid_mask", fiducialityMask);

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
    electron_result.addUserInt("type_mask", electronTypeMask);

    ///////////////////
    // Predefined ID //
    ///////////////////
    // electron_result.addUserInt("category", classify(gsfElectron)); // this is the sani classification
    setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Fall17IsoV2", "Fall17V2_Iso");
    setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Fall17NoIsoV2", "Fall17V2_NoIso");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-veto", "Fall17V2_Veto");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-loose", "Fall17V2_Loose");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-medium", "Fall17V2_Medium");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-tight", "Fall17V2_Tight");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-veto", "Fall17V1_Veto");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-loose", "Fall17V1_Loose");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-medium", "Fall17V1_Medium");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-tight", "Fall17V1_Tight");

    //////////////
    // Electron //
    //////////////
    //electron_result.addUserFloat("ecalEnergy", el->ecalEnergy());  // energy corrections and uncertainties
    //electron_result.addUserFloat("ecalEnergyError", el->ecalEnergyError());
    electron_result.addUserFloat("ecalEnergy", el->correctedEcalEnergy());  // energy corrections and uncertainties
    electron_result.addUserFloat("ecalEnergyError", el->correctedEcalEnergyError());
    electron_result.addUserFloat("trackMomentumError", el->trackMomentumError());

    //-- Scale and smearing corrections are now stored in the miniAOD https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    float uncorrected_pt = el->pt();
    float uncorrected_mass = el->mass();
    float uncorrected_energy = el->energy();
    electron_result.addUserFloat("scale_smear_corr", el->userFloat("ecalTrkEnergyPostCorr") / uncorrected_energy); // get scale/smear correction factor directly from miniAOD

    // the p4 of the electron is the uncorrected one
    electron_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, el->eta(), el->phi(), uncorrected_mass));

    // get scale uncertainties and their breakdown
    electron_result.addUserFloat("scale_smear_corr_scale_totalUp", el->userFloat("energyScaleUp") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_scale_statUp", el->userFloat("energyScaleStatUp") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_scale_systUp", el->userFloat("energyScaleSystUp") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_scale_gainUp", el->userFloat("energyScaleGainUp") / uncorrected_energy);

    electron_result.addUserFloat("scale_smear_corr_scale_totalDn", el->userFloat("energyScaleDown") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_scale_statDn", el->userFloat("energyScaleStatDown") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_scale_systDn", el->userFloat("energyScaleSystDown") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_scale_gainDn", el->userFloat("energyScaleGainDown") / uncorrected_energy);

    // get smearing uncertainties and their breakdown
    electron_result.addUserFloat("scale_smear_corr_smear_totalUp", el->userFloat("energySigmaUp") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_smear_rhoUp", el->userFloat("energySigmaRhoUp") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_smear_phiUp", el->userFloat("energySigmaPhiUp") / uncorrected_energy);

    electron_result.addUserFloat("scale_smear_corr_smear_totalDn", el->userFloat("energySigmaDown") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_smear_rhoDn", el->userFloat("energySigmaRhoDown") / uncorrected_energy);
    electron_result.addUserFloat("scale_smear_corr_smear_phiDn", el->userFloat("energySigmaPhiDown") / uncorrected_energy);

    /////////////
    // Vectors //
    /////////////
    electron_result.addUserFloat("trk_pt", gsfTrack->pt());
    electron_result.addUserFloat("trk_eta", gsfTrack->eta());
    electron_result.addUserFloat("trk_phi", gsfTrack->phi());

    math::XYZVectorF p3In  = el->trackMomentumAtVtx();
    electron_result.addUserFloat("trk_atvtx_pt", p3In.R());
    electron_result.addUserFloat("trk_atvtx_eta", p3In.Eta());
    electron_result.addUserFloat("trk_atvtx_phi", p3In.Phi());

    math::XYZVectorF p3Out = el->trackMomentumOut();
    electron_result.addUserFloat("trk_out_pt", p3Out.R());
    electron_result.addUserFloat("trk_out_eta", p3Out.Eta());
    electron_result.addUserFloat("trk_out_phi", p3Out.Phi());

    electron_result.addUserFloat("vtx_x", el->vx());
    electron_result.addUserFloat("vtx_y", el->vy());
    electron_result.addUserFloat("vtx_z", el->vz());

    electron_result.addUserFloat("trk_vtx_x", gsfTrack->vx());
    electron_result.addUserFloat("trk_vtx_y", gsfTrack->vy());
    electron_result.addUserFloat("trk_vtx_z", gsfTrack->vz());

    ///////////////
    // Isolation //
    ///////////////
    electron_result.addUserFloat("hcalDepth1TowerSumEt", el->dr03HcalDepth1TowerSumEt());
    electron_result.addUserFloat("ecalIso", el->dr03EcalRecHitSumEt());
    electron_result.addUserFloat("hcalIso", el->dr03HcalTowerSumEt());
    electron_result.addUserFloat("tkIso", el->dr03TkSumPt());

    electron_result.addUserFloat("ecalIso04", el->dr04EcalRecHitSumEt());
    electron_result.addUserFloat("hcalIso04", el->dr04HcalTowerSumEt());
    electron_result.addUserFloat("tkIso04", el->dr04TkSumPt());

    //////////////////
    // PF Isolation //
    //////////////////
    GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    electron_result.addUserFloat("pfChargedHadronIso", pfIso.sumChargedHadronPt);
    electron_result.addUserFloat("pfNeutralHadronIso", pfIso.sumNeutralHadronEt);
    electron_result.addUserFloat("pfPhotonIso", pfIso.sumPhotonEt);
    electron_result.addUserFloat("pfPUIso", pfIso.sumPUPt);

    //////////////////
    // Supercluster //
    //////////////////
    electron_result.addUserFloat("etaSC", el->superCluster()->eta());
    electron_result.addUserFloat("phiSC", el->superCluster()->phi());
    electron_result.addUserFloat("energySC", el->superCluster()->energy());
    electron_result.addUserFloat("rawEnergySC", el->superCluster()->rawEnergy());
    electron_result.addUserFloat("preshowerEnergySC", el->superCluster()->preshowerEnergy());
    electron_result.addUserFloat("sigmaIEtaIEta", el->sigmaIetaIeta());
    electron_result.addUserFloat("etaWidthSC", el->superCluster()->etaWidth());
    electron_result.addUserFloat("phiWidthSC", el->superCluster()->phiWidth());

    // We used to make these using the cluster tools, but now we can take them directly from RECO electron
    electron_result.addUserFloat("sigmaIPhiIPhi", el->sigmaIphiIphi());

    // Take these directly from the PAT electron of the miniAOD
    electron_result.addUserFloat("full5x5_sigmaIPhiIPhi", el->full5x5_sigmaIphiIphi());
    electron_result.addUserFloat("full5x5_sigmaEtaEta", el->full5x5_sigmaEtaEta());
    electron_result.addUserFloat("full5x5_sigmaIEtaIEta", el->full5x5_sigmaIetaIeta());
    electron_result.addUserFloat("full5x5_r9", el->full5x5_r9());
    electron_result.addUserFloat("r9", el->r9());
    electron_result.addUserFloat("full5x5_e1x5", el->full5x5_e1x5());
    electron_result.addUserFloat("full5x5_e5x5", el->full5x5_e5x5());
    electron_result.addUserFloat("full5x5_e2x5Max", el->full5x5_e2x5Max());

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
      //  electron_result.addUserFloat("scSeedE5x5", clusterTools_->e5x5(*(el->superCluster()->seed())));
      //	electron_result.addUserFloat("scSeedE2x5max", clusterTools_->e2x5Max(*(el->superCluster()->seed())));
      //  electron_result.addUserFloat("scSeedSigmaIetaIeta", see);
      //  electron_result.addUserFloat("scSeedSigmaIphiIphi", spp); 


      // The one below is kept for historical reasons
      electron_result.addUserFloat("SC_seed_energy", el->superCluster()->seed()->energy());
      electron_result.addUserFloat("SC_seed_eta", el->superCluster()->seed()->eta());
    }
    else{
      electron_result.addUserFloat("SC_seed_energy", -1.);
      electron_result.addUserFloat("SC_seed_eta", -1.);
    }
    //
    //            //
    //            const BasicCluster&  clRef              = *(el->superCluster()->seed());
    //            const vector<float>& covs               = clusterTools_->covariances(clRef);                         // get the covariances computed in 5x5 around the seed
    //            const vector<float>& lcovs              = clusterTools_->localCovariances(clRef);                    // get the local covariances computed in a 5x5 around the seed
    //            const vector<float>  localCovariancesSC = clusterTools_->scLocalCovariances(*(el->superCluster()));  // get the local covariances computed using all crystals in the SC
    //
    //            //
    //get from RECO            electron_result.addUserFloat("sigmaIPhiIPhi", isfinite(lcovs[2]) ? lcovs[2] > 0 ? sqrt(lcovs[2]) : -1 * sqrt(-1 * lcovs[2]) : -1. );
    //
    //            //


    ////////
    // ID //
    ////////
    //electron_result.addUserFloat("hOverE", el->hadronicOverEm());
    electron_result.addUserFloat("hOverE", el->hcalOverEcal());
    electron_result.addUserFloat("full5x5_hOverE", el->full5x5_hcalOverEcal());
    electron_result.addUserFloat("eOverPIn", el->eSuperClusterOverP());
    electron_result.addUserFloat("eOverPOut", el->eEleClusterOverPout());
    electron_result.addUserFloat("fbrem", el->fbrem());
    electron_result.addUserFloat("dEtaIn", el->deltaEtaSuperClusterTrackAtVtx());
    electron_result.addUserFloat("dEtaOut", el->deltaEtaSeedClusterTrackAtCalo());
    electron_result.addUserFloat("dPhiIn", el->deltaPhiSuperClusterTrackAtVtx());
    electron_result.addUserFloat("dPhiOut", el->deltaPhiSeedClusterTrackAtCalo());

    ////////////
    // Charge //
    ////////////
    bool isGsfCtfScPixChargeConsistent = el->isGsfCtfScPixChargeConsistent();
    bool isGsfScPixChargeConsistent = el->isGsfScPixChargeConsistent();
    bool isGsfCtfChargeConsistent = el->isGsfCtfChargeConsistent();
    int trk_q = gsfTrack->charge();
    electron_result.addUserInt("charge", el->charge());
    //electron_result.addUserInt("threeCharge", el->threeCharge());
    electron_result.addUserInt("trk_charge", trk_q);
    electron_result.addUserInt("SC_pixCharge", el->scPixCharge());
    electron_result.addUserInt("isGsfCtfScPixChargeConsistent", static_cast<int>(isGsfCtfScPixChargeConsistent)); // Used also in id
    electron_result.addUserInt("isGsfScPixChargeConsistent", static_cast<int>(isGsfScPixChargeConsistent));
    electron_result.addUserInt("isGsfCtfChargeConsistent", static_cast<int>(isGsfCtfChargeConsistent));

    ////////////
    // Tracks //
    ////////////
    {
      float pt       = gsfTrack->pt();
      float p        = gsfTrack->p();
      float pz       = gsfTrack->pz();
      float pterr = (trk_q!=0 ? sqrt(pt*pt*p*p/pow(trk_q, 2)*(gsfTrack->covariance(0, 0))+2*pt*p/trk_q*pz*(gsfTrack->covariance(0, 1))+ pz*pz*(gsfTrack->covariance(1, 1))) : -1.);
      electron_result.addUserFloat("chi2", gsfTrack->chi2());
      electron_result.addUserFloat("ndof", gsfTrack->ndof());
      electron_result.addUserFloat("d0Err", gsfTrack->d0Error());
      electron_result.addUserFloat("z0Err", gsfTrack->dzError());
      electron_result.addUserFloat("ptErr", pterr);
      electron_result.addUserFloat("ptErr_gsf", gsfTrack->ptError());
      electron_result.addUserFloat("validHits", gsfTrack->numberOfValidHits());
      electron_result.addUserFloat("lostHits", gsfTrack->numberOfLostHits());
    }

    // Vertex dxy and dz
    if (firstGoodVertex!=vertexCollection->end()){
      electron_result.addUserFloat("dxy_PV", gsfTrack->dxy(firstGoodVertex->position()));
      electron_result.addUserFloat("dz_PV", gsfTrack->dz(firstGoodVertex->position()));
    }
    else{
      electron_result.addUserFloat("dxy_PV", -1.);
      electron_result.addUserFloat("dz_PV", -1.);
    }
    if (vertexCollection->begin()!=vertexCollection->end()){
      electron_result.addUserFloat("dz_firstPV", gsfTrack->dz(vertexCollection->begin()->position()));
      electron_result.addUserFloat("dxy_firstPV", gsfTrack->dxy(vertexCollection->begin()->position()));
    }
    else{
      electron_result.addUserFloat("dz_firstPV", -1.);
      electron_result.addUserFloat("dxy_firstPV", -1.);
    }

    /////////
    // CTF //
    /////////
    if (ctfTrack.isNonnull()){
      //electron_result.addUserFloat("trkshFrac", static_cast<float>( el->shFracInnerHits() )                                  );
      electron_result.addUserFloat("trkshFrac", static_cast<float>(el->ctfGsfOverlap()));
      electron_result.addUserFloat("trkdr", deltaR(gsfTrack->eta(), gsfTrack->phi(), ctfTrack->eta(), ctfTrack->phi()));
      electron_result.addUserFloat("ckf_chi2", ctfTrack->chi2());
      electron_result.addUserFloat("ckf_ndof", ctfTrack->ndof());
      electron_result.addUserInt("ckf_laywithmeas", ctfTrack->hitPattern().trackerLayersWithMeasurement());
      electron_result.addUserInt("ckf_charge", ctfTrack->charge());
    }
    else{
      electron_result.addUserFloat("trkshFrac", -1.);
      electron_result.addUserFloat("trkdr", -1.);
      electron_result.addUserFloat("ckf_chi2", -1.);
      electron_result.addUserFloat("ckf_ndof", -1.);
      electron_result.addUserInt("ckf_laywithmeas", -1.);
      electron_result.addUserInt("ckf_charge", 0);
    }


    //        ////////////////////
    //        // Regular Vertex //
    //        ////////////////////        
    //        TransientTrack tt = theTTBuilder->build(el->gsfTrack());
    //    
    //        if ( firstGoodVertex!=vertexCollection->end() ) {
    //            Measurement1D ip3D_regular = IPTools::absoluteImpactParameter3D(tt, *firstGoodVertex).second;
    //            //
    //            electron_result.addUserFloat("ip3d", ip3D_regular.value() );
    //            electron_result.addUserFloat("ip3derr", ip3D_regular.error() );
    //        } else {
    //            //
    //            electron_result.addUserFloat("ip3d", -999. );
    //            electron_result.addUserFloat("ip3derr", -1. );
    //        }


    //Impact Parameters
    //electron_result.addUserFloat("IP3D", el->ip3d()); // miniAOD
    electron_result.addUserFloat("IP3D", el->dB(pat::Electron::PV3D));
    electron_result.addUserFloat("IP3Derr", el->edB(pat::Electron::PV3D));
    electron_result.addUserFloat("IP2D", el->dB(pat::Electron::PV2D));
    electron_result.addUserFloat("IP2Derr", el->edB(pat::Electron::PV2D));

    /////////////////
    // Hit Pattern //
    /////////////////
    // Redesign according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingHitPatternRedesign
    const HitPattern& pattern = gsfTrack->hitPattern();
    //const HitPattern& p_inner = gsfTrack->trackerExpectedHitsInner(); 
    //const HitPattern& p_outer = gsfTrack->trackerExpectedHitsOuter();
    electron_result.addUserInt("n_missing_inner_hits", pattern.numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    electron_result.addUserInt("n_missing_outer_hits", pattern.numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
    electron_result.addUserInt("n_valid_pixel_hits", pattern.numberOfValidPixelHits());
    electron_result.addUserInt("n_lost_pixel_hits", pattern.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)); // Not sure about this. Could be MISSING_INNER_HITS instead.
    electron_result.addUserInt("n_tracker_layers", pattern.trackerLayersWithMeasurement());
    electron_result.addUserInt("n_tracker_layers_3D", pattern.pixelLayersWithMeasurement() + pattern.numberOfValidStripLayersWithMonoAndStereo());
    electron_result.addUserInt("n_tracker_layers_lost", pattern.trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));

    /////////////////
    // Conversions //
    /////////////////
    reco::ConversionRef conv_ref = ConversionTools::matchedConversion(*el, convs_h, beamSpot);
    float vertexFitProbability = -1;
    if (!conv_ref.isNull()){
      const reco::Vertex &vtx = conv_ref.get()->conversionVertex();
      if (vtx.isValid()) vertexFitProbability = TMath::Prob(vtx.chi2(), vtx.ndof());
    }
    electron_result.addUserFloat("conv_vtx_prob", vertexFitProbability);
    //cout<<"Found electron with pt eta phi "<<el->p4().pt() <<" "<< el->p4().eta() <<" "<< el->p4().phi()<<" and vertexFitProbability "<<vertexFitProbability<<endl;

    //////////////////////////////////////////////
    // Flag For Vertex Fit Conversion Rejection //
    //////////////////////////////////////////////
    electron_result.addUserInt("conv_vtx_flag", static_cast<int>(!el->passConversionVeto())); // PAT variable: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc#467

    //////////////////////
    // genMatch miniAOD //
    //////////////////////
    LorentzVector mc_p4(0, 0, 0, 0);
    const reco::GenParticle* gen = el->genParticle();
    if (gen){
      mc_p4 = gen->p4();
      electron_result.addUserInt("mc_patMatch_id", gen->pdgId());
      electron_result.addUserFloat("mc_patMatch_pt", mc_p4.pt());
      electron_result.addUserFloat("mc_patMatch_eta", mc_p4.eta());
      electron_result.addUserFloat("mc_patMatch_phi", mc_p4.phi());
      electron_result.addUserFloat("mc_patMatch_mass", mc_p4.mass());
      electron_result.addUserFloat("mc_patMatch_dr", ROOT::Math::VectorUtil::DeltaR(gen->p4(), el->p4()));
    }
    else{
      electron_result.addUserInt("mc_patMatch_id", 0);
      electron_result.addUserFloat("mc_patMatch_pt", -1);
      electron_result.addUserFloat("mc_patMatch_eta", 0);
      electron_result.addUserFloat("mc_patMatch_phi", 0);
      electron_result.addUserFloat("mc_patMatch_mass", 0);
      electron_result.addUserFloat("mc_patMatch_dr", -1);
    }

    //////////////////////
    // mini-isolation   //
    //////////////////////
    // float minichiso     = 0.;
    // float mininhiso     = 0.;
    // float miniemiso     = 0.;
    // float minidbiso     = 0.;
    // elMiniIso(el, true, 0.0, minichiso, mininhiso, miniemiso, minidbiso);
    // electron_result.addUserFloat("miniIso_uncor   ",  minichiso + mininhiso + miniemiso );
    // electron_result.addUserFloat("miniIso_ch      ",  minichiso );
    // electron_result.addUserFloat("miniIso_nh      ",  mininhiso );
    // electron_result.addUserFloat("miniIso_em      ",  miniemiso );
    // electron_result.addUserFloat("miniIso_db      ",  minidbiso );
    auto el2 = el->clone();
    auto miniiso = el2->miniPFIsolation();
    electron_result.addUserFloat("miniIso_uncor", miniiso.chargedHadronIso() + miniiso.neutralHadronIso() + miniiso.photonIso());
    electron_result.addUserFloat("miniIso_ch", miniiso.chargedHadronIso());
    electron_result.addUserFloat("miniIso_nh", miniiso.neutralHadronIso());
    electron_result.addUserFloat("miniIso_em", miniiso.photonIso());
    electron_result.addUserFloat("miniIso_db", miniiso.puChargedHadronIso());
    delete el2;

    ///////////////////////////
    // PFCluster isolation   //
    ///////////////////////////
    electron_result.addUserFloat("ecalPFClusterIso", el->ecalPFClusterIso());
    electron_result.addUserFloat("hcalPFClusterIso", el->hcalPFClusterIso());

    // Put the object into the result collection
    result->emplace_back(electron_result);
  }

  // Put the result collection into the event
  iEvent.put(std::move(result));
}

//----------------------------------------------------------------------------
// Electron Id classification function (a flag for the Sani type class)
//----------------------------------------------------------------------------
int ElectronMaker::classify(RefToBase<pat::Electron> const& electron) {
  double eOverP = electron->eSuperClusterOverP();
  double fbrem = electron->fbrem();

  int cat;
  if ((electron->isEB() && fbrem<0.06) || (electron->isEE() && fbrem<0.1)) cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) cat=0;
  else cat=2;
  return cat;
}

//little labour-saving function to get the reference to the ValueMap
template<typename T> const ValueMap<T>& ElectronMaker::getValueMap(Event const& iEvent, InputTag const& inputTag){
  Handle< ValueMap<T> > handle;
  iEvent.getByLabel(inputTag, handle);
  return *(handle.product());
}

double ElectronMaker::electronIsoValuePF(GsfElectron const& el, Vertex const& vtx, float coner, float minptn, float dzcut, float footprintdr, float gammastripveto, float elestripveto, int filterId){

  float pfciso = 0.;
  float pfniso = 0.;
  float pffootprint = 0.;
  float pfjurveto = 0.;
  float pfjurvetoq = 0.;

  //TrackRef siTrack     = el.closestCtfTrackRef();
  TrackRef siTrack     = el.closestTrack();
  GsfTrackRef gsfTrack = el.gsfTrack();

  if (gsfTrack.isNull() && siTrack.isNull()) return -1.;

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

void ElectronMaker::elIsoCustomCone(edm::View<pat::Electron>::const_iterator const& el, float dr, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso){
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

void ElectronMaker::elMiniIso(edm::View<pat::Electron>::const_iterator const& el, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float &dbiso){
  float pt = el->p4().pt();
  float dr = 0.2;
  if (pt>50) dr = 10./pt;
  if (pt>200) dr = 0.05;
  elIsoCustomCone(el, dr, useVetoCones, ptthresh, chiso, nhiso, emiso, dbiso);
  return;
}

void ElectronMaker::setMVAIdUserVariables(edm::View<pat::Electron>::const_iterator const& el, pat::Electron& electron_result, std::string const& id_name, std::string const& id_identifier) const{
  if (el->hasUserInt(id_name+"Categories")){
    electron_result.addUserFloat("id_MVA_"+id_identifier+"_Val", el->userFloat(id_name+"Values"));
    electron_result.addUserInt("id_MVA_"+id_identifier+"_Cat", el->userInt(id_name+"Categories"));
  }
  else throw cms::Exception("ElectronMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}
void ElectronMaker::setCutBasedIdUserVariables(edm::View<pat::Electron>::const_iterator const& el, pat::Electron& electron_result, std::string const& id_name, std::string const& id_identifier) const{
  if (el->isElectronIDAvailable(id_name)){
    electron_result.addUserInt("id_cutBased_"+id_identifier+"_Bits", el->userInt(id_name));
  }
  else throw cms::Exception("ElectronMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}


//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
