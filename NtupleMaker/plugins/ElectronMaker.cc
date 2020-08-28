#include <memory>
#include <cmath>
#include <limits>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

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
#include "TString.h"

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

#include "CMS3/NtupleMaker/interface/plugins/ElectronMaker.h"
#include "CMS3/NtupleMaker/interface/plugins/MatchUtilities.h"
#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/ElectronSelectionHelpers.h"
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include "CMS3/Dictionaries/interface/EgammaFiduciality.h"
#include "CMS3/Dictionaries/interface/ElectronTriggerCutEnums.h"

#include "CMSDataTools/AnalysisTree/interface/HelperFunctions.h"

#include "MELAStreamHelpers.hh"


using namespace std;
using namespace reco;
using namespace edm;
using namespace MELAStreamHelpers;


typedef math::XYZPoint Point;
typedef Ref<edmNew::DetSetVector<SiStripCluster>, SiStripCluster > ClusterRef;
typedef Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > pixel_ClusterRef;


ElectronMaker::ElectronMaker(const ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasPrefix")),
  year_(iConfig.getParameter<int>("year")),

  MVACuts_(iConfig.getParameter<edm::VParameterSet>("MVACuts"))
{
  setupMVACuts();

  electronsToken  = consumes<edm::View<pat::Electron>  >(iConfig.getParameter<edm::InputTag>("electronsInputTag"));

  beamSpotToken  = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotInputTag"));

  vtxToken  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));

  recoConversionToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("recoConversionInputTag"));

  ebhitsToken = consumes< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("EBHitsInputTag"));
  eehitsToken = consumes< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("EEHitsInputTag"));
  //eshitsToken = consumes< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("ESHitsInputTag"));

  rhoToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoInputTag"));
  rhoCaloToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoCaloInputTag"));

  produces<pat::ElectronCollection>().setBranchAlias(aliasprefix_);
}
ElectronMaker::~ElectronMaker(){}

void ElectronMaker::beginRun(const edm::Run&, const EventSetup& es){}
void ElectronMaker::beginJob(){}
void ElectronMaker::endJob(){}

void ElectronMaker::produce(Event& iEvent, const EventSetup& iSetup){
  auto result = std::make_unique<pat::ElectronCollection>();

  // uncertainties and corrections
  // somewhat complicated: see 
  // http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_1_2/doc/html/d5/d4b/GsfElectron_8h-source.html
  // note that if ecalEnergy == eSC depends on if further ecal corrections have been applied to the electron
  // after its construction

  // Event rho
  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting rho from the event...");
  const double rho_event = *rhoHandle;

  // Event rho (calo), used in HLT
  edm::Handle< double > rhoCaloHandle;
  iEvent.getByToken(rhoCaloToken, rhoCaloHandle);
  if (!rhoCaloHandle.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting calo rho from the event...");
  const double rho_calo_event = *rhoCaloHandle;

  ///////////////
  // Electrons //
  ///////////////
  Handle<View<pat::Electron> > els_h;
  iEvent.getByToken(electronsToken, els_h);
  if (!els_h.isValid()) throw cms::Exception("ElectronMaker::produce: error getting electron collection from Event!");

  ////////////
  // Vertex //
  ////////////
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vtxToken, vertexHandle);
  if (!vertexHandle.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting vertex collection from Event!");

  const VertexCollection* vertexCollection = vertexHandle.product();
  VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
  for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); vtx++){
    if (VertexSelectionHelpers::testGoodVertex(*vtx)){
      firstGoodVertex = vtx;
      break;
    }
  }
  bool hasVertex = (!vertexCollection->empty());
  bool hasGoodVertex = (firstGoodVertex!=vertexCollection->end());

  /////////////////
  // Conversions //
  /////////////////
  edm::Handle<reco::ConversionCollection> convs_h;
  iEvent.getByToken(recoConversionToken, convs_h);
  if (!convs_h.isValid()) throw cms::Exception("ElectronMaker::produce: error getting conversion collection");

  /////////////////////
  // Hit collections //
  /////////////////////
  edm::Handle< EcalRecHitCollection > ebhitsHandle;
  iEvent.getByToken(ebhitsToken, ebhitsHandle);
  if (!ebhitsHandle.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting EB hits from the event...");

  edm::Handle< EcalRecHitCollection > eehitsHandle;
  iEvent.getByToken(eehitsToken, eehitsHandle);
  if (!eehitsHandle.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting EE hits from the event...");

  /*
  edm::Handle< EcalRecHitCollection > eshitsHandle;
  iEvent.getByToken(eshitsToken, eshitsHandle);
  if (!eshitsHandle.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting ES hits from the event...");
  */

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
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken, beamSpotH);
  if (!beamSpotH.isValid()) throw cms::Exception("ElectronMaker::produce: Error getting the beam spot from the event...");
  const Point beamSpot = beamSpotH->position();

  /////////////////////////
  // Loop over electrons //
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

    //////////////////////
    // Fiduciality Mask //
    //////////////////////
    cms3_egamma_fid_type_mask_t fiducialityMask = 0;  // The enums are in interface/EgammaFiduciality.h
    if (el->isEB()) fiducialityMask |= 1 << ISEB;
    if (el->isEE()) fiducialityMask |= 1 << ISEE;
    if (el->isEBEEGap()) fiducialityMask |= 1 << ISEBEEGAP;
    if (el->isEBEtaGap()) fiducialityMask |= 1 << ISEBETAGAP;
    if (el->isEBPhiGap()) fiducialityMask |= 1 << ISEBPHIGAP;
    if (el->isEBGap()) fiducialityMask |= 1 << ISEBGAP;
    if (el->isEEDeeGap()) fiducialityMask |= 1 << ISEEDEEGAP;
    if (el->isEERingGap()) fiducialityMask |= 1 << ISEERINGGAP;
    if (el->isEEGap()) fiducialityMask |= 1 << ISEEGAP;
    if (el->isGap()) fiducialityMask |= 1 << ISGAP;
    electron_result.addUserInt("fid_mask", fiducialityMask);

    ///////////////////////////
    // Corrections & Seeding //
    ///////////////////////////
    cms3_egamma_fid_type_mask_t electronTypeMask = 0;
    if (el->isEcalEnergyCorrected()) electronTypeMask |= 1 << ISECALENERGYCORRECTED;
    if (el->trackerDrivenSeed()) electronTypeMask |= 1 << ISTRACKERDRIVEN;
    if (el->ecalDrivenSeed()) electronTypeMask |= 1 << ISECALDRIVEN;
    if (el->passingPflowPreselection()) electronTypeMask |= 1 << ISPFPRESELECTED;
    if (el->passingCutBasedPreselection()) electronTypeMask |= 1 << ISCUTPRESELECTED;
    // el->ecalDriven() == el->ecalDrivenSeed() && el->passingCutBasedPreselection(), so no need to keep a bit for that
    if (el->passingMvaPreselection()) electronTypeMask |= 1 << ISMVAPRESELECTED;
    //if ( el->isMomentumCorrected() ) electronTypeMask |= 1 << ISMOMENTUMCORRECTED; // Depricated in CMSSW_4_2x ( DataFormats/EgammaCandidates/interface/GsfElectron.h )
    electron_result.addUserInt("type_mask", electronTypeMask);

    ///////////////////
    // Predefined ID //
    ///////////////////
    // electron_result.addUserInt("category", classify(gsfElectron)); // this is the sani classification
    setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Fall17IsoV2", "Fall17V2_Iso");
    setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Fall17NoIsoV2", "Fall17V2_NoIso");
    switch (this->year_){
    case 2016:
      setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Summer16IdIso", "HZZRun2Legacy_Iso");
      break;
    case 2017:
      setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Fall17IsoV2", "HZZRun2Legacy_Iso");
      break;
    case 2018:
      setMVAIdUserVariables(el, electron_result, "ElectronMVAEstimatorRun2Autumn18IdIso", "HZZRun2Legacy_Iso");
      break;
    default:
      break;
    }
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-veto", "Fall17V2_Veto");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-loose", "Fall17V2_Loose");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-medium", "Fall17V2_Medium");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V2-tight", "Fall17V2_Tight");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-veto", "Fall17V1_Veto");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-loose", "Fall17V1_Loose");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-medium", "Fall17V1_Medium");
    setCutBasedIdUserVariables(el, electron_result, "cutBasedElectronID-Fall17-94X-V1-tight", "Fall17V1_Tight");
    applyTriggerEmulationCuts(rho_calo_event, el, electron_result, 0);
    applyTriggerEmulationCuts(rho_calo_event, el, electron_result, 1);

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

    // the p4 of the electron is the uncorrected one
    electron_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, el->eta(), el->phi(), uncorrected_mass));

    if (el->hasUserFloat("ecalTrkEnergyPostCorr")){
      electron_result.addUserFloat("scale_smear_corr", el->userFloat("ecalTrkEnergyPostCorr") / uncorrected_energy); // get scale/smear correction factor directly from miniAOD

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
    }
    else{ // Ensure that the user floats still exist
      electron_result.addUserFloat("scale_smear_corr", 1);

      electron_result.addUserFloat("scale_smear_corr_scale_totalUp", 1);
      electron_result.addUserFloat("scale_smear_corr_scale_statUp", 1);
      electron_result.addUserFloat("scale_smear_corr_scale_systUp", 1);
      electron_result.addUserFloat("scale_smear_corr_scale_gainUp", 1);

      electron_result.addUserFloat("scale_smear_corr_scale_totalDn", 1);
      electron_result.addUserFloat("scale_smear_corr_scale_statDn", 1);
      electron_result.addUserFloat("scale_smear_corr_scale_systDn", 1);
      electron_result.addUserFloat("scale_smear_corr_scale_gainDn", 1);

      electron_result.addUserFloat("scale_smear_corr_smear_totalUp", 1);
      electron_result.addUserFloat("scale_smear_corr_smear_statUp", 1);
      electron_result.addUserFloat("scale_smear_corr_smear_systUp", 1);
      electron_result.addUserFloat("scale_smear_corr_smear_gainUp", 1);

      electron_result.addUserFloat("scale_smear_corr_smear_totalDn", 1);
      electron_result.addUserFloat("scale_smear_corr_smear_statDn", 1);
      electron_result.addUserFloat("scale_smear_corr_smear_systDn", 1);
      electron_result.addUserFloat("scale_smear_corr_smear_gainDn", 1);
    }

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
    electron_result.addUserFloat("hcalDepth1TowerSumEt03", el->dr03HcalDepth1TowerSumEt());
    electron_result.addUserFloat("ecalIso03", el->dr03EcalRecHitSumEt());
    electron_result.addUserFloat("hcalIso03", el->dr03HcalTowerSumEt());
    electron_result.addUserFloat("tkIso03", el->dr03TkSumPt());

    electron_result.addUserFloat("ecalIso04", el->dr04EcalRecHitSumEt());
    electron_result.addUserFloat("hcalIso04", el->dr04HcalTowerSumEt());
    electron_result.addUserFloat("tkIso04", el->dr04TkSumPt());

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

    //////////////////
    // PF Isolation //
    //////////////////
    GsfElectron::PflowIsolationVariables const& pfIso = el->pfIsolationVariables();
    electron_result.addUserFloat("pfChargedHadronIso", pfIso.sumChargedHadronPt);
    electron_result.addUserFloat("pfNeutralHadronIso", pfIso.sumNeutralHadronEt);
    electron_result.addUserFloat("pfPhotonIso", pfIso.sumPhotonEt);
    electron_result.addUserFloat("pfPUIso", pfIso.sumPUPt);
    double pfIso03_sum_charged_nofsr=0, pfIso03_sum_neutral_nofsr=0;
    electron_result.addUserFloat("pfIso03_comb_nofsr", ElectronSelectionHelpers::electronPFIsoComb(*el, year_, ElectronSelectionHelpers::PFIso03, rho_event, 0., &pfIso03_sum_charged_nofsr, &pfIso03_sum_neutral_nofsr));
    electron_result.addUserFloat("pfIso03_sum_charged_nofsr", pfIso03_sum_charged_nofsr);
    electron_result.addUserFloat("pfIso03_sum_neutral_nofsr", pfIso03_sum_neutral_nofsr);
    double pfIso04_sum_charged_nofsr=0, pfIso04_sum_neutral_nofsr=0;
    electron_result.addUserFloat("pfIso04_comb_nofsr", ElectronSelectionHelpers::electronPFIsoComb(*el, year_, ElectronSelectionHelpers::PFIso04, rho_event, 0., &pfIso04_sum_charged_nofsr, &pfIso04_sum_neutral_nofsr));
    electron_result.addUserFloat("pfIso04_sum_charged_nofsr", pfIso04_sum_charged_nofsr);
    electron_result.addUserFloat("pfIso04_sum_neutral_nofsr", pfIso04_sum_neutral_nofsr);


    ////////
    // ID //
    ////////
    //electron_result.addUserFloat("full5x5_sigmaEtaEta", el->full5x5_sigmaEtaEta());
    electron_result.addUserFloat("full5x5_sigmaIEtaIEta", el->full5x5_sigmaIetaIeta());
    electron_result.addUserFloat("full5x5_sigmaIPhiIPhi", el->full5x5_sigmaIphiIphi());
    electron_result.addUserFloat("full5x5_r9", el->full5x5_r9());
    electron_result.addUserFloat("full5x5_hOverE", el->full5x5_hcalOverEcal());
    electron_result.addUserFloat("full5x5_e1x5", el->full5x5_e1x5());
    electron_result.addUserFloat("full5x5_e5x5", el->full5x5_e5x5());
    electron_result.addUserFloat("full5x5_e2x5Max", el->full5x5_e2x5Max());
    //electron_result.addUserFloat("sigmaEtaEta", el->sigmaEtaEta());
    electron_result.addUserFloat("sigmaIEtaIEta", el->sigmaIetaIeta());
    electron_result.addUserFloat("sigmaIPhiIPhi", el->sigmaIphiIphi());
    electron_result.addUserFloat("r9", el->r9());
    //electron_result.addUserFloat("hOverE", el->hadronicOverEm());
    electron_result.addUserFloat("hOverE", el->hcalOverEcal());

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
    unsigned char charge_consistency_bits=0;
    if (el->isGsfCtfChargeConsistent()) charge_consistency_bits |= (1 << 0);
    if (el->isGsfScPixChargeConsistent()) charge_consistency_bits |= (1 << 1);
    if (el->isGsfCtfScPixChargeConsistent()) charge_consistency_bits |= (1 << 2); // Used also in some ids
    int trk_q = gsfTrack->charge();
    electron_result.addUserInt("charge", el->charge());
    //electron_result.addUserInt("threeCharge", el->threeCharge());
    electron_result.addUserInt("trk_charge", trk_q);
    electron_result.addUserInt("SC_pixCharge", el->scPixCharge());
    electron_result.addUserInt("charge_consistency_bits", charge_consistency_bits);

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
      electron_result.addUserFloat("d0err", gsfTrack->d0Error());
      electron_result.addUserFloat("z0err", gsfTrack->dzError());
      electron_result.addUserFloat("pterr_calc", pterr);
      electron_result.addUserFloat("pterr_gsf", gsfTrack->ptError());
      electron_result.addUserInt("n_valid_hits", gsfTrack->numberOfValidHits());
      electron_result.addUserInt("n_lost_hits", gsfTrack->numberOfLostHits());
    }

    // Vertex dxy and dz
    if (hasGoodVertex){
      electron_result.addUserFloat("dxy_PV", gsfTrack->dxy(firstGoodVertex->position()));
      electron_result.addUserFloat("dz_PV", gsfTrack->dz(firstGoodVertex->position()));
    }
    else{
      electron_result.addUserFloat("dxy_PV", -999.);
      electron_result.addUserFloat("dz_PV", -999.);
    }
    if (hasVertex){
      electron_result.addUserFloat("dxy_firstPV", gsfTrack->dxy(vertexCollection->begin()->position()));
      electron_result.addUserFloat("dz_firstPV", gsfTrack->dz(vertexCollection->begin()->position()));
    }
    else{
      electron_result.addUserFloat("dxy_firstPV", -999.);
      electron_result.addUserFloat("dz_firstPV", -999.);
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
      electron_result.addUserInt("ckf_n_layers_with_measurement", ctfTrack->hitPattern().trackerLayersWithMeasurement());
      electron_result.addUserInt("ckf_charge", ctfTrack->charge());
    }
    else{
      electron_result.addUserFloat("trkshFrac", -1.);
      electron_result.addUserFloat("trkdr", -1.);
      electron_result.addUserFloat("ckf_chi2", -1.);
      electron_result.addUserFloat("ckf_ndof", -1.);
      electron_result.addUserInt("ckf_n_layers_with_measurement", -1.);
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
    double IP3D = el->dB(pat::Electron::PV3D);
    double IP3Derr = el->edB(pat::Electron::PV3D);
    double IP2D = el->dB(pat::Electron::PV2D);
    double IP2Derr = el->edB(pat::Electron::PV2D);
    electron_result.addUserFloat("IP3D", IP3D);
    electron_result.addUserFloat("IP3Derr", IP3Derr);
    electron_result.addUserFloat("SIP3D", (IP3Derr>0. ? IP3D / IP3Derr : 9999.));
    electron_result.addUserFloat("IP2D", IP2D);
    electron_result.addUserFloat("IP2Derr", IP2Derr);
    electron_result.addUserFloat("SIP2D", (IP2Derr>0. ? IP2D / IP2Derr : 9999.));

    /////////////////
    // Hit Pattern //
    /////////////////
    // Redesign according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingHitPatternRedesign
    const HitPattern& pattern = gsfTrack->hitPattern();
    //const HitPattern& p_inner = gsfTrack->trackerExpectedHitsInner(); 
    //const HitPattern& p_outer = gsfTrack->trackerExpectedHitsOuter();
    electron_result.addUserInt("n_missing_inner_hits", pattern.numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    electron_result.addUserInt("n_missing_outer_hits", pattern.numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
    electron_result.addUserInt("n_missing_pixel_hits", pattern.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS)); // Not sure about this. Could be MISSING_INNER_HITS instead.
    electron_result.addUserInt("n_all_missing_inner_hits", pattern.numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
    electron_result.addUserInt("n_all_missing_outer_hits", pattern.numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS));
    electron_result.addUserInt("n_valid_pixel_hits", pattern.numberOfValidPixelHits());
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
    // PAT variable: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc#467
    electron_result.addUserInt("conv_vtx_flag", static_cast<int>(!el->passConversionVeto()));

    // Compute the swiss cross variable and add seed time
    float E4overE1 = -99;
    float seedTime = 0;
    reco::CaloClusterPtr electron_seed = el->seed();
    if (electron_seed.isNonnull()){
      auto const& seedHitsAndFractions = electron_seed->hitsAndFractions();
      if (!seedHitsAndFractions.empty()){
        const DetId seedId = seedHitsAndFractions.front().first;
        float e1 = getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 0, 0, &seedTime);
        if (e1>0.f){
          float s4 = 0;
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 1, 0);
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), -1, 0);
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 0, 1);
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 0, -1);
          E4overE1 = s4/e1;
        }
      }
    }
    electron_result.addUserFloat("E4overE1", E4overE1);
    electron_result.addUserFloat("seedTime", seedTime);

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
    pat::PFIsolation const& miniiso = el->miniPFIsolation();
    electron_result.addUserFloat("miniIso_ch", miniiso.chargedHadronIso());
    electron_result.addUserFloat("miniIso_nh", miniiso.neutralHadronIso());
    electron_result.addUserFloat("miniIso_em", miniiso.photonIso());
    electron_result.addUserFloat("miniIso_db", miniiso.puChargedHadronIso());
    double miniIso_sum_charged_nofsr=0, miniIso_sum_neutral_nofsr=0;
    electron_result.addUserFloat("miniIso_comb_nofsr", ElectronSelectionHelpers::electronMiniIsoComb(*el, year_, rho_event, 0., &miniIso_sum_charged_nofsr, &miniIso_sum_neutral_nofsr));
    electron_result.addUserFloat("miniIso_comb_nofsr_uncorrected", miniiso.chargedHadronIso() + miniiso.neutralHadronIso() + miniiso.photonIso());
    electron_result.addUserFloat("miniIso_sum_charged_nofsr", miniIso_sum_charged_nofsr);
    electron_result.addUserFloat("miniIso_sum_neutral_nofsr", miniIso_sum_neutral_nofsr);

    /////////////////////////
    // PFCluster isolation //
    /////////////////////////
    electron_result.addUserFloat("ecalPFClusterIso", el->ecalPFClusterIso());
    electron_result.addUserFloat("hcalPFClusterIso", el->hcalPFClusterIso());

    ///////////////////////////
    // Associated candidates //
    ///////////////////////////
    auto associated_pfcands = el->associatedPackedPFCandidates();
    pat::PackedCandidate const* closestPFElectron_associated = nullptr;
    float min_dR_electron_pfelectron_associated = -1;
    std::vector<pat::PackedCandidate const*> pfelectroncands;
    unsigned int n_associated_pfcands = associated_pfcands.size();
    double associated_pfcands_sum_sc_pt = 0, associated_pfcands_sum_px = 0, associated_pfcands_sum_py = 0;
    for (auto const& pfcand:associated_pfcands){
      associated_pfcands_sum_sc_pt += pfcand->pt();
      associated_pfcands_sum_px += pfcand->px();
      associated_pfcands_sum_py += pfcand->py();
      if (std::abs(pfcand->pdgId()) == 11) pfelectroncands.push_back(&(*pfcand));
    }
    unsigned int n_associated_pfelectrons = pfelectroncands.size();
    // Do electron - PFcand matching
    {
      std::vector<pat::Electron const*> dummy_electron_list; dummy_electron_list.push_back(&(*el));
      std::unordered_map<pat::Electron const*, pat::PackedCandidate const*> patelectron_pfelectron_map;
      CMS3ObjectHelpers::matchParticles(
        CMS3ObjectHelpers::kMatchBy_DeltaR,
        dummy_electron_list.begin(), dummy_electron_list.end(),
        pfelectroncands.begin(), pfelectroncands.end(),
        patelectron_pfelectron_map
      );
      auto it_match = patelectron_pfelectron_map.find(&(*el));
      if (it_match != patelectron_pfelectron_map.end() && it_match->second){
        closestPFElectron_associated = it_match->second;
        min_dR_electron_pfelectron_associated = reco::deltaR(closestPFElectron_associated->p4(), el->p4());
      }
    }
    // Record
    electron_result.addUserInt("n_associated_pfcands", n_associated_pfcands);
    electron_result.addUserInt("n_associated_pfelectrons", n_associated_pfelectrons);
    electron_result.addUserFloat("associated_pfcands_sum_sc_pt", associated_pfcands_sum_sc_pt);
    electron_result.addUserFloat("associated_pfcands_sum_px", associated_pfcands_sum_px);
    electron_result.addUserFloat("associated_pfcands_sum_py", associated_pfcands_sum_py);
    electron_result.addUserFloat("min_dR_electron_pfelectron_associated", min_dR_electron_pfelectron_associated);
    // Only for checks
    /*
    {
      double dRmax_ch=-1;
      double dRmax_nh=-1;
      double dRmax_em=-1;
      double dRmax_mu=-1;
      double dRmax_e=-1;
      double dRmin_ch=-1;
      double dRmin_nh=-1;
      double dRmin_em=-1;
      double dRmin_mu=-1;
      double dRmin_e=-1;
      for (auto const& pfcand:associated_pfcands){
        //if (pfcand->vertexRef().key()!=0) continue;
        unsigned int abs_id = std::abs(pfcand->pdgId());
        double dr = reco::deltaR(el->p4(), pfcand->p4());
        if (abs_id == 211){
          dRmax_ch = std::max(dRmax_ch, dr);
          if (dRmin_ch<0. || dRmin_ch>dr) dRmin_ch = dr;
        }
        else if (abs_id == 130){
          dRmax_nh = std::max(dRmax_nh, dr);
          if (dRmin_nh<0. || dRmin_nh>dr) dRmin_nh = dr;
        }
        else if (abs_id == 22){
          dRmax_em = std::max(dRmax_em, dr);
          if (dRmin_em<0. || dRmin_em>dr) dRmin_em = dr;
        }
        else if (abs_id == 11){
          dRmax_e = std::max(dRmax_e, dr);
          if (dRmin_e<0. || dRmin_e>dr) dRmin_e = dr;
        }
        else if (abs_id == 13){
          dRmax_mu = std::max(dRmax_mu, dr);
          if (dRmin_mu<0. || dRmin_mu>dr) dRmin_mu = dr;
        }
        //if (pfcand->vertexRef().key()!=0){
        //  MELAout << "PF cand (id=" << pfcand->pdgId() << ", dR=" << dr << ", dz=" << pfcand->dzAssociatedPV() << ") has vtx " << pfcand->vertexRef().key() << endl;
        //}
      }
      electron_result.addUserFloat("dRmax_ch", dRmax_ch);
      electron_result.addUserFloat("dRmax_nh", dRmax_nh);
      electron_result.addUserFloat("dRmax_em", dRmax_em);
      electron_result.addUserFloat("dRmax_mu", dRmax_mu);
      electron_result.addUserFloat("dRmax_e", dRmax_e);
      electron_result.addUserFloat("dRmin_ch", dRmin_ch);
      electron_result.addUserFloat("dRmin_nh", dRmin_nh);
      electron_result.addUserFloat("dRmin_em", dRmin_em);
      electron_result.addUserFloat("dRmin_mu", dRmin_mu);
      electron_result.addUserFloat("dRmin_e", dRmin_e);
    }
    */

    // Not needed because momentum of electron and PF electron are within a 100 MeV of each other (much less at small pT)
    /*
    {
      float closestPFElectron_associated_px=0, closestPFElectron_associated_py=0;
      if (closestPFElectron_associated){
        closestPFElectron_associated_px = closestPFElectron_associated->px();
        closestPFElectron_associated_py = closestPFElectron_associated->py();
      }
      electron_result.addUserFloat("closestPFElectron_associated_px", closestPFElectron_associated_px);
      electron_result.addUserFloat("closestPFElectron_associated_py", closestPFElectron_associated_py);
    }
    */

    // Add EGamma PFPhoton ID
    // Use the associated PF photon candidate for MET safety checks
    setEGammaPFElectronIdSelectionBits(el, closestPFElectron_associated, electron_result);

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

void ElectronMaker::setupMVACuts(){
  for (edm::ParameterSet const& pset:MVACuts_){
    std::string wpLabel = pset.getParameter<std::string>("mvaLabel");
    {
      // Fix label in case 'wp' is not present
      std::vector<std::string> tmplist;
      HelperFunctions::splitOptionRecursive(wpLabel, tmplist, '_', false);
      if (tmplist.back().find("wp")==std::string::npos) tmplist.back().insert(0, "wp");
      for (size_t i=0; i<tmplist.size(); i++){
        if (i==0) wpLabel = tmplist.at(i);
        else wpLabel = wpLabel + "_" + tmplist.at(i);
      }
    }
    HelperFunctions::replaceString<std::string, const std::string>(wpLabel, "RawValues", "");
    HelperFunctions::replaceString<std::string, const std::string>(wpLabel, "Values", "");

    //std::cout << "ElectronMaker::setupMVACuts: Adding MVA WPs for " << wpLabel << std::endl;

    std::string mvaLabel = wpLabel;
    {
      std::vector<std::string> tmplist;
      HelperFunctions::splitOptionRecursive(mvaLabel, tmplist, '_', false);
      mvaLabel="";
      for (size_t i=0; i<tmplist.size()-1; i++){
        if (i==0) mvaLabel = tmplist.at(i);
        else mvaLabel = mvaLabel + "_" + tmplist.at(i);
      }
    }
    std::vector<std::string> wpCuts = pset.getParameter<std::vector<std::string>>("mvaCuts");
    std::vector< StringCutObjectSelector<pat::Electron, true> > cutobjs; cutobjs.reserve(wpCuts.size());
    for (std::string strcut:wpCuts){
      strcut = strcut + " )";
      strcut.insert(0, "Values') > ( "); strcut.insert(0, mvaLabel); strcut.insert(0, "userFloat('");
      cutobjs.emplace_back(strcut);
    }
    MVACutObjects[wpLabel] = cutobjs;
    //std::cout << "ElectronMaker::setupMVACuts: Adding " << wpLabel << " cuts on " << mvaLabel << ":" << std::endl;
    //for (auto const& c:wpCuts) std::cout << "\t" << c << std::endl;
  }
}
void ElectronMaker::setMVAIdUserVariables(edm::View<pat::Electron>::const_iterator const& el, pat::Electron& electron_result, std::string const& id_name, std::string const& id_identifier) const{
  if (el->hasUserInt(id_name+"Categories")){
    // Set MVA category and value
    float val = el->userFloat(id_name+"Values");
    int cat = el->userInt(id_name+"Categories");
    if (cat<0 || static_cast<unsigned int>(cat)>std::numeric_limits<cms3_electron_mvacat_t>::max()){
      throw cms::Exception(Form("ElectronMaker::setMVAIdUserVariables: Id %s has an out-of-bounds category label %i.", id_name.data(), cat));
    }

    std::string strRecord = "id_MVA_"+id_identifier;
    electron_result.addUserFloat(strRecord+"_Val", val);
    electron_result.addUserInt(strRecord+"_Cat", cat);

    // Embed MVA cuts
    for (auto const& it_cuts:MVACutObjects){
      std::string const& label = it_cuts.first;
      std::string strpass = label;
      HelperFunctions::replaceString<std::string, const std::string>(strpass, id_name, (strRecord+"_pass"));
      if (label.find(id_name)!=std::string::npos && cat>=0){
        auto const& cutobj = it_cuts.second.at(cat);
        electron_result.addUserInt(strpass, static_cast<int>(cutobj(*el)));
      }
    }
  }
  else throw cms::Exception("ElectronMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}
void ElectronMaker::setCutBasedIdUserVariables(edm::View<pat::Electron>::const_iterator const& el, pat::Electron& electron_result, std::string const& id_name, std::string const& id_identifier) const{
  if (el->isElectronIDAvailable(id_name)){
    unsigned int cutbits = el->userInt(id_name);
    if (cutbits>std::numeric_limits<cms3_electron_cutbasedbits_t>::max()){
      throw cms::Exception(Form("ElectronMaker::setCutBasedIdUserVariables: Id %s has an out-of-bounds cut bits value %u.", id_name.data(), cutbits));
    }

    electron_result.addUserInt("id_cutBased_"+id_identifier+"_Bits", cutbits);
  }
  else throw cms::Exception("ElectronMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}

void ElectronMaker::setEGammaPFElectronIdSelectionBits(edm::View<pat::Electron>::const_iterator const& electron, pat::PackedCandidate const* pfCand, pat::Electron& electron_result) const{
  // Selection flow follows RecoParticleFlow/PFProducer/src/PFEGammaFilters.cc
  cms3_electron_cutbasedbits_egPFElectron_t id_egamma_pfElectron_Bits = 0;

  if (ElectronSelectionHelpers::testEGammaPFElectronSelection(*electron, year_)) id_egamma_pfElectron_Bits |= 1 << ISEGAMMAPFELECTRON_BASE;
  if (ElectronSelectionHelpers::testEGammaPFElectronPrimarySelection(*electron, year_)) id_egamma_pfElectron_Bits |= 1 << ISEGAMMAPFELECTRON_PRIMARY;
  if (ElectronSelectionHelpers::testEGammaPFElectronMETSafetySelection(*electron, pfCand, year_)) id_egamma_pfElectron_Bits |= 1 << ISEGAMMAPFELECTRON_METSAFE;

  electron_result.addUserInt("id_egamma_pfElectron_Bits", id_egamma_pfElectron_Bits);
}

void ElectronMaker::applyTriggerEmulationCuts(double const& rho, edm::View<pat::Electron>::const_iterator const& el, pat::Electron& electron_result, unsigned int const& trigVersion) const{
  const GsfTrackRef gsfTrack = el->gsfTrack();
  const HitPattern& pattern = gsfTrack->hitPattern();
  //int n_all_missing_inner_hits = pattern.numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
  int n_missing_inner_hits = pattern.numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  double const gsftrk_chisq = gsfTrack->chi2();
  double const ecalEnergy = el->correctedEcalEnergy(); // Same as GsfElectron::ecalEnergy()
  double const etaSC = el->superCluster()->eta();
  double const abs_etaSC = std::abs(etaSC);
  double const pt_uncorrected = el->p4().pt();
  double const& energy_uncorrected = ecalEnergy;
  // ID variables
  double const sigmaIEtaIEta_full5x5 = el->full5x5_sigmaIetaIeta();
  //double const sigmaIEtaIEta = el->sigmaIetaIeta();
  double const hOverE = el->hcalOverEcal(); // Same as GsfElectron::hadronicOverEm()
  //double const hOverE_full5x5 = el->full5x5_hcalOverEcal();
  double const eOverPIn = el->eSuperClusterOverP();
  //double const eOverPOut = el->eEleClusterOverPout();
  //double const eSeedOverPIn = el->eSeedClusterOverP();
  //double const fbrem = el->fbrem();
  double const dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
  double const abs_dEtaIn = std::abs(dEtaIn);
  //double const dEtaInSeed = el->deltaEtaSeedClusterTrackAtVtx();
  //double const abs_dEtaInSeed = std::abs(dEtaInSeed);
  double const dEtaOut = el->deltaEtaSeedClusterTrackAtCalo();
  double const abs_dEtaOut = std::abs(dEtaOut);
  double const dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
  double const abs_dPhiIn = std::abs(dPhiIn);
  //double const dPhiOut = el->deltaPhiSeedClusterTrackAtCalo();
  //double const abs_dPhiOut = std::abs(dPhiOut);
  double const Einv_minus_Pinv_in = std::abs(1./ecalEnergy - eOverPIn/ecalEnergy);
  // Isolation variables
  double const trkIso03 = el->dr03TkSumPt();
  double const ecalPFClusterIso = el->ecalPFClusterIso();
  double const hcalPFClusterIso = el->hcalPFClusterIso();
  double const ea_ecalPFClusterIso = ElectronSelectionHelpers::electronEffArea_ECALPFCluster(*el, this->year_, trigVersion);
  double const ea_hcalPFClusterIso = ElectronSelectionHelpers::electronEffArea_HCALPFCluster(*el, this->year_, trigVersion);
  double const ecalPFClusterIsoCorr = std::max(0., ecalPFClusterIso - rho*ea_ecalPFClusterIso);
  double const hcalPFClusterIsoCorr = std::max(0., hcalPFClusterIso - rho*ea_hcalPFClusterIso);
  //double const ecalIso04 = el->dr04EcalRecHitSumEt();
  //double const hcalIso04 = el->dr04HcalTowerSumEt();

  using namespace ElectronTriggerCutEnums;

  unsigned int trigcutbits = 0;
  for (unsigned int itype=0; itype<(int) nElectronTriggerCutTypes; itype++){
    bool res = true;

    ElectronTriggerCutTypes type = (ElectronTriggerCutTypes) itype;
    if (type == kCuts_CustomTriggerSafe_Id_v1){
      if (abs_etaSC < 1.479){
        if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
        if (abs_dEtaIn > 0.01) res = false;
        if (abs_dPhiIn > 0.04) res = false;
        if (Einv_minus_Pinv_in > 0.01) res = false;
      }
      else if (abs_etaSC < 2.5){
        if (sigmaIEtaIEta_full5x5 > 0.031) res = false;
        if (abs_dEtaIn > 0.01) res = false;
        if (abs_dPhiIn > 0.08) res = false;
        if (Einv_minus_Pinv_in > 0.01) res = false;
      }
      if (hOverE > 0.08) res = false;
    }
    else if (type == kCuts_CustomTriggerSafe_Iso_v1){
      if (ecalPFClusterIso/pt_uncorrected > 0.45) res = false;
      if (hcalPFClusterIso/pt_uncorrected > 0.25) res = false;
      if (trkIso03/pt_uncorrected > 0.2) res = false;
    }
    else if (type == kCuts_CustomTriggerSafe_Id_v2){
      if (abs_etaSC < 1.479){
        if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
        if (abs_dEtaOut > 0.004) res = false;
        if (abs_dPhiIn > 0.020) res = false;
        if (Einv_minus_Pinv_in > 0.013) res = false;
      }
      else if (abs_etaSC < 2.5){
        if (sigmaIEtaIEta_full5x5 > 0.031) res = false;
        if (Einv_minus_Pinv_in > 0.013) res = false;
        if (gsftrk_chisq > 3) res = false;
      }
      if (hOverE > 0.06) res = false;
    }
    else if (type == kCuts_CustomTriggerSafe_Iso_v2){
      if (abs_etaSC < 1.479){
        if (ecalPFClusterIsoCorr/pt_uncorrected > 0.16) res = false;
        if (hcalPFClusterIsoCorr/pt_uncorrected > 0.12) res = false;
        if (trkIso03/pt_uncorrected > 0.080) res = false;
      }
      else if (abs_etaSC < 2.5){
        if (ecalPFClusterIsoCorr/pt_uncorrected > 0.12) res = false;
        if (hcalPFClusterIsoCorr/pt_uncorrected > 0.12) res = false;
        if (trkIso03/pt_uncorrected > 0.08) res = false;
      }
    }
    else{
      if (this->year_ == 2016){
        switch (type){
        case kCuts_CaloIdL_TrackIdL:
        {
          // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*
          // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*
          // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.015) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.1) res = false;
          }
          if (hOverE > 0.13) res = false;
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v1:
        {
          // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.015) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.1) res = false;
          }
          if (hOverE > 0.13) res = false;
          if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5) res = false;
          if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3) res = false;
          if (trkIso03/pt_uncorrected > 0.2) res = false;
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v2:
        {
          // HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
          }
          if (hOverE > 0.13) res = false;
          if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5) res = false;
          if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3) res = false;
          if (trkIso03/pt_uncorrected > 0.2) res = false;
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v3:
        {
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.014) res = false;
            if (hOverE > 0.15) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (hOverE > 0.1) res = false;
          }
          if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5) res = false;
          if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3) res = false;
          if (trkIso03/pt_uncorrected > 0.2) res = false;
          // There is also PM S2 requirement, but we ignore that.
          break;
        }

        case kCuts_CaloIdL_GsfTrkIdVL:
        {
          // HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.014) res = false;
            if (hOverE > 0.15) res = false;
            if (abs_dEtaOut > 0.02) res = false;
            if (abs_dPhiIn > 0.15) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (hOverE > 0.1) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.02) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.15) res = false;
          }
          break;
        }

        case kCuts_CaloIdL_MW:
        {
          // HLT_DoubleEle33_CaloIdL_MW_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.014) res = false;
            if (hOverE > 0.15) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (hOverE > 0.1) res = false;
          }
          break;
        }

        case kCuts_WPLoose:
        {
          // HLT_Ele27_eta2p1_WPLoose_Gsf_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
            if (hOverE > 0.04) res = false;
            if (abs_dEtaOut > (trigVersion == 0 ? 0.004 : 0.0044)) res = false;
            if (abs_dPhiIn > (trigVersion == 0 ? 0.02 : 0.025)) res = false;
            if (Einv_minus_Pinv_in > 0.012) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.15 : 0.04)) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.15 : 0.06)) res = false;
            if (trkIso03/pt_uncorrected > (trigVersion == 0 ? 0.08 : 0.085)) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > (trigVersion == 0 ? 0.032 : 0.031)) res = false;
            if (hOverE > (trigVersion == 0 ? 0.065 : 0.04)) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.006) res = false;
            if (Einv_minus_Pinv_in > 0.01) res = false;
            if (gsftrk_chisq > (trigVersion == 0 ? 2.8 : 4.0)) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.15 : 0.045)) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.16 : 0.08)) res = false;
            if (trkIso03/pt_uncorrected > (trigVersion == 0 ? 0.08 : 0.07)) res = false;
            if (n_missing_inner_hits > 2) res = false;
          }
          break;
        }

        case kCuts_WPTight_v1:
        {
          // HLT_Ele25_eta2p1_WPTight_Gsf_v*
          // HLT_Ele27_WPTight_Gsf_v*
          // HLT_Ele32_eta2p1_WPTight_Gsf_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 >(trigVersion == 0 ? 0.011 : 0.0105)) res = false;
            if (hOverE > (trigVersion == 0 ? 0.05 : 0.07)) res = false;
            if (abs_dEtaOut > (trigVersion == 0 ? 0.0035 : 0.004)) res = false;
            if (abs_dPhiIn > (trigVersion == 0 ? 0.015 : 0.02)) res = false;
            if (Einv_minus_Pinv_in > 0.01) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.145 : 0.032)) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.13 : 0.055)) res = false;
            if (trkIso03/pt_uncorrected > (trigVersion == 0 ? 0.06 : 0.07)) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > (trigVersion == 0 ? 0.029 : 0.0285)) res = false;
            if (hOverE > (trigVersion == 0 ? 0.05 : 0.035)) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.0045) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.022) res = false;
            if (Einv_minus_Pinv_in > 0.008) res = false;
            if (gsftrk_chisq > (trigVersion == 0 ? 2.5 : 3.5)) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.105 : 0.04)) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.095 : 0.05)) res = false;
            if (trkIso03/pt_uncorrected > (trigVersion == 0 ? 0.06 : 0.05)) res = false;
            if (n_missing_inner_hits > 1) res = false;
          }
          break;
        }

        case kCuts_DoublePhoton:
        {
          // HLT_DoublePhoton60_v*
          if (abs_etaSC < 1.479){
            if (hOverE > 0.15) res = false;
          }
          else{
            if (hOverE > 0.1) res = false;
          }
          break;
        }

        case kCuts_WPTight_v2:
        {
          res = false;
          break;
        }
        default:
          throw cms::Exception("UndefinedElectronTrigger") << "ElectronMaker::applyTriggerCuts: Trigger enum" << itype << " is not implemented for year " << this->year_ << endl;
          break;
        }
      }
      else if (this->year_ == 2017){
        switch (type){
        case kCuts_CaloIdL_TrackIdL:
        {
          // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v* (version 1 only)
          // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v* (version 1 only)
          // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v* (version 0 or 1)
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.015) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.1) res = false;
          }
          if (hOverE > 0.13) res = false;
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v1:
        {
          // HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.5 : 0.5 + 0.29*rho)) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.3 : 0.3 + 0.2*rho)) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.015) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.1) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.5 : 0.5 + 0.21*rho)) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > (trigVersion == 0 ? 0.3 : 0.3 + 0.25*rho)) res = false;
          }
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v2:
        {
          // HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5 + 0.29*rho) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3 + 0.2*rho) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (abs_dEtaOut > 0.015) res = false;
            if (abs_dPhiIn > 0.1) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5 + 0.21*rho) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3 + 0.25*rho) res = false;
          }
          break;
        }

        case kCuts_CaloIdL_MW:
        {
          // HLT_DoubleEle33_CaloIdL_MW_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.014) res = false;
            if (hOverE > 0.15) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (hOverE > 0.1) res = false;
          }
          // There is also PM S2<250 requirement, but we ignore that.
          break;
        }


        // HLT_Ele35_WPTight_Gsf_v*: Progression is kCuts_WPTight_v2 V1, kCuts_WPTight_v1 V2, kCuts_WPTight_v2 V2
        // HLT_Ele38_WPTight_Gsf_v*: Progression is kCuts_WPTight_v1 V1, kCuts_WPTight_v1 V2, kCuts_WPTight_v2 V2
        // HLT_Ele40_WPTight_Gsf_v*: Progression is kCuts_WPTight_v1 V1, kCuts_WPTight_v1 V2, kCuts_WPTight_v2 V2
        case kCuts_WPTight_v1:
        {
          if (trigVersion == 0){
            if (abs_etaSC < 1.479){
              if (sigmaIEtaIEta_full5x5 > 0.0105) res = false;
              if (hOverE > 0.07) res = false;
              if (abs_dEtaOut > 0.004) res = false;
              if (abs_dPhiIn > 0.02) res = false;
              if (Einv_minus_Pinv_in > 0.01) res = false;
              if (ecalPFClusterIsoCorr/pt_uncorrected > 0.032) res = false;
              if (hcalPFClusterIsoCorr/pt_uncorrected > 0.055) res = false;
              if (trkIso03/pt_uncorrected > 0.07) res = false;
            }
            else{
              if (sigmaIEtaIEta_full5x5 > 0.0285) res = false;
              if (hOverE > 0.035) res = false;
              if (Einv_minus_Pinv_in > 0.008) res = false;
              if (gsftrk_chisq > 3.5) res = false;
              if (ecalPFClusterIsoCorr/pt_uncorrected > 0.04) res = false;
              if (hcalPFClusterIsoCorr/pt_uncorrected > 0.05) res = false;
              if (trkIso03/pt_uncorrected > 0.05) res = false;
              if (n_missing_inner_hits > 1) res = false;
            }
          }
          else{
            if (abs_etaSC < 1.){
              if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
              if (hOverE > (0.887 + 0.02*energy_uncorrected + 0.098*rho)/energy_uncorrected) res = false;
              if (abs_dEtaOut > 0.004) res = false;
              if (abs_dPhiIn > 0.02) res = false;
              if (Einv_minus_Pinv_in > 0.012) res = false;
              if (ecalPFClusterIsoCorr > (-0.581 + 0.03*pt_uncorrected + 0.325*rho)) res = false;
              if (hcalPFClusterIsoCorr > (0.786 + 0.03*pt_uncorrected + 0.259*rho)) res = false;
              if (trkIso03 > (0.838 + 0.03*pt_uncorrected + 0.029*rho)) res = false;
            }
            else if (abs_etaSC < 1.479){
              if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
              if (hOverE > (1.476 + 0.02*energy_uncorrected + 0.159*rho)/energy_uncorrected) res = false;
              if (abs_dEtaOut > 0.004) res = false;
              if (abs_dPhiIn > 0.02) res = false;
              if (Einv_minus_Pinv_in > 0.012) res = false;
              if (ecalPFClusterIsoCorr > (-0.698 + 0.03*pt_uncorrected + 0.296*rho)) res = false;
              if (hcalPFClusterIsoCorr > (0.298 + 0.03*pt_uncorrected + 0.328*rho)) res = false;
              if (trkIso03 > (-0.385 + 0.03*pt_uncorrected + 0.111*rho)) res = false;
            }
            else if (abs_etaSC < 2.1){
              if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
              if (hOverE > (2.672 + 0.015*energy_uncorrected + 0.353*rho)/energy_uncorrected) res = false;
              if (Einv_minus_Pinv_in > 0.011) res = false;
              if (ecalPFClusterIsoCorr > (-0.892 + 0.025*pt_uncorrected + 0.283*rho)) res = false;
              if (hcalPFClusterIsoCorr > (0.402 + 0.025*pt_uncorrected + 0.414*rho)) res = false;
              if (trkIso03 > (-0.363 + 0.025*pt_uncorrected + 0.114*rho)) res = false;
              if (n_missing_inner_hits > 1) res = false;
            }
            else{
              if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
              if (hOverE > (5.095 + 0.015*energy_uncorrected + 0.423*rho)/energy_uncorrected) res = false;
              if (Einv_minus_Pinv_in > 0.011) res = false;
              if (gsftrk_chisq > 2.5) res = false;
              if (ecalPFClusterIsoCorr > (-0.885 + 0.025*pt_uncorrected + 0.438*rho)) res = false;
              if (hcalPFClusterIsoCorr > (-0.061 + 0.025*pt_uncorrected + 0.456*rho)) res = false;
              if (trkIso03 > (0.702 + 0.025*pt_uncorrected + 0.032*rho)) res = false;
              if (n_missing_inner_hits > 1) res = false;
            }
          }
          break;
        }
        case kCuts_WPTight_v2:
        {
          if (trigVersion == 0){
            if (abs_etaSC < 1.479){
              if (sigmaIEtaIEta_full5x5 > 0.0105) res = false;
              if (hOverE > 0.07) res = false;
              if (abs_dEtaOut > 0.004) res = false;
              if (abs_dPhiIn > 0.02) res = false;
              if (Einv_minus_Pinv_in > 0.01) res = false;
              if (ecalPFClusterIsoCorr/pt_uncorrected > 0.032) res = false;
              if (hcalPFClusterIsoCorr/pt_uncorrected > 0.055) res = false;
              if (trkIso03/pt_uncorrected > 0.07) res = false;
            }
            else{
              if (sigmaIEtaIEta_full5x5 > 0.0285) res = false;
              if (hOverE > 0.035) res = false;
              if (abs_dEtaOut > 0.0045) res = false;
              if (abs_dPhiIn > 0.022) res = false;
              if (Einv_minus_Pinv_in > 0.008) res = false;
              if (gsftrk_chisq > 3.5) res = false;
              if (ecalPFClusterIsoCorr/pt_uncorrected > 0.04) res = false;
              if (hcalPFClusterIsoCorr/pt_uncorrected > 0.05) res = false;
              if (trkIso03/pt_uncorrected > 0.05) res = false;
              if (n_missing_inner_hits > 1) res = false;
            }
          }
          else{
            if (abs_etaSC < 1.){
              if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
              if (hOverE > (0.887 + 0.02*energy_uncorrected + 0.098*rho)/energy_uncorrected) res = false;
              if (abs_dEtaOut > 0.004) res = false;
              if (abs_dPhiIn > 0.02) res = false;
              if (Einv_minus_Pinv_in > 0.012) res = false;
              if (ecalPFClusterIsoCorr > (-0.581 + 0.03*pt_uncorrected + 0.325*rho)) res = false;
              if (hcalPFClusterIsoCorr > (0.786 + 0.03*pt_uncorrected + 0.259*rho)) res = false;
              if (trkIso03 > (0.838 + 0.03*pt_uncorrected + 0.029*rho)) res = false;
            }
            else if (abs_etaSC < 1.479){
              if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
              if (hOverE > (1.476 + 0.02*energy_uncorrected + 0.159*rho)/energy_uncorrected) res = false;
              if (abs_dEtaOut > 0.004) res = false;
              if (abs_dPhiIn > 0.02) res = false;
              if (Einv_minus_Pinv_in > 0.012) res = false;
              if (ecalPFClusterIsoCorr > (-0.698 + 0.03*pt_uncorrected + 0.296*rho)) res = false;
              if (hcalPFClusterIsoCorr > (0.298 + 0.03*pt_uncorrected + 0.328*rho)) res = false;
              if (trkIso03 > (-0.385 + 0.03*pt_uncorrected + 0.111*rho)) res = false;
            }
            else if (abs_etaSC < 2.1){
              if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
              if (hOverE > (2.672 + 0.015*energy_uncorrected + 0.353*rho)/energy_uncorrected) res = false;
              if (abs_dEtaOut > 0.005) res = false;
              if (abs_dPhiIn > 0.023) res = false;
              if (Einv_minus_Pinv_in > 0.011) res = false;
              if (ecalPFClusterIsoCorr > (-0.892 + 0.025*pt_uncorrected + 0.283*rho)) res = false;
              if (hcalPFClusterIsoCorr > (0.402 + 0.025*pt_uncorrected + 0.414*rho)) res = false;
              if (trkIso03 > (-0.363 + 0.025*pt_uncorrected + 0.114*rho)) res = false;
              if (n_missing_inner_hits > 1) res = false;
            }
            else{
              if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
              if (hOverE > (5.095 + 0.015*energy_uncorrected + 0.423*rho)/energy_uncorrected) res = false;
              if (abs_dEtaOut > 0.005) res = false;
              if (abs_dPhiIn > 0.023) res = false;
              if (Einv_minus_Pinv_in > 0.011) res = false;
              if (gsftrk_chisq > 2.5) res = false;
              if (ecalPFClusterIsoCorr > (-0.885 + 0.025*pt_uncorrected + 0.438*rho)) res = false;
              if (hcalPFClusterIsoCorr > (-0.061 + 0.025*pt_uncorrected + 0.456*rho)) res = false;
              if (trkIso03 > (0.702 + 0.025*pt_uncorrected + 0.032*rho)) res = false;
              if (n_missing_inner_hits > 1) res = false;
            }
          }
          break;
        }

        case kCuts_DoublePhoton:
        {
          // HLT_DoublePhoton70_v*
          if (abs_etaSC < 1.479){
            if (hOverE > 0.15) res = false;
          }
          else{
            if (hOverE > 0.1) res = false;
          }
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v3:
        case kCuts_CaloIdL_GsfTrkIdVL:
        case kCuts_WPLoose:
        {
          res = false;
          break;
        }

        default:
          throw cms::Exception("UndefinedElectronTrigger") << "ElectronMaker::applyTriggerCuts: Trigger enum" << itype << " is not implemented for year " << this->year_ << endl;
          break;
        }

      }
      else if (this->year_ == 2018){
        switch (type){
        case kCuts_CaloIdL_TrackIdL:
        {
          // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*
          // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*
          // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.015) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.1) res = false;
          }
          if (hOverE > 0.13) res = false;
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v1:
        {
          // HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*
          // HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*
          // HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.013) res = false;
            if (abs_dEtaOut > 0.01) res = false;
            if (abs_dPhiIn > 0.07) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5 + 0.29*rho) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3 + 0.2*rho) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.015) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.1) res = false;
            if (ecalPFClusterIsoCorr/pt_uncorrected > 0.5 + 0.21*rho) res = false;
            if (hcalPFClusterIsoCorr/pt_uncorrected > 0.3 + 0.25*rho) res = false;
          }
          if (hOverE > 0.13) res = false;
          if (trkIso03/pt_uncorrected > 0.2) res = false;
          break;
        }

        case kCuts_CaloIdL_MW:
        {
          // HLT_DoubleEle25_CaloIdL_MW_v*
          if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.014) res = false;
            if (hOverE > 0.15) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.035) res = false;
            if (hOverE > 0.1) res = false;
          }
          // There is also PM S2<250 requirement, but we ignore that.
          break;
        }

        case kCuts_WPTight_v1:
        {
          // HLT_Ele32_WPTight_Gsf_v*
          if (abs_etaSC < 1.){
            if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
            if (hOverE > (0.75 + 0.03*energy_uncorrected + 0.1*rho)/energy_uncorrected) res = false;
            if (abs_dEtaOut > 0.004) res = false;
            if (abs_dPhiIn > 0.02) res = false;
            if (Einv_minus_Pinv_in > 0.012) res = false;
            if (ecalPFClusterIsoCorr > (1.75 + 0.03*pt_uncorrected + 0.2*rho)) res = false;
            if (hcalPFClusterIsoCorr > (2.5 + 0.03*pt_uncorrected + 0.2*rho)) res = false;
            if (trkIso03 > (0.838 + 0.03*pt_uncorrected + 0.029*rho)) res = false;
          }
          else if (abs_etaSC < 1.479){
            if (sigmaIEtaIEta_full5x5 > 0.011) res = false;
            if (hOverE > (2.25 + 0.03*energy_uncorrected + 0.1*rho)/energy_uncorrected) res = false;
            if (abs_dEtaOut > 0.004) res = false;
            if (abs_dPhiIn > 0.02) res = false;
            if (Einv_minus_Pinv_in > 0.012) res = false;
            if (ecalPFClusterIsoCorr > (1.75 + 0.03*pt_uncorrected + 0.2*rho)) res = false;
            if (hcalPFClusterIsoCorr > (3.0 + 0.03*pt_uncorrected + 0.2*rho)) res = false;
            if (trkIso03 > (-0.385 + 0.03*pt_uncorrected + 0.111*rho)) res = false;
          }
          else if (abs_etaSC < 2.0){
            if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
            if (hOverE > (3. + 0.03*energy_uncorrected + 0.3*rho)/energy_uncorrected) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.005) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.023) res = false;
            if (Einv_minus_Pinv_in > 0.011) res = false;
            if (ecalPFClusterIsoCorr > (1. + 0.025*pt_uncorrected + 0.25*rho)) res = false;
            if (hcalPFClusterIsoCorr > (1. + 0.03*pt_uncorrected + 0.4*rho)) res = false;
            if (trkIso03 > (-0.363 + 0.025*pt_uncorrected + 0.114*rho)) res = false;
            if (n_missing_inner_hits > 1) res = false;
          }
          else if (abs_etaSC < 2.1){
            if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
            if (hOverE > (3. + 0.03*energy_uncorrected + 0.3*rho)/energy_uncorrected) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.005) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.023) res = false;
            if (Einv_minus_Pinv_in > 0.011) res = false;
            if (ecalPFClusterIsoCorr > (1. + 0.025*pt_uncorrected + 0.25*rho)) res = false;
            if (hcalPFClusterIsoCorr > (1. + 0.03*pt_uncorrected + 0.5*rho)) res = false;
            if (trkIso03 > (-0.363 + 0.025*pt_uncorrected + 0.114*rho)) res = false;
            if (n_missing_inner_hits > 1) res = false;
          }
          else{
            if (sigmaIEtaIEta_full5x5 > 0.0305) res = false;
            if (hOverE > (5. + 0.03*energy_uncorrected + 0.5*rho)/energy_uncorrected) res = false;
            if (trigVersion == 1 && abs_dEtaOut > 0.005) res = false;
            if (trigVersion == 1 && abs_dPhiIn > 0.023) res = false;
            if (Einv_minus_Pinv_in > 0.011) res = false;
            if (gsftrk_chisq > 2.5) res = false;
            if (ecalPFClusterIsoCorr > (2. + 0.025*pt_uncorrected + 0.3*rho)) res = false;
            if (hcalPFClusterIsoCorr > (2. + 0.03*pt_uncorrected + 0.5*rho)) res = false;
            if (trkIso03 > (0.702 + 0.025*pt_uncorrected + 0.032*rho)) res = false;
            if (n_missing_inner_hits > 1) res = false;
          }
          break;
        }

        case kCuts_DoublePhoton:
        {
          // HLT_DoublePhoton70_v*
          if (abs_etaSC < 1.479){
            if (hOverE > 0.15) res = false;
          }
          else{
            if (hOverE > 0.1) res = false;
          }
          break;
        }

        case kCuts_CaloIdL_TrackIdL_IsoVL_v2:
        case kCuts_CaloIdL_TrackIdL_IsoVL_v3:
        case kCuts_CaloIdL_GsfTrkIdVL:
        case kCuts_WPLoose:
        case kCuts_WPTight_v2:
        {
          res = false;
          break;
        }

        default:
          throw cms::Exception("UndefinedElectronTrigger") << "ElectronMaker::applyTriggerCuts: Trigger enum" << itype << " is not implemented for year " << this->year_ << endl;
          break;
        }

      }
      else throw cms::Exception("UndefinedElectronTrigger") << "ElectronMaker::applyTriggerEmulationCuts: Year " << this->year_ << " is not implemented." << endl;
    }

    HelperFunctions::set_bit(trigcutbits, itype, res);
  }

  electron_result.addUserInt(Form("id_cutBased_triggerEmulationV%u_Bits", trigVersion+1), trigcutbits);
}

float ElectronMaker::getRecHitEnergyTime(DetId const& id, EcalRecHitCollection const* ebhits, EcalRecHitCollection const* eehits, unsigned short di, unsigned short dj, float* outtime){
  if (!ebhits || !eehits) return 0;

  DetId nid = DetId(0);
  bool isEB = false;
  if (id.subdetId() == EcalBarrel){
    nid = EBDetId::offsetBy(id, di, dj);
    isEB = true;
  }
  else if (id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy(id, di, dj);

  if (nid == DetId(0)) return 0;

  EcalRecHitCollection const& rechits = (isEB ? *ebhits : *eehits);
  for (auto const& rechit:rechits){
    if (rechit.detid() == nid){
      if (outtime) *outtime = rechit.time();
      return rechit.energy();
    }
  }
  return 0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
