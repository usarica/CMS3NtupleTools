#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <CMS3/NtupleMaker/interface/plugins/CMS3Ntuplizer.h>
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace edm;


CMS3Ntuplizer::CMS3Ntuplizer(const edm::ParameterSet& pset_) :
  pset(pset_),
  outtree(nullptr),
  firstEvent(true),

  year(pset.getParameter<int>("year")),
  treename(pset.getUntrackedParameter<std::string>("treename")),
  isMC(pset.getParameter<bool>("isMC"))
{
  if (year!=2016 && year!=2017 && year!=2018) throw cms::Exception("CMS3Ntuplizer::CMS3Ntuplizer: Year is undefined!");

  electronsToken  = consumes< edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("electronSrc"));
  photonsToken  = consumes< edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("photonSrc"));
  //muonsToken  = consumes< edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonSrc"));

  genInfoToken = consumes<GenInfo>(pset.getParameter<edm::InputTag>("genInfoSrc"));

}
CMS3Ntuplizer::~CMS3Ntuplizer(){
  //clearMELABranches(); // Cleans LHE branches
  //delete pileUpReweight;
  //delete metCorrHandler;
}


void CMS3Ntuplizer::beginJob(){
  edm::Service<TFileService> fs;
  TTree* tout = fs->make<TTree>(treename, "Selected event summary");
  outtree = new BaseTree(nullptr, tout, nullptr, nullptr, false);
  outtree->setAcquireTreePossession(false);
  //buildMELABranches();
}
void CMS3Ntuplizer::endJob(){
  delete outtree;
}

void CMS3Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}
void CMS3Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

void CMS3Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&){}
void CMS3Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&){}

void CMS3Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  /*****************/
  /* Input handles */
  /*****************/

  // Electrons
  edm::Handle< edm::View<pat::Electron> > electronsHandle;
  iEvent.getByToken(electronsToken, electronsHandle);
  if (!electronsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::analyze: Error getting the electron collection from the event...");

  // Photons
  edm::Handle< edm::View<pat::Photon> > photonsHandle;
  iEvent.getByToken(photonsToken, photonsHandle);
  if (!photonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::analyze: Error getting the photon collection from the event...");

  // Muons
  edm::Handle< edm::View<pat::Muon> > muonsHandle;
  iEvent.getByToken(muonsToken, muonsHandle);
  if (!muonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::analyze: Error getting the muon collection from the event...");


  /********************************/
  /* Set the communicator entries */
  /********************************/
  /*
  When naeing variables, try to be conscious of the nanoAOD naming conventions, but do not make a big fuss about them either!
  The latest list of variables are documented at https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html
  */

  // Convenience macros to easily make and push vector values
#define MAKE_VECTOR_WITH_RESERVE(type_, name_, size_) std::vector<type_> name_; name_.reserve(size_);
#define PUSH_VECTOR_WITH_NAME(name_, var_) commonEntry.setNamedVal(TString(name_)+"_"+#var_, var_);

  // Electrons
  {
    //const char colName[] = "Electron"; // nanoAOD
    const char colName[] = "electrons";

    size_t n_electrons = electronsHandle->size();

    MAKE_VECTOR_WITH_RESERVE(float, pt, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, eta, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, phi, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, mass, n_electrons);

    MAKE_VECTOR_WITH_RESERVE(int, charge, n_electrons);

    // Has no convention correspondence in nanoAOD
    MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalUp, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalDn, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalUp, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalDn, n_electrons);

    MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_Iso_Val, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_Iso_Cat, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_NoIso_Val, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_NoIso_Cat, n_electrons);

    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Veto_Bits, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Loose_Bits, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Medium_Bits, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Tight_Bits, n_electrons);

    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Veto_Bits, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Loose_Bits, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Medium_Bits, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Tight_Bits, n_electrons);

    MAKE_VECTOR_WITH_RESERVE(unsigned int, fid_mask, n_electrons);
    MAKE_VECTOR_WITH_RESERVE(unsigned int, type_mask, n_electrons);

    for (View<pat::Electron>::const_iterator obj = electronsHandle->begin(); obj != electronsHandle->end(); obj++){
      // Core particle quantities
      // Uncorrected p4
      pt.push_back(obj->pt());
      eta.push_back(obj->eta());
      phi.push_back(obj->phi());
      mass.push_back(obj->mass());

      // Charge: Can obtain pdgId from this, so no need to record pdgId again
      charge.push_back(obj->userInt("charge"));

      // Scale and smear
      // Nominal value: Needs to multiply the uncorrected p4 at analysis level
      scale_smear_corr.push_back(obj->userFloat("scale_smear_corr"));
      // Uncertainties: Only store total up/dn for the moment
      scale_smear_corr_scale_totalUp.push_back(obj->userFloat("scale_smear_corr_scale_totalUp"));
      scale_smear_corr_scale_totalDn.push_back(obj->userFloat("scale_smear_corr_scale_totalDn"));
      scale_smear_corr_smear_totalUp.push_back(obj->userFloat("scale_smear_corr_smear_totalUp"));
      scale_smear_corr_smear_totalDn.push_back(obj->userFloat("scale_smear_corr_smear_totalDn"));

      // Id variables
      // Fall17V2_Iso MVA id
      id_MVA_Fall17V2_Iso_Val.push_back(obj->userFloat("id_MVA_Fall17V2_Iso_Val"));
      id_MVA_Fall17V2_Iso_Cat.push_back(obj->userInt("id_MVA_Fall17V2_Iso_Cat"));

      // Fall17V2_NoIso MVA id
      id_MVA_Fall17V2_NoIso_Val.push_back(obj->userFloat("id_MVA_Fall17V2_NoIso_Val"));
      id_MVA_Fall17V2_NoIso_Cat.push_back(obj->userInt("id_MVA_Fall17V2_NoIso_Cat"));

      // Fall17V2 cut-based ids
      id_cutBased_Fall17V2_Veto_Bits.push_back(obj->userInt("id_cutBased_Fall17V2_Veto_Bits"));
      id_cutBased_Fall17V2_Loose_Bits.push_back(obj->userInt("id_cutBased_Fall17V2_Loose_Bits"));
      id_cutBased_Fall17V2_Medium_Bits.push_back(obj->userInt("id_cutBased_Fall17V2_Medium_Bits"));
      id_cutBased_Fall17V2_Tight_Bits.push_back(obj->userInt("id_cutBased_Fall17V2_Tight_Bits"));

      // Fall17V1 cut-based ids
      id_cutBased_Fall17V1_Veto_Bits.push_back(obj->userInt("id_cutBased_Fall17V1_Veto_Bits"));
      id_cutBased_Fall17V1_Loose_Bits.push_back(obj->userInt("id_cutBased_Fall17V1_Loose_Bits"));
      id_cutBased_Fall17V1_Medium_Bits.push_back(obj->userInt("id_cutBased_Fall17V1_Medium_Bits"));
      id_cutBased_Fall17V1_Tight_Bits.push_back(obj->userInt("id_cutBased_Fall17V1_Tight_Bits"));

      // Masks
      fid_mask.push_back(obj->userInt("fid_mask"));
      type_mask.push_back(obj->userInt("type_mask"));
    }

    // Pass collections to the communicator
    PUSH_VECTOR_WITH_NAME(colName, pt);
    PUSH_VECTOR_WITH_NAME(colName, eta);
    PUSH_VECTOR_WITH_NAME(colName, phi);
    PUSH_VECTOR_WITH_NAME(colName, mass);

    PUSH_VECTOR_WITH_NAME(colName, charge);

    // Has no convention correspondence in nanoAOD
    PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr);
    PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalUp);
    PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalDn);
    PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalUp);
    PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalDn);

    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Val);
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Cat);
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Val);
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Cat);

    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Veto_Bits);
    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);

    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Veto_Bits);
    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Loose_Bits);
    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Medium_Bits);
    PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Tight_Bits);

    PUSH_VECTOR_WITH_NAME(colName, fid_mask);
    PUSH_VECTOR_WITH_NAME(colName, type_mask);
  }

  // Muons
  {
    //const char colName[] = "Muon"; // nanoAOD
    const char colName[] = "muons";

    size_t n_muons = muonsHandle->size();

    MAKE_VECTOR_WITH_RESERVE(float, pt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, eta, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, phi, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, mass, n_muons);

    MAKE_VECTOR_WITH_RESERVE(int, charge, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, type, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, caloCompatibility, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, segmentCompatibility, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, numberOfMatchedStations, n_muons);

    MAKE_VECTOR_WITH_RESERVE(int, selectors, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, simType, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, simExtType, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, gfit_chi2, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, gfit_ndof, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, gfit_validSTAHits, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, gfit_pt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, gfit_eta, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, gfit_phi, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, gfit_algo, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, gfit_ptErr, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, bfit_pt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, bfit_eta, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, bfit_phi, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, bfit_algo, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, bfit_ptErr, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, trkKink, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, chi2LocalPosition, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, chi2LocalMomentum, n_muons);

    MAKE_VECTOR_WITH_RESERVE(int, pid_TMLastStationLoose, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pid_TMLastStationTight, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pid_TM2DCompatibilityLoose, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pid_TM2DCompatibilityTight, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pid_TMOneStationTight, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pid_PFMuon, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, ecal_time, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, hcal_time, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, iso_trckvetoDep, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, iso_ecalvetoDep, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, iso_hcalvetoDep, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, iso_hovetoDep, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, iso03_sumPt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, iso03_emEt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, iso03_hadEt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, iso03_ntrk, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, trk_pt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, trk_eta, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, trk_phi, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, validHits, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, lostHits, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, d0Err, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, z0Err, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, ptErr, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, algo, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, algoOrig, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, nlayers, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, validPixelHits, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, exp_innerlayers, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, exp_outerlayers, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, dxyPV, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, dzPV, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, dz_firstPV, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, dxy_firstPV, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_ChargedHadronPt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_ChargedParticlePt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_NeutralHadronEt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_PhotonEt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_sumNeutralHadronEtHighThreshold, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_sumPhotonEtHighThreshold, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR03_pf_PUPt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_ChargedHadronPt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_ChargedParticlePt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_NeutralHadronEt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_PhotonEt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_sumNeutralHadronEtHighThreshold, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_sumPhotonEtHighThreshold, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, isoR04_pf_PUPt, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, pfpt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, pfeta, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, pfphi, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, pfmass, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pfcharge, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pfparticleId, n_muons);
    MAKE_VECTOR_WITH_RESERVE(int, pfidx, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, ip3d, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, ip3derr, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, ip2d, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, ip2derr, n_muons);

    MAKE_VECTOR_WITH_RESERVE(int, mc_patMatch_id, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, mc_patMatch_pt, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, mc_patMatch_eta, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, mc_patMatch_phi, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, mc_patMatch_mass, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, mc_patMatch_dr, n_muons);

    MAKE_VECTOR_WITH_RESERVE(float, miniIso_uncor, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, miniIso_ch, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, miniIso_nh, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, miniIso_em, n_muons);
    MAKE_VECTOR_WITH_RESERVE(float, miniIso_db, n_muons);

    for (View<pat::Muon>::const_iterator obj = muonsHandle->begin(); obj != muonsHandle->end(); obj++){
      // Core particle quantities
      pt.push_back(obj->pt());
      eta.push_back(obj->eta());
      phi.push_back(obj->phi());
      mass.push_back(obj->mass());
      charge.push_back(obj->userInt("charge"));
      type.push_back(obj->userInt("type"));
      caloCompatibility.push_back(obj->userFloat("caloCompatibility"));
      segmentCompatibility.push_back(obj->userFloat("segmentCompatibility"));
      numberOfMatchedStations.push_back(obj->userInt("numberOfMatchedStations"));

      // Selectors
      selectors.push_back(obj->userInt("selectors"));
      simType.push_back(obj->userInt("simType"));
      simExtType.push_back(obj->userInt("simExtType"));

      // Global fit track
      gfit_chi2.push_back(obj->userFloat("gfit_chi2"));
      gfit_ndof.push_back(obj->userInt("gfit_ndof"));
      gfit_validSTAHits.push_back(obj->userInt("gfit_validSTAHits"));
      gfit_pt.push_back(obj->userFloat("gfit_pt"));
      gfit_eta.push_back(obj->userFloat("gfit_eta"));
      gfit_phi.push_back(obj->userFloat("gfit_phi"));
      gfit_algo.push_back(obj->userInt("gfit_algo"));
      gfit_ptErr.push_back(obj->userFloat("gfit_ptErr"));

      // Best fit track
      bfit_pt.push_back(obj->userFloat("bfit_pt"));
      bfit_eta.push_back(obj->userFloat("bfit_eta"));
      bfit_phi.push_back(obj->userFloat("bfit_phi"));
      bfit_algo.push_back(obj->userInt("bfit_algo"));
      bfit_ptErr.push_back(obj->userFloat("bfit_ptErr"));

      // Muon Quality
      trkKink.push_back(obj->userFloat("trkKink"));
      chi2LocalPosition.push_back(obj->userFloat("chi2LocalPosition"));
      chi2LocalMomentum.push_back(obj->userFloat("chi2LocalMomentum"));

      // ID
      pid_TMLastStationLoose.push_back(obj->userInt("pid_TMLastStationLoose"));
      pid_TMLastStationTight.push_back(obj->userInt("pid_TMLastStationTight"));
      pid_TM2DCompatibilityLoose.push_back(obj->userInt("pid_TM2DCompatibilityLoose"));
      pid_TM2DCompatibilityTight.push_back(obj->userInt("pid_TM2DCompatibilityTight"));
      pid_TMOneStationTight.push_back(obj->userInt("pid_TMOneStationTight"));
      pid_PFMuon.push_back(obj->userInt("pid_PFMuon"));

      // Energy
      ecal_time.push_back(obj->userFloat("ecal_time"));
      hcal_time.push_back(obj->userFloat("hcal_time"));

      // Isolation
      iso_trckvetoDep.push_back(obj->userFloat("iso_trckvetoDep"));
      iso_ecalvetoDep.push_back(obj->userFloat("iso_ecalvetoDep"));
      iso_hcalvetoDep.push_back(obj->userFloat("iso_hcalvetoDep"));
      iso_hovetoDep.push_back(obj->userFloat("iso_hovetoDep"));
      iso03_sumPt.push_back(obj->userFloat("iso03_sumPt"));
      iso03_emEt.push_back(obj->userFloat("iso03_emEt"));
      iso03_hadEt.push_back(obj->userFloat("iso03_hadEt"));
      iso03_ntrk.push_back(obj->userInt("iso03_ntrk"));

      // Tracks
      trk_pt.push_back(obj->userFloat("trk_pt"));
      trk_eta.push_back(obj->userFloat("trk_eta"));
      trk_phi.push_back(obj->userFloat("trk_phi"));
      validHits.push_back(obj->userInt("validHits"));
      lostHits.push_back(obj->userInt("lostHits"));
      d0Err.push_back(obj->userFloat("d0Err"));
      z0Err.push_back(obj->userFloat("z0Err"));
      ptErr.push_back(obj->userFloat("ptErr"));
      algo.push_back(obj->userInt("algo"));
      algoOrig.push_back(obj->userInt("algoOrig"));
      nlayers.push_back(obj->userInt("nlayers"));
      validPixelHits.push_back(obj->userInt("validPixelHits"));
      exp_innerlayers.push_back(obj->userInt("exp_innerlayers"));
      exp_outerlayers.push_back(obj->userInt("exp_outerlayers"));
      dxyPV.push_back(obj->userFloat("dxyPV"));
      dzPV.push_back(obj->userFloat("dzPV"));
      dz_firstPV.push_back(obj->userFloat("dz_firstPV"));
      dxy_firstPV.push_back(obj->userFloat("dxy_firstPV"));

      // PF isolation
      isoR03_pf_ChargedHadronPt.push_back(obj->userFloat("isoR03_pf_ChargedHadronPt"));
      isoR03_pf_ChargedParticlePt.push_back(obj->userFloat("isoR03_pf_ChargedParticlePt"));
      isoR03_pf_NeutralHadronEt.push_back(obj->userFloat("isoR03_pf_NeutralHadronEt"));
      isoR03_pf_PhotonEt.push_back(obj->userFloat("isoR03_pf_PhotonEt"));
      isoR03_pf_sumNeutralHadronEtHighThreshold.push_back(obj->userFloat("isoR03_pf_sumNeutralHadronEtHighThreshold"));
      isoR03_pf_sumPhotonEtHighThreshold.push_back(obj->userFloat("isoR03_pf_sumPhotonEtHighThreshold"));
      isoR03_pf_PUPt.push_back(obj->userFloat("isoR03_pf_PUPt"));

      isoR04_pf_ChargedHadronPt.push_back(obj->userFloat("isoR04_pf_ChargedHadronPt"));
      isoR04_pf_ChargedParticlePt.push_back(obj->userFloat("isoR04_pf_ChargedParticlePt"));
      isoR04_pf_NeutralHadronEt.push_back(obj->userFloat("isoR04_pf_NeutralHadronEt"));
      isoR04_pf_PhotonEt.push_back(obj->userFloat("isoR04_pf_PhotonEt"));
      isoR04_pf_sumNeutralHadronEtHighThreshold.push_back(obj->userFloat("isoR04_pf_sumNeutralHadronEtHighThreshold"));
      isoR04_pf_sumPhotonEtHighThreshold.push_back(obj->userFloat("isoR04_pf_sumPhotonEtHighThreshold"));
      isoR04_pf_PUPt.push_back(obj->userFloat("isoR04_pf_PUPt"));

      // PF particle
      pfpt.push_back(obj->userFloat("pfpt"));
      pfeta.push_back(obj->userFloat("pfeta"));
      pfphi.push_back(obj->userFloat("pfphi"));
      pfmass.push_back(obj->userFloat("pfmass"));
      pfcharge.push_back(obj->userInt("pfcharge"));
      pfparticleId.push_back(obj->userInt("pfparticleId"));
      pfidx.push_back(obj->userInt("pfidx"));

      // IP 3D
      ip3d.push_back(obj->userFloat("ip3d"));
      ip3derr.push_back(obj->userFloat("ip3derr"));
      ip2d.push_back(obj->userFloat("ip2d"));
      ip2derr.push_back(obj->userFloat("ip2derr"));

      // genMatch miniAOD
      mc_patMatch_id.push_back(obj->userInt("mc_patMatch_id"));
      mc_patMatch_pt.push_back(obj->userFloat("mc_patMatch_pt"));
      mc_patMatch_eta.push_back(obj->userFloat("mc_patMatch_eta"));
      mc_patMatch_phi.push_back(obj->userFloat("mc_patMatch_phi"));
      mc_patMatch_mass.push_back(obj->userFloat("mc_patMatch_mass"));
      mc_patMatch_dr.push_back(obj->userFloat("mc_patMatch_dr"));

      // mini-isolation
      miniIso_uncor.push_back(obj->userFloat("miniIso_uncor"));
      miniIso_ch.push_back(obj->userFloat("miniIso_ch"));
      miniIso_nh.push_back(obj->userFloat("miniIso_nh"));
      miniIso_em.push_back(obj->userFloat("miniIso_em"));
      miniIso_db.push_back(obj->userFloat("miniIso_db"));

    }

    // Pass collections to the communicator
    PUSH_VECTOR_WITH_NAME(colName, pt);
    PUSH_VECTOR_WITH_NAME(colName, eta);
    PUSH_VECTOR_WITH_NAME(colName, phi);
    PUSH_VECTOR_WITH_NAME(colName, mass);

    PUSH_VECTOR_WITH_NAME(colName, charge);
    PUSH_VECTOR_WITH_NAME(colName, type);

    PUSH_VECTOR_WITH_NAME(colName, caloCompatibility);
    PUSH_VECTOR_WITH_NAME(colName, segmentCompatibility);
    PUSH_VECTOR_WITH_NAME(colName, numberOfMatchedStations);

    PUSH_VECTOR_WITH_NAME(colName, selectors);
    PUSH_VECTOR_WITH_NAME(colName, simType);
    PUSH_VECTOR_WITH_NAME(colName, simExtType);

    PUSH_VECTOR_WITH_NAME(colName, gfit_chi2);
    PUSH_VECTOR_WITH_NAME(colName, gfit_ndof);
    PUSH_VECTOR_WITH_NAME(colName, gfit_validSTAHits);
    PUSH_VECTOR_WITH_NAME(colName, gfit_pt);
    PUSH_VECTOR_WITH_NAME(colName, gfit_eta);
    PUSH_VECTOR_WITH_NAME(colName, gfit_phi);
    PUSH_VECTOR_WITH_NAME(colName, gfit_algo);
    PUSH_VECTOR_WITH_NAME(colName, gfit_ptErr);

    PUSH_VECTOR_WITH_NAME(colName, bfit_pt);
    PUSH_VECTOR_WITH_NAME(colName, bfit_eta);
    PUSH_VECTOR_WITH_NAME(colName, bfit_phi);
    PUSH_VECTOR_WITH_NAME(colName, bfit_algo);
    PUSH_VECTOR_WITH_NAME(colName, bfit_ptErr);

    PUSH_VECTOR_WITH_NAME(colName, trkKink);
    PUSH_VECTOR_WITH_NAME(colName, chi2LocalPosition);
    PUSH_VECTOR_WITH_NAME(colName, chi2LocalMomentum);

    PUSH_VECTOR_WITH_NAME(colName, pid_TMLastStationLoose);
    PUSH_VECTOR_WITH_NAME(colName, pid_TMLastStationTight);
    PUSH_VECTOR_WITH_NAME(colName, pid_TM2DCompatibilityLoose);
    PUSH_VECTOR_WITH_NAME(colName, pid_TM2DCompatibilityTight);
    PUSH_VECTOR_WITH_NAME(colName, pid_TMOneStationTight);
    PUSH_VECTOR_WITH_NAME(colName, pid_PFMuon);

    PUSH_VECTOR_WITH_NAME(colName, ecal_time);
    PUSH_VECTOR_WITH_NAME(colName, hcal_time);

    PUSH_VECTOR_WITH_NAME(colName, iso_trckvetoDep);
    PUSH_VECTOR_WITH_NAME(colName, iso_ecalvetoDep);
    PUSH_VECTOR_WITH_NAME(colName, iso_hcalvetoDep);
    PUSH_VECTOR_WITH_NAME(colName, iso_hovetoDep);
    PUSH_VECTOR_WITH_NAME(colName, iso03_sumPt);
    PUSH_VECTOR_WITH_NAME(colName, iso03_emEt);
    PUSH_VECTOR_WITH_NAME(colName, iso03_hadEt);
    PUSH_VECTOR_WITH_NAME(colName, iso03_ntrk);

    PUSH_VECTOR_WITH_NAME(colName, trk_pt);
    PUSH_VECTOR_WITH_NAME(colName, trk_eta);
    PUSH_VECTOR_WITH_NAME(colName, trk_phi);
    PUSH_VECTOR_WITH_NAME(colName, validHits);
    PUSH_VECTOR_WITH_NAME(colName, lostHits);
    PUSH_VECTOR_WITH_NAME(colName, d0Err);
    PUSH_VECTOR_WITH_NAME(colName, z0Err);
    PUSH_VECTOR_WITH_NAME(colName, ptErr);
    PUSH_VECTOR_WITH_NAME(colName, algo);
    PUSH_VECTOR_WITH_NAME(colName, algoOrig);
    PUSH_VECTOR_WITH_NAME(colName, nlayers);
    PUSH_VECTOR_WITH_NAME(colName, validPixelHits);
    PUSH_VECTOR_WITH_NAME(colName, exp_innerlayers);
    PUSH_VECTOR_WITH_NAME(colName, exp_outerlayers);
    PUSH_VECTOR_WITH_NAME(colName, dxyPV);
    PUSH_VECTOR_WITH_NAME(colName, dzPV);
    PUSH_VECTOR_WITH_NAME(colName, dz_firstPV);
    PUSH_VECTOR_WITH_NAME(colName, dxy_firstPV);

    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_ChargedHadronPt);
    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_ChargedParticlePt);
    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_NeutralHadronEt);
    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_PhotonEt);
    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_sumNeutralHadronEtHighThreshold);
    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_sumPhotonEtHighThreshold);
    PUSH_VECTOR_WITH_NAME(colName, isoR03_pf_PUPt);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_ChargedHadronPt);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_ChargedParticlePt);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_NeutralHadronEt);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_PhotonEt);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_sumNeutralHadronEtHighThreshold);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_sumPhotonEtHighThreshold);
    PUSH_VECTOR_WITH_NAME(colName, isoR04_pf_PUPt);

    PUSH_VECTOR_WITH_NAME(colName, pfpt);
    PUSH_VECTOR_WITH_NAME(colName, pfeta);
    PUSH_VECTOR_WITH_NAME(colName, pfphi);
    PUSH_VECTOR_WITH_NAME(colName, pfmass);
    PUSH_VECTOR_WITH_NAME(colName, pfcharge);
    PUSH_VECTOR_WITH_NAME(colName, pfparticleId);
    PUSH_VECTOR_WITH_NAME(colName, pfidx);
    PUSH_VECTOR_WITH_NAME(colName, ip3d);
    PUSH_VECTOR_WITH_NAME(colName, ip3derr);
    PUSH_VECTOR_WITH_NAME(colName, ip2d);
    PUSH_VECTOR_WITH_NAME(colName, ip2derr);

    PUSH_VECTOR_WITH_NAME(colName, mc_patMatch_id);
    PUSH_VECTOR_WITH_NAME(colName, mc_patMatch_pt);
    PUSH_VECTOR_WITH_NAME(colName, mc_patMatch_eta);
    PUSH_VECTOR_WITH_NAME(colName, mc_patMatch_phi);
    PUSH_VECTOR_WITH_NAME(colName, mc_patMatch_mass);
    PUSH_VECTOR_WITH_NAME(colName, mc_patMatch_dr);

    PUSH_VECTOR_WITH_NAME(colName, miniIso_uncor);
    PUSH_VECTOR_WITH_NAME(colName, miniIso_ch);
    PUSH_VECTOR_WITH_NAME(colName, miniIso_nh);
    PUSH_VECTOR_WITH_NAME(colName, miniIso_em);
    PUSH_VECTOR_WITH_NAME(colName, miniIso_db);

  }


  // GenInfo
  if (isMC){
    edm::Handle< GenInfo > genInfo;
    iEvent.getByToken(genInfoToken, genInfo);
    if (!genInfo.isValid()) throw cms::Exception("CMS3Ntuplizer::analyze: Error getting the gen. info. from the event...");

    recordGenInfo(*genInfo);

  }

  // Undefine the convenience macros
#undef PUSH_VECTOR_WITH_NAME
#undef MAKE_VECTOR_WITH_RESERVE

  /************************************************/
  /* Record the communicator values into the tree */
  /************************************************/

  // If this is the first event, create the tree branches based on what is available in the commonEntry.
  if (firstEvent){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) outtree->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) outtree->putBranch(itb->first, &(itb->second));
    SIMPLE_DATA_OUTPUT_DIRECTIVES
    VECTOR_DATA_OUTPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
  }

  // Record whatever is in commonEntry into the tree.
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) outtree->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) outtree->setVal(itb->first, &(itb->second));
  SIMPLE_DATA_OUTPUT_DIRECTIVES
  VECTOR_DATA_OUTPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

  // Fill the tree
  outtree->fill();

  // No longer the first event...
  if (firstEvent) firstEvent = false;
}

void CMS3Ntuplizer::recordGenInfo(GenInfo const& genInfo){

#define SET_GENINFO_VARIABLE(var) commonEntry.setNamedVal(#var, genInfo.var);

  SET_GENINFO_VARIABLE(xsec)
  SET_GENINFO_VARIABLE(xsecerr)

  SET_GENINFO_VARIABLE(qscale)
  SET_GENINFO_VARIABLE(alphaS)

  SET_GENINFO_VARIABLE(genMET)
  SET_GENINFO_VARIABLE(genMETPhi)

  SET_GENINFO_VARIABLE(sumEt)
  SET_GENINFO_VARIABLE(pThat)

  // LHE variations
  SET_GENINFO_VARIABLE(genHEPMCweight_default)
  SET_GENINFO_VARIABLE(genHEPMCweight_NNPDF30)

  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF1)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF2)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF0p5)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF1)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF2)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF0p5)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF1)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF2)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF0p5)

  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Up_default)
  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Dn_default)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Up_default)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Dn_default)

  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Up_NNPDF30)
  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Dn_NNPDF30)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Up_NNPDF30)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Dn_NNPDF30)

  // Pythis PS weights
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muRoneoversqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muRoneoversqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muRsqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muRsqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR0p5)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR0p5)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR2)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR2)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR0p25)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR0p25)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR4)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR4)

#undef SET_GENINFO_VARIABLE

  for (auto const it:genInfo.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);

}


//define this as a plug-in
DEFINE_FWK_MODULE(CMS3Ntuplizer);
// vim: ts=2:sw=2
