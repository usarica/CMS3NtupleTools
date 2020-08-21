#include <limits>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctionsCore.h>
#include <CMS3/NtupleMaker/interface/plugins/GenMaker.h>
#include <CMS3/NtupleMaker/interface/MCUtilities.h>

#include <JHUGenMELA/MELA/interface/PDGHelpers.h>

#include "MELAStreamHelpers.hh"


typedef math::XYZTLorentzVectorF LorentzVector;

using namespace std;
using namespace edm;
using namespace reco;
using namespace MELAStreamHelpers;


GenMaker::GenMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasprefix")),
  year(iConfig.getParameter<int>("year")),

  xsecOverride(static_cast<float>(iConfig.getParameter<double>("xsec"))),
  brOverride(static_cast<float>(iConfig.getParameter<double>("BR"))),

  LHEInputTag_(iConfig.getParameter<edm::InputTag>("LHEInputTag")),
  genEvtInfoInputTag_(iConfig.getParameter<edm::InputTag>("genEvtInfoInputTag")),
  prunedGenParticlesInputTag_(iConfig.getParameter<edm::InputTag>("prunedGenParticlesInputTag")),
  packedGenParticlesInputTag_(iConfig.getParameter<edm::InputTag>("packedGenParticlesInputTag")),
  genJetsInputTag_(iConfig.getParameter<edm::InputTag>("genJetsInputTag")),
  genMETInputTag_(iConfig.getParameter<edm::InputTag>("genMETInputTag")),

  ntuplePackedGenParticles_(iConfig.getParameter<bool>("ntuplePackedGenParticles")),

  superMH(static_cast<float>(iConfig.getParameter<double>("superMH"))),

  doHiggsKinematics(iConfig.getParameter<bool>("doHiggsKinematics")),
  candVVmode(MELAEvent::getCandidateVVModeFromString(iConfig.getUntrackedParameter<string>("candVVmode"))),
  decayVVmode(iConfig.getParameter<int>("decayVVmode")),
  lheMElist(iConfig.getParameter< std::vector<std::string> >("lheMElist")),

  KFactor_QCD_ggVV_Sig_handle(nullptr),
  KFactor_QCD_qqVV_Bkg_handle(nullptr),
  KFactor_EW_qqVV_Bkg_handle(nullptr)
{
  consumesMany<LHEEventProduct>();
  LHERunInfoToken = consumes<LHERunInfoProduct, edm::InRun>(LHEInputTag_);
  genEvtInfoToken = consumes<GenEventInfoProduct>(genEvtInfoInputTag_);
  prunedGenParticlesToken = consumes<reco::GenParticleCollection>(prunedGenParticlesInputTag_);
  packedGenParticlesToken = consumes<pat::PackedGenParticleCollection>(packedGenParticlesInputTag_);
  genJetsToken = consumes< edm::View<reco::GenJet> >(genJetsInputTag_);
  genMETToken = consumes< edm::View<pat::MET> >(genMETInputTag_);

  consumesMany<HepMCProduct>();

  if (year==2016) LHEHandler::set_maxlines_print_header(1000);
  else LHEHandler::set_maxlines_print_header(-1);
  lheHandler_default = std::make_shared<LHEHandler>(
    candVVmode, decayVVmode,
    ((candVVmode!=MELAEvent::nCandidateVVModes && (!lheMElist.empty() || doHiggsKinematics)) ? LHEHandler::doHiggsKinematics : LHEHandler::noKinematics),
    year, LHEHandler::keepDefaultPDF, LHEHandler::keepDefaultQCDOrder
    );
  lheHandler_NNPDF30_NLO = std::make_shared<LHEHandler>(
    MELAEvent::nCandidateVVModes, -1,
    LHEHandler::noKinematics,
    year, LHEHandler::tryNNPDF30, LHEHandler::tryNLO
    );

  // Setup ME computation
  setupMELA();

  // Setup K factor handles
  setupKFactorHandles(iConfig);

  produces<GenInfo>();
}

GenMaker::~GenMaker(){
  // Clear ME computations
  cleanMELA();
}

void GenMaker::beginJob(){}
void GenMaker::endJob(){}

void GenMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& /*iSetup*/){
  static bool firstRun = true;
  if (firstRun){ // Do these only at the first run
    // Extract LHE header
    edm::Handle<LHERunInfoProduct> lhe_runinfo;
    iRun.getByLabel(LHEInputTag_, lhe_runinfo);
    lheHandler_default->setHeaderFromRunInfo(&lhe_runinfo);
    lheHandler_NNPDF30_NLO->setHeaderFromRunInfo(&lhe_runinfo);
    firstRun=false;
  }
}

void GenMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique<GenInfo>();

  edm::Handle<reco::GenParticleCollection> prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesToken, prunedGenParticlesHandle);
  if (!prunedGenParticlesHandle.isValid()){
    edm::LogError("GenMaker") << "GenMaker::produce: Failed to retrieve pruned gen. particle collection...";
    return;
  }
  std::vector<reco::GenParticle> const* prunedGenParticles = prunedGenParticlesHandle.product();

  /*
  edm::Handle<pat::PackedGenParticleCollection> packedGenParticlesHandle;
  iEvent.getByToken(packedGenParticlesToken, packedGenParticlesHandle);
  if (!packedGenParticlesHandle.isValid()){
    edm::LogError("GenMaker") << "GenMaker::produce: Failed to retrieve packed gen. particle collection...";
    return;
  }
  std::vector<pat::PackedGenParticle> const* packedGenParticles = packedGenParticlesHandle.product();
  */

  // Gen. jets
  edm::Handle< edm::View<reco::GenJet> > genJetsHandle;
  iEvent.getByToken(genJetsToken, genJetsHandle);
  if (!genJetsHandle.isValid()) throw cms::Exception("GenMaker::produce: Error getting the gen. jets handle...");

  // Gen. MET
  edm::Handle< edm::View<pat::MET> > genMETHandle;
  iEvent.getByToken(genMETToken, genMETHandle);
  if (!genMETHandle.isValid()) throw cms::Exception("GenMaker::produce: Error getting the gen. MET handle...");
  result->genmet_met = (genMETHandle->front()).genMET()->pt();
  result->genmet_metPhi = (genMETHandle->front()).genMET()->phi();

  // Gen. and LHE weights
  edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
  iEvent.getByToken(genEvtInfoToken, genEvtInfoHandle);
  if (!genEvtInfoHandle.isValid()){
    edm::LogError("GenMaker") << "GenMaker::produce: Failed to retrieve gen. event info...";
    return;
  }
  else{
    result->qscale = genEvtInfoHandle->qScale();
    result->alphaS = genEvtInfoHandle->alphaQCD();
  }

  if (xsecOverride>0.f && brOverride>0.f) result->xsec = xsecOverride*brOverride;
  else if (xsecOverride>0.f) result->xsec = xsecOverride;

  std::vector< edm::Handle<HepMCProduct> > HepMCProductHandles;
  iEvent.getManyByType(HepMCProductHandles);

  std::vector< edm::Handle<LHEEventProduct> > LHEEventInfoList;
  iEvent.getManyByType(LHEEventInfoList);

  float genps_weight = 0;
  if (!HepMCProductHandles.empty()){
    const HepMC::GenEvent* HepMCGenEvent = HepMCProductHandles.front()->GetEvent();

    result->pThat = HepMCGenEvent->event_scale();

    auto HepMCweights = HepMCGenEvent->weights();
    if (!HepMCweights.empty()) genps_weight = HepMCweights.front();
  }

  if (!LHEEventInfoList.empty() && LHEEventInfoList.front().isValid()){
    edm::Handle<LHEEventProduct>& LHEEventInfo = LHEEventInfoList.front();

    lheHandler_default->setHandle(&LHEEventInfo);
    lheHandler_default->extract();
    // Special case for the default:
    // Compute MEs
    if (!lheMElist.empty()){
      MELACandidate* cand = lheHandler_default->getBestCandidate();
      doMELA(cand, *result);
    }
    // Record the LHE-level particles (filled if lheHandler_default->doKinematics>=LHEHandler::doBasicKinematics)
    std::vector<MELAParticle*> const& basicParticleList = lheHandler_default->getParticleList();
    for (auto it_part = basicParticleList.cbegin(); it_part != basicParticleList.cend(); it_part++){
      MELAParticle* part_i = *it_part;
      result->lheparticles_px.push_back(part_i->x());
      result->lheparticles_py.push_back(part_i->y());
      result->lheparticles_pz.push_back(part_i->z());
      result->lheparticles_E.push_back(part_i->t());
      if (part_i->id < std::numeric_limits<cms3_id_t>::min() || part_i->id > std::numeric_limits<cms3_id_t>::max()){
        cms::Exception e("NumericLimits");
        e << "Particle id " << part_i->id << " exceeds numerical bounds.";
        throw e;
      }
      if (part_i->genStatus < std::numeric_limits<cms3_genstatus_t>::min() || part_i->genStatus > std::numeric_limits<cms3_genstatus_t>::max()){
        cms::Exception e("NumericLimits");
        e << "Particle status " << part_i->genStatus << " exceeds numerical bounds.";
        throw e;
      }
      result->lheparticles_id.push_back(part_i->id);
      result->lheparticles_status.push_back(part_i->genStatus);

      int lheparticles_mother0_index=-1; int lheparticles_mother1_index=-1;
      for (int imom=0; imom<std::min(part_i->getNMothers(), 2); imom++){
        MELAParticle* mother = part_i->getMother(imom);
        int& lheparticles_mother_index = (imom==0 ? lheparticles_mother0_index : lheparticles_mother1_index);
        for (auto jt_part = basicParticleList.cbegin(); jt_part != basicParticleList.cend(); jt_part++){
          MELAParticle* part_j = *jt_part;
          if (part_j==mother){
            lheparticles_mother_index = (jt_part - basicParticleList.cbegin());
            break;
          }
        }
      }
      result->lheparticles_mother0_index.push_back(lheparticles_mother0_index);
      result->lheparticles_mother1_index.push_back(lheparticles_mother1_index);
    }


    lheHandler_NNPDF30_NLO->setHandle(&LHEEventInfo);
    lheHandler_NNPDF30_NLO->extract();

    if (genEvtInfoHandle.isValid()){
      result->genHEPMCweight_default = genEvtInfoHandle->weight();
      if (result->genHEPMCweight_default==1.) result->genHEPMCweight_default = lheHandler_default->getLHEOriginalWeight();
    }
    else result->genHEPMCweight_default = lheHandler_default->getLHEOriginalWeight(); // Default was also 1, so if !genEvtInfoHandle.isValid(), the statement still passes
    result->genHEPMCweight_NNPDF30 = result->genHEPMCweight_default; // lheHandler_NNPDF30_NLO->getLHEOriginalWeight() should give the same value
    result->genHEPMCweight_default *= lheHandler_default->getWeightRescale();
    result->genHEPMCweight_NNPDF30 *= lheHandler_NNPDF30_NLO->getWeightRescale();

    result->LHEweight_scaledOriginalWeight_default = result->LHEweight_scaledOriginalWeight_NNPDF30 = lheHandler_default->getLHEOriginalWeight();
    result->LHEweight_scaledOriginalWeight_default *= lheHandler_default->getWeightRescale();
    result->LHEweight_scaledOriginalWeight_NNPDF30 *= lheHandler_NNPDF30_NLO->getWeightRescale();

    result->LHEweight_defaultMemberZero = lheHandler_default->getMemberZeroWeight(); // Same as lheHandler_NNPDF30_NLO->getMemberZeroWeight()

    result->LHEweight_QCDscale_muR1_muF1 = lheHandler_default->getLHEWeight(0, 1.);
    result->LHEweight_QCDscale_muR1_muF2 = lheHandler_default->getLHEWeight(1, 1.);
    result->LHEweight_QCDscale_muR1_muF0p5 = lheHandler_default->getLHEWeight(2, 1.);
    result->LHEweight_QCDscale_muR2_muF1 = lheHandler_default->getLHEWeight(3, 1.);
    result->LHEweight_QCDscale_muR2_muF2 = lheHandler_default->getLHEWeight(4, 1.);
    result->LHEweight_QCDscale_muR2_muF0p5 = lheHandler_default->getLHEWeight(5, 1.);
    result->LHEweight_QCDscale_muR0p5_muF1 = lheHandler_default->getLHEWeight(6, 1.);
    result->LHEweight_QCDscale_muR0p5_muF2 = lheHandler_default->getLHEWeight(7, 1.);
    result->LHEweight_QCDscale_muR0p5_muF0p5 = lheHandler_default->getLHEWeight(8, 1.);

    result->LHEweight_PDFVariation_Up_default = lheHandler_default->getLHEWeight_PDFVariationUpDn(1, 1.);
    result->LHEweight_PDFVariation_Dn_default = lheHandler_default->getLHEWeight_PDFVariationUpDn(-1, 1.);
    result->LHEweight_AsMZ_Up_default = lheHandler_default->getLHEWeigh_AsMZUpDn(1, 1.);
    result->LHEweight_AsMZ_Dn_default = lheHandler_default->getLHEWeigh_AsMZUpDn(-1, 1.);

    result->LHEweight_PDFVariation_Up_NNPDF30 = lheHandler_NNPDF30_NLO->getLHEWeight_PDFVariationUpDn(1, 1.);
    result->LHEweight_PDFVariation_Dn_NNPDF30 = lheHandler_NNPDF30_NLO->getLHEWeight_PDFVariationUpDn(-1, 1.);
    result->LHEweight_AsMZ_Up_NNPDF30 = lheHandler_NNPDF30_NLO->getLHEWeigh_AsMZUpDn(1, 1.);
    result->LHEweight_AsMZ_Dn_NNPDF30 = lheHandler_NNPDF30_NLO->getLHEWeigh_AsMZUpDn(-1, 1.);

    result->xsec_lhe = LHEEventInfo->originalXWGTUP();

    lheHandler_default->clear();
    lheHandler_NNPDF30_NLO->clear();
  }
  else{
    if (genEvtInfoHandle.isValid()) result->genHEPMCweight_default = genEvtInfoHandle->weight();
    else result->genHEPMCweight_default = genps_weight;
    result->genHEPMCweight_NNPDF30 = result->genHEPMCweight_default; // No other choice really
  }

  // Extract PS weights
  {
    const auto& genweights = genEvtInfoHandle->weights();
    if (genweights.size() > 1){
      if ((genweights.size() != 14 && genweights.size() != 46) || genweights.at(0) != genweights.at(1)){
        cms::Exception e("GenWeights");
        e << "Expected to find 1 gen weight, or 14 or 46 with the first two the same, found " << genweights.size() << ":\n";
        for (auto w : genweights) e << w << " ";
        throw e;
      }
      auto const& nominal = genweights[0];
      result->PythiaWeight_isr_muRoneoversqrt2 = genweights[2] / nominal;
      result->PythiaWeight_fsr_muRoneoversqrt2 = genweights[3] / nominal;
      result->PythiaWeight_isr_muRsqrt2 = genweights[4] / nominal;
      result->PythiaWeight_fsr_muRsqrt2 = genweights[5] / nominal;

      result->PythiaWeight_isr_muR0p5 = genweights[6] / nominal;
      result->PythiaWeight_fsr_muR0p5 = genweights[7] / nominal;
      result->PythiaWeight_isr_muR2 = genweights[8] / nominal;
      result->PythiaWeight_fsr_muR2 = genweights[9] / nominal;

      result->PythiaWeight_isr_muR0p25 = genweights[10] / nominal;
      result->PythiaWeight_fsr_muR0p25 = genweights[11] / nominal;
      result->PythiaWeight_isr_muR4 = genweights[12] / nominal;
      result->PythiaWeight_fsr_muR4 = genweights[13] / nominal;
    }
  }

  {
    auto& n_shower_gluons_to_bottom = result->n_shower_gluons_to_bottom; n_shower_gluons_to_bottom = 0;
    auto& n_shower_gluons_to_charm = result->n_shower_gluons_to_charm; n_shower_gluons_to_charm = 0;
    std::vector<reco::GenParticle const*> shower_gluons_to_bottom;
    std::vector<reco::GenParticle const*> shower_gluons_to_bottom_daughters;
    std::vector<reco::GenParticle const*> shower_gluons_to_charm;
    std::vector<reco::GenParticle const*> shower_gluons_to_charm_daughters;

    float& sumEt = result->sumEt; sumEt=0;
    LorentzVector tempvect(0, 0, 0, 0);

    // Compute HT and MHT from hard gen. partons
    float& genhardpartons_HT = result->genhardpartons_HT;
    float& genhardpartons_MHT = result->genhardpartons_MHT;
    genhardpartons_HT = genhardpartons_MHT = 0;
    LorentzVector tempvect_hardpartons(0, 0, 0, 0);

    for (std::vector<reco::GenParticle>::const_iterator genps_it = prunedGenParticles->begin(); genps_it != prunedGenParticles->end(); genps_it++){
      auto const* genps = &(*genps_it);
      int id = genps->pdgId();
      if (PDGHelpers::isANeutrino(id) && genps->status()==1) tempvect += genps->p4();
      if (PDGHelpers::isAKnownJet(id) && genps->isHardProcess() && (genps->status()==1 || genps->status()==23 || genps->status()==24)){
        tempvect_hardpartons += genps->p4();
        genhardpartons_HT += genps->pt();
      }
      if (PDGHelpers::isAGluon(id) && !genps->isHardProcess()){
        std::vector<reco::GenParticle const*> tmp_daughters;
        MCUtilities::getAllDaughters(genps, tmp_daughters, false);
        for (auto const& dau:tmp_daughters){
          // Skip the duaghter if it is already examined
          if (
            HelperFunctions::checkListVariable(shower_gluons_to_bottom_daughters, dau)
            ||
            HelperFunctions::checkListVariable(shower_gluons_to_charm_daughters, dau)
            ) continue;

          unsigned int const dau_id = std::abs(dau->pdgId());
          if (
            // Bottom baryons
            (dau_id>=5000 && dau_id<6000)
            ||
            // Bottom mesons
            (dau_id>=500 && dau_id<600)
            ||
            (dau_id>10000 && (dau_id/100 % 10 == 5))
            ){
            if (!HelperFunctions::checkListVariable(shower_gluons_to_bottom, genps)){
              shower_gluons_to_bottom.push_back(genps);
              shower_gluons_to_bottom_daughters.push_back(dau);
              //MELAout << "\t- Adding daughter " << dau->pdgId() << " with p4 = " << dau->p4() << " mother gluon p4 = " << genps->p4() << endl;
            }
          }
          else if (
            // Charm baryons
            (dau_id>=4000 && dau_id<5000)
            ||
            // Charm mesons
            (dau_id>=400 && dau_id<500)
            ||
            (dau_id>10000 && (dau_id/100 % 10 == 4))
            ){
            if (!HelperFunctions::checkListVariable(shower_gluons_to_charm, genps)){
              shower_gluons_to_charm.push_back(genps);
              shower_gluons_to_charm_daughters.push_back(dau);
              //MELAout << "\t- Adding daughter " << dau->pdgId() << " with p4 = " << dau->p4() << " mother gluon p4 = " << genps->p4() << endl;
            }
          }
          //else if (dau_id==21) MELAout << "\t- Another daughter is a gluon" << endl;
        } // End loop over gluon daughters
      } // End if-statemetn for gluons
    } // End loop over gen particles

    sumEt = tempvect.pt();
    genhardpartons_MHT = tempvect.Pt();

    n_shower_gluons_to_bottom = shower_gluons_to_bottom.size();
    n_shower_gluons_to_charm = shower_gluons_to_charm.size();
    //MELAout << "Number of g->bb, g->cc: " << n_shower_gluons_to_bottom << ", " << n_shower_gluons_to_charm << endl;
  }

  // Compute gen. jet HT and MHT
  {
    result->genjets_HT = result->genjets_MHT = 0;
    LorentzVector tempvect(0, 0, 0, 0);
    for (edm::View<reco::GenJet>::const_iterator genjet_it = genJetsHandle->begin(); genjet_it != genJetsHandle->end(); genjet_it++){
      tempvect += genjet_it->p4();
      result->genjets_HT += genjet_it->pt();
    }
    result->genjets_MHT = tempvect.Pt();
  }

  // Compute K factors
  if (KFactor_QCD_ggVV_Sig_handle){
    //std::vector<reco::GenParticle const*> hardparticles;
    std::vector<reco::GenParticle const*> Higgses;
    std::vector<reco::GenParticle const*> Vbosons;

    for (std::vector<reco::GenParticle>::const_iterator it_part = prunedGenParticles->begin(); it_part != prunedGenParticles->end(); it_part++){
      reco::GenParticle const* part = &(*it_part);
      if (part->isHardProcess()){
        if (PDGHelpers::isAHiggs(part->pdgId())) Higgses.push_back(part);
        else if (PDGHelpers::isAZBoson(part->pdgId()) || PDGHelpers::isAWBoson(part->pdgId()) || PDGHelpers::isAPhoton(part->pdgId())) Vbosons.push_back(part);
        //hardparticles.push_back(part);
      }
    }

    if (Higgses.size()==1) result->Kfactors[KFactorHelpers::KFactorHandler_QCD_ggVV_Sig::KFactorArgName] = Higgses.front()->mass();
    else if (Vbosons.size()==2) result->Kfactors[KFactorHelpers::KFactorHandler_QCD_ggVV_Sig::KFactorArgName] = (Vbosons.front()->p4() + Vbosons.back()->p4()).M();
    else throw cms::Exception("GenMaker::produce: No single Higgs candidate or two intermediate V bosons are found to pass to KFactor_QCD_ggVV_Sig_handle.");
    for (auto const& kfpair:kfactor_num_denum_list){
      if (kfpair.first == KFactorHelpers::kf_QCD_NNLO_GGVV_SIG || kfpair.first == KFactorHelpers::kf_QCD_NLO_GGVV_SIG) KFactor_QCD_ggVV_Sig_handle->eval(kfpair.first, kfpair.second, result->Kfactors);
    }
  }
  if (KFactor_QCD_qqVV_Bkg_handle || KFactor_EW_qqVV_Bkg_handle){
    KFactorHelpers::KFactorType corr_type = KFactorHelpers::nKFactorTypes;
    if (KFactor_EW_qqVV_Bkg_handle) corr_type = KFactor_EW_qqVV_Bkg_handle->getType();
    else{
      for (auto const& kfpair:kfactor_num_denum_list){
        if (
          kfpair.first == KFactorHelpers::kf_QCD_NNLO_QQZZ_BKG
          ||
          kfpair.first == KFactorHelpers::kf_QCD_NNLO_QQWZ_BKG
          ||
          kfpair.first == KFactorHelpers::kf_QCD_NNLO_QQWW_BKG
          ){
          corr_type = kfpair.first;
          break;
        }
      }
    }

    KFactorHelpers::VVFinalStateType final_state_type = KFactorHelpers::nVVFinalStateTypes;
    switch (corr_type){
    case KFactorHelpers::kf_QCD_NNLO_QQZZ_BKG:
    case KFactorHelpers::kf_EW_NLO_QQZZ_BKG:
      final_state_type = KFactorHelpers::kZZ;
      break;
    case KFactorHelpers::kf_QCD_NNLO_QQWZ_BKG:
    case KFactorHelpers::kf_EW_NLO_QQWZ_BKG:
      final_state_type = KFactorHelpers::kWZ;
      break;
    case KFactorHelpers::kf_QCD_NNLO_QQWW_BKG:
    case KFactorHelpers::kf_EW_NLO_QQWW_BKG:
      final_state_type = KFactorHelpers::kWW;
      break;
    default:
      throw cms::Exception("UnknownFinalState") << "GenMaker::produce: No known VV final state for the correction type " << corr_type << ".";
      break;
    }

    //MELAout << "Extracting LHE info" << endl;
    std::vector<reco::GenParticle const*> incomingQuarks;
    std::vector<reco::GenParticle const*> incomingGluons;
    std::vector<reco::GenParticle const*> outgoingQuarks;
    std::vector<reco::GenParticle const*> outgoingGluons;
    std::pair<reco::GenParticle const*, reco::GenParticle const*> V1pair;
    std::pair<reco::GenParticle const*, reco::GenParticle const*> V2pair;
    KFactorHelpers::getVVTopology(
      final_state_type, *prunedGenParticles,
      incomingQuarks, incomingGluons,
      outgoingQuarks, outgoingGluons,
      V1pair, V2pair
    );

    for (auto const& kfpair:kfactor_num_denum_list){
      //MELAout << "Calculating K factor " << kfpair.first << ", " << kfpair.second << endl;
      if (
        kfpair.first == KFactorHelpers::kf_QCD_NNLO_QQZZ_BKG
        ||
        kfpair.first == KFactorHelpers::kf_QCD_NNLO_QQWZ_BKG
        ||
        kfpair.first == KFactorHelpers::kf_QCD_NNLO_QQWW_BKG
        ) KFactor_QCD_qqVV_Bkg_handle->eval(kfpair.first, V1pair, V2pair, result->Kfactors);
      if (
        genEvtInfoHandle.isValid()
        && (
          kfpair.first == KFactorHelpers::kf_EW_NLO_QQZZ_BKG
          ||
          kfpair.first == KFactorHelpers::kf_EW_NLO_QQWZ_BKG
          ||
          kfpair.first == KFactorHelpers::kf_EW_NLO_QQWW_BKG
          )
        ) KFactor_EW_qqVV_Bkg_handle->eval(
          *genEvtInfoHandle,
          incomingQuarks, V1pair, V2pair,
          result->Kfactors
        );
      //MELAout << "\t- Done" << endl;
    }
  }

  iEvent.put(std::move(result));
}


/******************/
/* ME COMPUTATION */
/******************/
void GenMaker::setupMELA(){
  if (candVVmode==MELAEvent::nCandidateVVModes || lheMElist.empty()){
    this->usesResource("GenMakerNoMELA");
    return;
  }
  this->usesResource("MELA");

  using namespace CMS3MELAHelpers;

  setupMela(year, superMH, TVar::ERROR); // Sets up MELA only once

  lheMEblock.buildMELABranches(lheMElist, true);
}
void GenMaker::doMELA(MELACandidate* cand, GenInfo& genInfo){
  using namespace CMS3MELAHelpers;
  if (melaHandle && cand){
    melaHandle->setCurrentCandidate(cand);

    lheMEblock.computeMELABranches();
    lheMEblock.pushMELABranches();
    lheMEblock.getBranchValues(genInfo.LHE_ME_weights); // Record the MEs into the EDProducer product
    genInfo.LHE_ME_weights["LHECandMass"] = cand->m(); // When LHE MEs are present, you must include this variable to combine different samples.

    melaHandle->resetInputEvent();
  }
  else if (melaHandle && !cand){
    MELAout
      << "GenMaker::doMELA: Default LHEHandler options were\n"
      << "candVVmode: " << candVVmode << '\n'
      << "decayVVmode: " << decayVVmode << '\n'
      << "kinFlag: " << ((candVVmode!=MELAEvent::nCandidateVVModes && (!lheMElist.empty() || doHiggsKinematics)) ? LHEHandler::doHiggsKinematics : LHEHandler::noKinematics) << '\n'
      << "year: " << year << '\n'
      << "but MELACandidate construction failed."
      << endl;
    throw cms::Exception("GenMaker::doMELA: The MELA handle is set but the candidate is null.");
  }
}
void GenMaker::cleanMELA(){
  using namespace CMS3MELAHelpers;
  // Shared pointer should be able to clear itself
  //clearMela();
}


/************************/
/* K FACTOR COMPUTATION */
/************************/
void GenMaker::setupKFactorHandles(edm::ParameterSet const& iConfig){
  edm::VParameterSet kfactorsets = iConfig.getParameter<edm::VParameterSet>("kfactors");
  if (!kfactorsets.empty()) kfactor_num_denum_list.reserve(kfactorsets.size());
  for (edm::ParameterSet const& pset:kfactorsets){
    // Get K factor specifications
    std::string strkfactor_num = pset.getParameter<std::string>("numerator");
    std::string strkfactor_den;
    if (pset.exists("denominator")) strkfactor_den = pset.getParameter<std::string>("denominator");
    // Use lowercase letters for comparison
    std::string strkfactor_num_lower, strkfactor_den_lower;
    HelperFunctions::lowercase(strkfactor_num, strkfactor_num_lower);
    HelperFunctions::lowercase(strkfactor_den, strkfactor_den_lower);

    // Build K factors
    KFactorHelpers::KFactorType numerator = KFactorHelpers::nKFactorTypes;
    KFactorHelpers::KFactorType denominator = KFactorHelpers::nKFactorTypes;
    if (strkfactor_num_lower == "kfactor_qcd_nnlo_ggvv_sig") numerator = KFactorHelpers::kf_QCD_NNLO_GGVV_SIG;
    else if (strkfactor_num_lower == "kfactor_qcd_nlo_ggvv_sig") numerator = KFactorHelpers::kf_QCD_NLO_GGVV_SIG;
    else if (strkfactor_num_lower == "kfactor_qcd_nnlo_qqzz_bkg") numerator = KFactorHelpers::kf_QCD_NNLO_QQZZ_BKG;
    else if (strkfactor_num_lower == "kfactor_qcd_nnlo_qqwz_bkg") numerator = KFactorHelpers::kf_QCD_NNLO_QQWZ_BKG;
    else if (strkfactor_num_lower == "kfactor_qcd_nnlo_qqww_bkg") numerator = KFactorHelpers::kf_QCD_NNLO_QQWW_BKG;
    else if (strkfactor_num_lower == "kfactor_ew_nlo_qqzz_bkg") numerator = KFactorHelpers::kf_EW_NLO_QQZZ_BKG;
    else if (strkfactor_num_lower == "kfactor_ew_nlo_qqwz_bkg") numerator = KFactorHelpers::kf_EW_NLO_QQWZ_BKG;
    else if (strkfactor_num_lower == "kfactor_ew_nlo_qqww_bkg") numerator = KFactorHelpers::kf_EW_NLO_QQWW_BKG;
    else throw cms::Exception(Form("GenMaker::setupKFactorHandles: Cannot identify the numerator of the K factor pair (%s, %s).", strkfactor_num.data(), strkfactor_den.data()));

    bool doBuild_KFactor_QCD_ggVV_Sig_handle = (numerator==KFactorHelpers::kf_QCD_NNLO_GGVV_SIG || numerator==KFactorHelpers::kf_QCD_NLO_GGVV_SIG);
    if (doBuild_KFactor_QCD_ggVV_Sig_handle){
      if (strkfactor_den_lower == "kfactor_qcd_nnlo_ggvv_sig") denominator = KFactorHelpers::kf_QCD_NNLO_GGVV_SIG;
      else if (strkfactor_den_lower == "kfactor_qcd_nlo_ggvv_sig") denominator = KFactorHelpers::kf_QCD_NLO_GGVV_SIG;
      else throw cms::Exception(Form("GenMaker::setupKFactorHandles: K factor pair (%s, %s) is not implemented.", strkfactor_num.data(), strkfactor_den.data()));
      //MELAout << "GenMaker::setupKFactorHandles: Building ggVV QCD K factors with num=" << strkfactor_num << " and den=" << strkfactor_den << endl;
      KFactor_QCD_ggVV_Sig_handle = std::make_shared<KFactorHelpers::KFactorHandler_QCD_ggVV_Sig>(this->year);
      //MELAout << "\t- Built!" << endl;
    }

    bool doBuild_KFactor_QCD_qqVV_Bkg_handle = (numerator==KFactorHelpers::kf_QCD_NNLO_QQZZ_BKG || numerator==KFactorHelpers::kf_QCD_NNLO_QQWZ_BKG || numerator==KFactorHelpers::kf_QCD_NNLO_QQWW_BKG);
    if (doBuild_KFactor_QCD_qqVV_Bkg_handle){
      KFactor_QCD_qqVV_Bkg_handle = std::make_shared<KFactorHelpers::KFactorHandler_QCD_qqVV_Bkg>(this->year);
    }

    bool doBuild_KFactor_EW_qqVV_Bkg_handle = (numerator==KFactorHelpers::kf_EW_NLO_QQZZ_BKG || numerator==KFactorHelpers::kf_EW_NLO_QQWZ_BKG || numerator==KFactorHelpers::kf_EW_NLO_QQWW_BKG);
    if (doBuild_KFactor_EW_qqVV_Bkg_handle){
      //MELAout << "GenMaker::setupKFactorHandles: Building qqVV EW K factors with num=" << strkfactor_num << " and den=" << strkfactor_den << endl;
      KFactor_EW_qqVV_Bkg_handle = std::make_shared<KFactorHelpers::KFactorHandler_EW_qqVV_Bkg>(this->year, numerator);
      //MELAout << "\t- Built!" << endl;
    }

    if (numerator!=KFactorHelpers::nKFactorTypes) kfactor_num_denum_list.emplace_back(numerator, denominator);
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenMaker);
