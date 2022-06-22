#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
#include "PlottingHelpers.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TRandom3.h"


// Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  /* Pythia weight adjustments are independent of PDF choice */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_PythiaScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PythiaScaleUp) \
  /* Factorization and renormalization scale weight adjustments are independent of PDF choice (because they are only done for the default PDF set) */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_QCDScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_QCDScaleUp) \
  /* a_s(mZ) and PDF replica weight adjustments come from the specific PDF set, so they are split between 'default' vs 'NNPDF3.0'. */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFReplicaUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaUp) \
  /* External adjustment factors from the registered histograms */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_HardJetsDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_HardJetsUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PythiaScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PythiaScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_QCDScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_QCDScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFReplicaUp) \
  /* ME adjustment weights */ \
  BRANCH_COMMAND(bool, invalidReweightingWgts) \
  BRANCH_COMMAND(float, sample_wgt) \
  BRANCH_COMMAND(float, sample_wgt_pairwiseComponent) \
  /* Other variables */ \
  BRANCH_COMMAND(bool, passGenCompSelection) \
  BRANCH_COMMAND(float, lheHiggs_mass) \
  BRANCH_COMMAND(float, lheHiggs_pt) \
  BRANCH_COMMAND(float, lheLeptonicDecay_pt) \
  BRANCH_COMMAND(float, lheLeptonicDecay_mass) \
  BRANCH_COMMAND(float, genLeptonicDecay_pt) \
  BRANCH_COMMAND(float, genLeptonicDecay_mass)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, lheparticles_id) \
  BRANCH_COMMAND(cms3_id_t, lheparticles_status) \
  BRANCH_COMMAND(cms3_id_t, lheQCDLOparticles_id) \
  BRANCH_COMMAND(cms3_id_t, lheQCDLOparticles_status) \
  BRANCH_COMMAND(cms3_id_t, genparticles_id) \
  BRANCH_COMMAND(float, lheQCDLOparticles_px) \
  BRANCH_COMMAND(float, lheQCDLOparticles_py) \
  BRANCH_COMMAND(float, lheQCDLOparticles_pz) \
  BRANCH_COMMAND(float, lheQCDLOparticles_E) \
  BRANCH_COMMAND(float, lheparticles_px) \
  BRANCH_COMMAND(float, lheparticles_py) \
  BRANCH_COMMAND(float, lheparticles_pz) \
  BRANCH_COMMAND(float, lheparticles_E) \
  BRANCH_COMMAND(float, genparticles_pt) \
  BRANCH_COMMAND(float, genparticles_eta) \
  BRANCH_COMMAND(float, genparticles_phi) \
  BRANCH_COMMAND(float, genparticles_mass) \
  BRANCH_COMMAND(float, genak4jets_pt) \
  BRANCH_COMMAND(float, genak4jets_eta) \
  BRANCH_COMMAND(float, genak4jets_phi) \
  BRANCH_COMMAND(float, genak4jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


std::vector<std::string> getMEList(){
  return std::vector<std::string>{
    "Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen Options:AddPConst=1",
    "Name:GG_SIG_kappaTopBot_1_ghz1_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Options:AddPConst=1",
    "Name:GG_BSI_kappaTopBot_1_ghz1_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM",
    "Name:GG_BKG_MCFM Process:bkgZZ Production:ZZGG MatrixElement:MCFM Options:AddPConst=1",
    "Name:QQB_BKG_MCFM Alias:<Name> Process:bkgZZ Production:ZZQQB MatrixElement:MCFM Options:AddPConst=1",
    "Name:JJVBF_SIG_ghv1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:JJQCD_SIG_ghg2_1_JHUGen Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:HadZH_SIG_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:HadWH_SIG_ghw1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

    "Name:JJVBF_S_SIG_ghv1_1_MCFM Process:HSMHiggs Production:JJVBF_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:JJVBF_S_BSI_ghv1_1_MCFM Process:bkgZZ_SMHiggs Production:JJVBF_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
    "Name:JJVBF_BKG_MCFM Process:bkgZZ Production:JJVBF MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

    "Name:HadZH_S_SIG_ghz1_1_MCFM Process:HSMHiggs Production:Had_ZH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:HadZH_S_BSI_ghz1_1_MCFM Process:bkgZZ_SMHiggs Production:Had_ZH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
    "Name:HadZH_BKG_MCFM Process:bkgZZ Production:Had_ZH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

    "Name:HadWH_S_SIG_ghw1_1_MCFM Process:HSMHiggs Production:Had_WH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:HadWH_S_BSI_ghw1_1_MCFM Process:bkgZZ_SMHiggs Production:Had_WH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
    "Name:HadWH_BKG_MCFM Process:bkgZZ Production:Had_WH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

    "Name:JJQCD_BKG_MCFM Process:bkgZZ Production:JJQCD MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

    "Name:HadZH_mavjj_true Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJTrue:1",
    "Name:HadWH_mavjj_true Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJTrue:1",
    "Name:HadZH_mavjj Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJ:1",
    "Name:HadWH_mavjj Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJ:1"
  };
}


void compute4lCorrection(
  TString strSampleSet, TString period, TString strdate
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  bool const isPowheg = strSampleSet.Contains("POWHEG");
  bool const acquireMEs = isPowheg;

  // Set output directory
  TString cinput_main = "output/ReweightedGenTrees/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  curdir->cd();

  SampleHelpers::configure(period, Form("store:201221_%s", period.Data()));

  IvyMELAHelpers::setupMela(SampleHelpers::getDataYear(), 125., TVar::ERROR);
  std::vector<std::string> recoMElist = getMEList();
  IvyMELAHelpers::GMECBlock MEblock;

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;
  PUJetIdScaleFactorHandler pujetidSFHandler;

  TFile* foutput = TFile::Open(strSampleSet + "_corrhists.root", "recreate");
  TTree* tout_correct = new TTree("SkimTree_Correct", "");
  TTree* tout_wrong = new TTree("SkimTree_Wrong", "");
  MEblock.addRefTree(tout_correct);
  MEblock.addRefTree(tout_wrong);
  MEblock.buildMELABranches(recoMElist, false);

  float wgt_bkg; tout_correct->Branch("wgt_bkg", &wgt_bkg); tout_wrong->Branch("wgt_bkg", &wgt_bkg);
  float wgt_sig; tout_correct->Branch("wgt_sig", &wgt_sig); tout_wrong->Branch("wgt_sig", &wgt_sig);
  float wgt_bsi; tout_correct->Branch("wgt_bsi", &wgt_bsi); tout_wrong->Branch("wgt_bsi", &wgt_bsi);
  float mZZ; tout_correct->Branch("mZZ", &mZZ); tout_wrong->Branch("mZZ", &mZZ);
  float mZ1; tout_correct->Branch("mZ1", &mZ1); tout_wrong->Branch("mZ1", &mZ1);
  float mZ2; tout_correct->Branch("mZ2", &mZ2); tout_wrong->Branch("mZ2", &mZ2);
  int id_Z1; tout_correct->Branch("id_Z1", &id_Z1); tout_wrong->Branch("id_Z1", &id_Z1);
  int id_Z2; tout_correct->Branch("id_Z2", &id_Z2); tout_wrong->Branch("id_Z2", &id_Z2);
  unsigned int nleps_reco; tout_correct->Branch("nleps_reco", &nleps_reco); tout_wrong->Branch("nleps_reco", &nleps_reco);
  unsigned int njets; tout_correct->Branch("njets", &njets); tout_wrong->Branch("njets", &njets);

  curdir->cd();

  // Create the output file
  TString cinput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(cinput, "_MINIAOD", "");
  TString strinput = cinput_main + "/" + cinput + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
  IVYout << "Acquiring input " << strinput << "..." << endl;
  BaseTree* tin = new BaseTree(strinput, "SkimTree", "", "");
  if (!tin->isValid()){
    IVYerr << "\t- Failed to acquire." << endl;
    delete tin;
    exit(1);
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE* NAME = nullptr; tin->bookBranch<TYPE>(#NAME, 0);
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>** NAME=nullptr; tin->bookBranch<std::vector<TYPE>*>(#NAME, nullptr);
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) tin->getValRef(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float*> ME_Kfactor_values;
  std::vector<TString> allbranchnames;
  tin->getValidBranchNamesWithoutAlias(allbranchnames, false);
  for (auto const& bname:allbranchnames){
    if (
      (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
      ||
      bname.BeginsWith("p_Gen")
      ||
      bname.Contains("LHECandMass")
      ||
      bname.BeginsWith("KFactor")
      ){
      tin->bookBranch<float>(bname, -1.f);
      ME_Kfactor_values[bname] = nullptr;
      tin->getValRef(bname, ME_Kfactor_values[bname]);
    }
  }

  tin->silenceUnused();

  float* val_Kfactor_QCD = nullptr;
  float* val_ME_BKG = nullptr;
  float* val_ME_SIG = nullptr;
  float* val_ME_BSI = nullptr;
  float* val_ME_CPS = nullptr;
  if (isGG){
    val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second;
    if (acquireMEs){
      val_ME_BKG = ME_Kfactor_values.find("p_Gen_GG_BKG_MCFM")->second;
      val_ME_SIG = ME_Kfactor_values.find("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM")->second;
      val_ME_BSI = ME_Kfactor_values.find("p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM")->second;
      val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
    }
  }
  else{
    if (acquireMEs){
      val_ME_BKG = ME_Kfactor_values.find("p_Gen_JJEW_BKG_MCFM")->second;
      val_ME_SIG = ME_Kfactor_values.find("p_Gen_JJEW_SIG_ghv1_1_MCFM")->second;
      val_ME_BSI = ME_Kfactor_values.find("p_Gen_JJEW_BSI_ghv1_1_MCFM")->second;
      val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
    }
  }

  int nEntries = tin->getNEvents();
  IVYout << "Looping over " << nEntries << " events:" << endl;
  for (int ev=0; ev<nEntries; ev++){
    if (SampleHelpers::doSignalInterrupt==1) break;

    tin->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    wgt_bkg =
      (*event_wgt) * (*event_wgt_adjustment_NNPDF30)
      * (*sample_wgt) * (*invalidReweightingWgts ? 0.f : 1.f)
      * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
      * (val_ME_BKG ? *val_ME_BKG : 1.f) * (val_ME_CPS ? *val_ME_CPS : 1.f);
    wgt_sig =
      (*event_wgt) * (*event_wgt_adjustment_NNPDF30)
      * (*sample_wgt) * (*invalidReweightingWgts ? 0.f : 1.f)
      * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
      * (val_ME_SIG ? *val_ME_SIG : 1.f) * (val_ME_CPS ? *val_ME_CPS : 1.f);
    wgt_bsi =
      (*event_wgt) * (*event_wgt_adjustment_NNPDF30)
      * (*sample_wgt) * (*invalidReweightingWgts ? 0.f : 1.f)
      * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
      * (val_ME_BSI ? *val_ME_BSI : 1.f) * (val_ME_CPS ? *val_ME_CPS : 1.f);

    unsigned char nleps_pt10 = 0;
    unsigned char nleps_pt20 = 0;
    std::vector<ParticleObject> mothers; mothers.reserve(2);
    std::vector<ParticleObject> leptons; leptons.reserve((*lheQCDLOparticles_id)->size());
    std::vector<ParticleObject> leptons_all; leptons_all.reserve((*lheQCDLOparticles_id)->size());
    std::vector<ParticleObject> quarks; quarks.reserve((*lheQCDLOparticles_id)->size());
    for (unsigned int ip=0; ip<(*lheQCDLOparticles_id)->size(); ip++){
      auto const& st = (*lheQCDLOparticles_status)->at(ip);
      auto const& pid = (*lheQCDLOparticles_id)->at(ip);
      ParticleObject::LorentzVector_t pp4((*lheQCDLOparticles_px)->at(ip), (*lheQCDLOparticles_py)->at(ip), (*lheQCDLOparticles_pz)->at(ip), (*lheQCDLOparticles_E)->at(ip));

      if (st<0){
        mothers.emplace_back(pid, pp4);
        continue;
      }

      if (PDGHelpers::isALepton(pid) || PDGHelpers::isANeutrino(pid)){
        leptons_all.emplace_back(pid, pp4);
        if (!(std::abs(pid)==15 || PDGHelpers::isANeutrino(pid))){
          bool const isMuon = (std::abs(pid)==13);
          double pt = pp4.Pt();
          double eta = pp4.Eta();
          if (pt<(isMuon ? 5. : 7.) || std::abs(eta)>=(isMuon ? 2.4 : 2.5)) continue;

          if (pt>=10.){
            nleps_pt10++;
            if (pt>=20.){
              nleps_pt20++;
            }
          }

          double phi = pp4.Phi();
          float dummy_SF = 1;
          float theEff_lnt = 1;
          float theEff_t = 1;
          float theEff = 1;
          if (isMuon){
            muonSFHandler.getIdIsoSFAndEff(theGlobalSyst, static_cast<float>(pt), static_cast<float>(eta), true, true, true, dummy_SF, &theEff_t);
            muonSFHandler.getIdIsoSFAndEff(theGlobalSyst, static_cast<float>(pt), static_cast<float>(eta), true, true, false, dummy_SF, &theEff_lnt);

          }
          else{
            electronSFHandler.getIdIsoSFAndEff(theGlobalSyst, static_cast<float>(pt), static_cast<float>(eta), 2, true, true, true, dummy_SF, &theEff_t);
            electronSFHandler.getIdIsoSFAndEff(theGlobalSyst, static_cast<float>(pt), static_cast<float>(eta), 2, true, true, false, dummy_SF, &theEff_lnt);
          }
          theEff = theEff_lnt + theEff_t;
          if (theEff>1.) cerr << "Eff = " << theEff << ">1" << endl;
          TRandom3 rand;
          rand.SetSeed(static_cast<unsigned long long>(std::abs(std::sin(phi)*100000.)) + static_cast<unsigned long long>(ev));
          double rndnum = rand.Uniform();
          if (rndnum>theEff) continue;

          leptons.emplace_back(pid, pp4);
        }
      }
      else if (PDGHelpers::isAJet(pid)) quarks.emplace_back(pid, pp4);
    }

    nleps_reco = leptons.size();
    if (nleps_reco<4) continue;
    if (!(nleps_pt10>=2 && nleps_pt20>=1)) continue;

    bool has_Hpeak = false;
    {
      std::vector<std::pair<ParticleObject*, ParticleObject*>> dilepton_pairs;
      std::vector<std::pair<ParticleObject*, ParticleObject*>> dineutrino_pairs;
      std::vector<std::pair<ParticleObject*, ParticleObject*>> diquark_pairs;
      for (auto it_l1=leptons_all.begin(); it_l1!=leptons_all.end();it_l1++){
        for (auto it_l2=it_l1+1; it_l2!=leptons_all.end(); it_l2++){
          if (it_l1->pdgId()+it_l2->pdgId()==0){
            if (PDGHelpers::isANeutrino(it_l1->pdgId())) dineutrino_pairs.emplace_back(&(*it_l1), &(*it_l2));
            else dilepton_pairs.emplace_back(&(*it_l1), &(*it_l2));
          }
        }
      }
      for (auto it_j1=quarks.begin(); it_j1!=quarks.end(); it_j1++){
        for (auto it_j2=it_j1+1; it_j2!=quarks.end(); it_j2++){
          if (it_j1->pdgId()+it_j2->pdgId()==0) diquark_pairs.emplace_back(&(*it_j1), &(*it_j2));
        }
      }
      bool has_4l = false;
      bool has_llqq = false;
      bool has_llnn = false;
      for (auto it_p1=dilepton_pairs.begin(); it_p1!=dilepton_pairs.end(); it_p1++){
        auto const& p1 = *it_p1;
        for (auto it_p2=it_p1+1; it_p2!=dilepton_pairs.end(); it_p2++){
          auto const& p2 = *it_p2;
          if (p1.first==p2.first || p1.first==p2.second || p1.second==p2.first || p1.second==p2.second) continue;
          if (std::abs((p1.first->p4()+p1.second->p4()+p2.first->p4()+p2.second->p4()).M()-125.)<0.05){
            has_4l = true;
            break;
          }
        }
      }
      for (auto const& pl:dilepton_pairs){
        for (auto const& pj:diquark_pairs){
          if (std::abs((pl.first->p4()+pl.second->p4()+pj.first->p4()+pj.second->p4()).M()-125.)<0.05){
            has_llqq = true;
            break;
          }
        }
      }
      for (auto const& pl:dilepton_pairs){
        for (auto const& pn:dineutrino_pairs){
          if (std::abs((pl.first->p4()+pl.second->p4()+pn.first->p4()+pn.second->p4()).M()-125.)<0.05){
            has_llnn = true;
            break;
          }
        }
      }
      has_Hpeak = (has_llqq || has_llnn || has_4l);
    }

    double min_mll = -1;
    double smallestZdiff = -1;
    double mll_smallestZdiff = -1;
    char idx_lepZ11 = -1;
    char idx_lepZ12 = -1;
    char idx_lepZ21 = -1;
    char idx_lepZ22 = -1;
    for (unsigned int il=0; il<leptons.size(); il++){
      auto const& lep1 = leptons.at(il);
      for (unsigned int jl=il+1; jl<leptons.size(); jl++){
        auto const& lep2 = leptons.at(jl);
        if (lep1.pdgId()==-lep2.pdgId()){
          double mll = (lep1.p4() + lep2.p4()).M();
          if (min_mll<0. || mll<min_mll) min_mll = mll;
          double Zdiff = std::abs(mll - PDGHelpers::Zmass);
          if (smallestZdiff<0. || Zdiff<smallestZdiff){
            smallestZdiff = Zdiff;
            mll_smallestZdiff = mll;
            idx_lepZ11 = il;
            idx_lepZ12 = jl;
            if (lep1.pdgId()<lep2.pdgId()) std::swap(idx_lepZ11, idx_lepZ12);
          }
        }
      }
    }

    if (min_mll<4.) continue;

    {
      std::vector<std::pair<unsigned int, unsigned int>> Z2pairs;
      for (unsigned int il=0; il<leptons.size(); il++){
        if (static_cast<char>(il)==idx_lepZ11 || static_cast<char>(il)==idx_lepZ12) continue;
        auto const& lep1 = leptons.at(il);
        for (unsigned int jl=il+1; jl<leptons.size(); jl++){
          if (static_cast<char>(jl)==idx_lepZ11 || static_cast<char>(jl)==idx_lepZ12) continue;
          auto const& lep2 = leptons.at(jl);
          if (lep1.pdgId()==-lep2.pdgId() && std::abs(lep1.pdgId())%2==1){
            double mll = (lep1.p4() + lep2.p4()).M();
            if (mll>12.){
              if (lep1.pdgId()>lep2.pdgId()) Z2pairs.emplace_back(il, jl);
              else Z2pairs.emplace_back(jl, il);
            }
          }
        }
      }
      if (Z2pairs.empty()) continue;
      double highest_ptsum = -1;
      for (auto const& pp:Z2pairs){
        double ptsum = leptons.at(pp.first).pt() + leptons.at(pp.second).pt();
        if (highest_ptsum<ptsum){
          highest_ptsum = ptsum;
          idx_lepZ21 = pp.first;
          idx_lepZ22 = pp.second;
        }
      }
    }

    if (!(idx_lepZ12>=0 && idx_lepZ22>=0 && idx_lepZ11>=0 && idx_lepZ21>=0 && mll_smallestZdiff>=40.)) continue;
    mZ2 = (leptons.at(idx_lepZ21).p4()+leptons.at(idx_lepZ22).p4()).M();
    if (mZ2<12.) continue;
    if (leptons.at(idx_lepZ11).pdgId()==leptons.at(idx_lepZ21).pdgId()){
      double mZ2_alt = std::min(
        (leptons.at(idx_lepZ11).p4()+leptons.at(idx_lepZ22).p4()).M(),
        (leptons.at(idx_lepZ21).p4()+leptons.at(idx_lepZ12).p4()).M()
      );
      if (mZ2_alt<12.) continue;
    }

    std::vector<ParticleObject*> leptons_ZZ; leptons_ZZ.reserve(4);
    leptons_ZZ.push_back(&(leptons.at(idx_lepZ11)));
    leptons_ZZ.push_back(&(leptons.at(idx_lepZ12)));
    leptons_ZZ.push_back(&(leptons.at(idx_lepZ21)));
    leptons_ZZ.push_back(&(leptons.at(idx_lepZ22)));

    double const invmass = (leptons.at(idx_lepZ11).p4()+leptons.at(idx_lepZ12).p4()+leptons.at(idx_lepZ21).p4()+leptons.at(idx_lepZ22).p4()).M();
    if (invmass<220.) continue;

    std::vector<ParticleObject> jets; jets.reserve((*genak4jets_pt)->size());
    for (unsigned int ip=0; ip<(*genak4jets_pt)->size(); ip++){
      ParticleObject::LorentzVector_t pp4;
      pp4 = ParticleObject::PolarLorentzVector_t((*genak4jets_pt)->at(ip), (*genak4jets_eta)->at(ip), (*genak4jets_phi)->at(ip), (*genak4jets_mass)->at(ip));

      double pt = pp4.Pt();
      double eta = pp4.Eta();
      if (pt<30. || std::abs(eta)>=4.7) continue;
      bool pass_dR = true;
      for (auto const& lep:leptons){
        double dRjl = lep.deltaR(pp4);
        if (dRjl<0.4){
          pass_dR = false;
          break;
        }
      }
      if (pass_dR){
        if (pt<50.){
          double phi = pp4.Phi();
          float dummy_SF=1, theEff=1;
          pujetidSFHandler.getSFAndEff(theGlobalSyst, static_cast<float>(pt), static_cast<float>(eta), true, true, true, true, dummy_SF, &theEff);
          TRandom3 rand;
          rand.SetSeed(static_cast<unsigned long long>(std::abs(std::sin(phi)*100000.)) + static_cast<unsigned long long>(ev));
          double rndnum = rand.Uniform();
          if (rndnum>theEff) continue;
        }
        jets.emplace_back(0, pp4);
      }
    }
    njets = jets.size();
    std::vector<ParticleObject*> jet_ptrs; jet_ptrs.reserve(njets);
    for (auto& jet:jets) jet_ptrs.push_back(&jet);
    ParticleObjectHelpers::sortByGreaterPt(jet_ptrs);

    {
      SimpleParticleCollection_t daughters;
      for(auto const& lep:leptons_ZZ) daughters.push_back(SimpleParticle_t(lep->pdgId(), ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(lep->p4())));

      SimpleParticleCollection_t associated;
      for (auto const& jet:jet_ptrs) associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(jet->p4())));

      IvyMELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_ZZ);
      IvyMELAHelpers::melaHandle->setInputEvent(&daughters, &associated, nullptr, false);
      MEblock.computeMELABranches();
      MEblock.pushMELABranches();

      IvyMELAHelpers::melaHandle->resetInputEvent();
    }


    TTree* tout = (!has_Hpeak ? tout_correct : tout_wrong);
    mZZ = invmass;
    id_Z1 = leptons_ZZ.at(0)->pdgId() * leptons_ZZ.at(1)->pdgId();
    id_Z2 = leptons_ZZ.at(2)->pdgId() * leptons_ZZ.at(3)->pdgId();
    mZ1 = (leptons_ZZ.at(0)->p4() + leptons_ZZ.at(1)->p4()).M();
    mZ2 = (leptons_ZZ.at(2)->p4() + leptons_ZZ.at(3)->p4()).M();

    tout->Fill();
  }

  delete tin;

  foutput->WriteTObject(tout_correct); delete tout_correct;
  foutput->WriteTObject(tout_wrong); delete tout_wrong;
  foutput->Close();
}


#include "TemplateHelpers.h"
#include "HistogramKernelDensitySmoothener.h"


void makeTemplates(
  TString strSampleSet, TString period, TString strdate,
  unsigned int icat, unsigned int ichannel
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TDirectory* curdir = gDirectory;
  curdir->cd();

  SampleHelpers::configure(period, Form("store:201221_%s", period.Data()));

  unsigned int idx_channel;
  TString strChannel;
  switch (ichannel){
  case 2:
    strChannel = "4e";
    idx_channel = 121*121;
    break;
  case 1:
    strChannel = "4mu";
    idx_channel = 169*169;
    break;
  default:
    strChannel = "2e2mu";
    idx_channel = 121*169;
    break;
  };

  ACHypothesisHelpers::ProductionType target_prod;
  TString strCategory;
  switch (icat){
  case 2:
    strCategory = "JJVBFTagged";
    target_prod = ACHypothesisHelpers::kVBF;
    break;
  case 1:
    strCategory = "HadVHTagged";
    target_prod = ACHypothesisHelpers::kHadVH;
    break;
  default:
    strCategory = "Untagged";
    target_prod = ACHypothesisHelpers::kGG;
    break;
  };

  std::vector<DiscriminantClasses::Type> KDtypes_obs;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kSM, target_prod, ACHypothesisHelpers::kZZ4l_offshell)){
    if (!HelperFunctions::checkListVariable(KDtypes_obs, KDtype)) KDtypes_obs.push_back(KDtype);
  }
  std::vector<DiscriminantClasses::KDspecs> KDlist_obs; KDlist_obs.reserve(KDtypes_obs.size());
  for (auto const& KDtype:KDtypes_obs) KDlist_obs.emplace_back(KDtype);
  DiscriminantClasses::constructDiscriminants(KDlist_obs, idx_channel, strCategory);
  if (KDlist_obs.size()!=2) IVYerr << "Size of KDlist_obs is " << KDlist_obs.size() << " != 2." << endl;
  std::vector<ExtendedBinning> binning_KDvars; binning_KDvars.reserve(KDlist_obs.size()+1);
  binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("m4l", ACHypothesisHelpers::kSM, target_prod, ACHypothesisHelpers::kZZ4l_offshell));
  for (auto const& KDspec:KDlist_obs){
    binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning(KDspec.KDname, ACHypothesisHelpers::kSM, target_prod, ACHypothesisHelpers::kZZ4l_offshell));
  }
  std::vector<float> smearingStrengthCoeffs(KDlist_obs.size()+1, 1.5);


  std::vector<DiscriminantClasses::Type> KDtypes_cat_VH{
    DiscriminantClasses::kDjjZH, DiscriminantClasses::kDjjWH
  };
  std::vector<DiscriminantClasses::Type> KDtypes_cat_VBF{
    DiscriminantClasses::kDjjVBF
  };
  std::vector<DiscriminantClasses::KDspecs> KDlist_cat_VH; KDlist_cat_VH.reserve(KDtypes_cat_VH.size());
  for (auto const& KDtype:KDtypes_cat_VH) KDlist_cat_VH.emplace_back(KDtype);
  DiscriminantClasses::constructDiscriminants(KDlist_cat_VH, idx_channel, "HadVHTagged");
  std::vector<DiscriminantClasses::KDspecs> KDlist_cat_VBF; KDlist_cat_VBF.reserve(KDtypes_cat_VBF.size());
  for (auto const& KDtype:KDtypes_cat_VBF) KDlist_cat_VBF.emplace_back(KDtype);
  DiscriminantClasses::constructDiscriminants(KDlist_cat_VBF, idx_channel, "JJVBFTagged");

  std::unordered_map<TString, float> ME_values;
  {
    std::vector<TString> MEnames;
    for (auto const& KDspec:KDlist_obs){
      for (auto const& strme:KDspec.KDvars){
        if (!HelperFunctions::checkListVariable(MEnames, strme)){
          MEnames.push_back(strme);
          ME_values[strme]=-1;
        }
      }
    }
    for (auto const& KDspec:KDlist_cat_VH){
      for (auto const& strme:KDspec.KDvars){
        if (!HelperFunctions::checkListVariable(MEnames, strme)){
          MEnames.push_back(strme);
          ME_values[strme]=-1;
        }
      }
    }
    for (auto const& KDspec:KDlist_cat_VBF){
      for (auto const& strme:KDspec.KDvars){
        if (!HelperFunctions::checkListVariable(MEnames, strme)){
          MEnames.push_back(strme);
          ME_values[strme]=-1;
        }
      }
    }
  }


  curdir->cd();

  TFile* foutput = TFile::Open(strSampleSet + "_finalhists_" + strChannel + "_" + strCategory + ".root", "recreate");
  curdir->cd();
  TFile* finput = TFile::Open(strSampleSet + "_corrhists.root", "read");
  curdir->cd();

  std::vector<TString> strHpeaklist{ "Correct", "Wrong" };
  for (auto const& strHpeak:strHpeaklist){
    finput->cd();
    TTree* tin = (TTree*) finput->Get(Form("SkimTree_%s", strHpeak.Data()));

    float wgt_bkg; tin->SetBranchAddress("wgt_bkg", &wgt_bkg);
    float wgt_sig; tin->SetBranchAddress("wgt_sig", &wgt_sig);
    float wgt_bsi; tin->SetBranchAddress("wgt_bsi", &wgt_bsi);
    float mZZ; tin->SetBranchAddress("mZZ", &mZZ);
    float mZ1; tin->SetBranchAddress("mZ1", &mZ1);
    float mZ2; tin->SetBranchAddress("mZ2", &mZ2);
    int id_Z1; tin->SetBranchAddress("id_Z1", &id_Z1);
    int id_Z2; tin->SetBranchAddress("id_Z2", &id_Z2);
    unsigned int nleps_reco; tin->SetBranchAddress("nleps_reco", &nleps_reco);
    unsigned int njets; tin->SetBranchAddress("njets", &njets);
    for (auto& pp:ME_values) tin->SetBranchAddress(pp.first, &pp.second);

    foutput->cd();

    TTree* tout = new TTree(Form("SelectedEvents_%s", strHpeak.Data()), "");
    float KD1, KD2;
    float& wgt = wgt_sig;
    bool selflag = true;
    tout->Branch("wgt", &wgt);
    tout->Branch("mZZ", &mZZ);
    tout->Branch("KD1", &KD1);
    tout->Branch("KD2", &KD2);

    int nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      tin->GetEntry(ev);

      if (mZ1<40. || mZ1>=120. || mZ2>=120.) continue;

      int id_Z1Z2 = id_Z1*id_Z2;
      if (static_cast<int>(idx_channel)!=id_Z1Z2) continue;

      unsigned int icat_sel = 0;
      if (nleps_reco==4 && njets>=2){
        float vmax_KDVH = -1;
        float vKDVBF = -1;
        for (auto const& KDspec:KDlist_cat_VH){
          std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
          for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(ME_values[strKDvar]);
          KDspec.KD->update(KDvars, mZZ);
          float KDval_tmp = *(KDspec.KD);
          vmax_KDVH = std::max(vmax_KDVH, KDval_tmp);
        }
        for (auto const& KDspec:KDlist_cat_VBF){
          std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
          for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(ME_values[strKDvar]);
          KDspec.KD->update(KDvars, mZZ);
          vKDVBF = *(KDspec.KD);
        }

        if (vKDVBF>=0.5) icat_sel = 2;
        else if (vmax_KDVH>=0.5) icat_sel = 1;
      }

      if (icat_sel!=icat) continue;

      {
        unsigned short iKD_obs = 0;
        for (auto const& KDspec:KDlist_obs){
          std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
          for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(ME_values[strKDvar]);
          KDspec.KD->update(KDvars, mZZ);
          if (iKD_obs==0) KD1 = *(KDspec.KD);
          else KD2 = *(KDspec.KD);
          iKD_obs++;
        }
      }
      tout->Fill();
    }

    foutput->WriteTObject(tout);

    using namespace HistogramKernelDensitySmoothener;

    std::vector<TreeHistogramAssociation_3D> tree_hist_assoc; tree_hist_assoc.reserve(1);
    tree_hist_assoc.emplace_back(
      Form("Sig_%s", strHpeak.Data()), Form("Sig_%s", strHpeak.Data()),
      tout, mZZ, KD1, KD2, wgt, selflag
    );
    std::vector<TH3F*> hSmooth = getSimultaneousSmoothHistograms(
      binning_KDvars.at(0), binning_KDvars.at(1), binning_KDvars.at(2),
      tree_hist_assoc,
      smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
      nullptr, nullptr, nullptr
    );

    for (auto const& hh:hSmooth){
      hh->GetXaxis()->SetTitle("m_{4l} (GeV)");
      hh->GetYaxis()->SetTitle(KDlist_obs.at(0).KDlabel);
      hh->GetZaxis()->SetTitle(KDlist_obs.at(1).KDlabel);
      foutput->WriteTObject(hh);
      delete hh;
    }
    delete tout;
  }

  finput->Close();
  foutput->Close();
}

template<typename T> void getCorrRatio(T* hCorrect, T* hWrong, T*& res, TString strappend=""){
  res = (T*) hCorrect->Clone("Sig_CorrRatio" + strappend);
  T* tmp_sum = (T*) hCorrect->Clone("Sig_Sum"); tmp_sum->Add(hWrong);
  res->Divide(tmp_sum);
  delete tmp_sum;
}

void getFullCorrFactors(){
  std::vector<TString> const inputhypos{
    "ZH_ZZ_4LFilter_POWHEG", "WminusH_ZZTo4L_POWHEG", "WplusH_ZZTo4L_POWHEG", "VBF_ZZTo4L_POWHEG"
  };
  std::vector<TString> const channames{
    "2e2mu",
    "4e",
    "4mu"
  };
  std::vector<TString> const catnames{
    "Untagged",
    "HadVHTagged",
    "JJVBFTagged"
  };

  for (auto const& strChannel:channames){
    for (auto const& strCategory:catnames){
      TFile* foutput = TFile::Open(TString("combined_finalhists_") + strChannel + "_" + strCategory + ".root", "recreate");
      std::pair<TH3F*, TH3F*> hout(nullptr, nullptr);

      unsigned short nacc = 0;
      for (auto const& inputhypo:inputhypos){
        TString cinput = inputhypo + "_finalhists_" + strChannel + "_" + strCategory + ".root";
        if (!HostHelpers::FileReadable(cinput)) continue;
        TFile* finput = TFile::Open(cinput, "read");
        finput->cd();
        TH3F* hCorrect = (TH3F*) finput->Get("Sig_Correct"); hCorrect->SetName("hCorrect_tmp");
        TH3F* hWrong = (TH3F*) finput->Get("Sig_Wrong"); hWrong->SetName("hWrong_tmp");
        foutput->cd();
        if (!hout.first){
          hout.first = (TH3F*) hCorrect->Clone("Sig_Correct");
          hout.second = (TH3F*) hWrong->Clone("Sig_Wrong");
        }
        else{
          hout.first->Add(hCorrect);
          hout.second->Add(hWrong);
        }
        finput->Close();
        nacc++;
      }
      if (nacc==0){
        foutput->Close();
        continue;
      }
      foutput->cd();
      TH3F* hRatio = nullptr;
      getCorrRatio(hout.first, hout.second, hRatio, "");
      foutput->WriteTObject(hout.first);
      foutput->WriteTObject(hout.second);
      foutput->WriteTObject(hRatio);

      for (unsigned int ikd=0; ikd<3; ikd++){
        TH1F* hC_proj = HelperFunctions::getHistogramSlice(hout.first, ikd, 1, 99, 1, 99, Form("%s_proj_KD%i", hout.first->GetName(), ikd+1));
        TH1F* hW_proj = HelperFunctions::getHistogramSlice(hout.second, ikd, 1, 99, 1, 99, Form("%s_proj_KD%i", hout.second->GetName(), ikd+1));
        TH1F* hRatio_proj = nullptr;
        getCorrRatio(hC_proj, hW_proj, hRatio_proj, Form("_proj_KD%i", ikd+1));
        hC_proj->GetXaxis()->SetTitle((ikd==0 ? hRatio->GetXaxis() : (ikd==1 ? hRatio->GetYaxis() : hRatio->GetZaxis()))->GetTitle());
        hC_proj->GetYaxis()->SetTitle("Events / bin / lumi.");
        hW_proj->GetXaxis()->SetTitle((ikd==0 ? hRatio->GetXaxis() : (ikd==1 ? hRatio->GetYaxis() : hRatio->GetZaxis()))->GetTitle());
        hW_proj->GetYaxis()->SetTitle("Events / bin / lumi.");
        hRatio_proj->GetXaxis()->SetTitle((ikd==0 ? hRatio->GetXaxis() : (ikd==1 ? hRatio->GetYaxis() : hRatio->GetZaxis()))->GetTitle());
        hRatio_proj->GetYaxis()->SetTitle("off-shell / (off-shell + on-shell)");
        foutput->WriteTObject(hC_proj); delete hC_proj;
        foutput->WriteTObject(hW_proj); delete hW_proj;
        foutput->WriteTObject(hRatio_proj); delete hRatio_proj;
      }

      delete hRatio;
      delete hout.first;
      delete hout.second;
      foutput->Close();
    }
  }
}

void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts="hist",
  bool useLogX=false,
  bool adjustYLow=false,
  float factorYHigh=-1,
  bool addRatioPanel=false
);

void plotFullCorrFactorProjs(){
  TDirectory* curdir = gDirectory;
  curdir->cd();

  std::vector<TString> const catnames{
    "Untagged",
    "HadVHTagged",
    "JJVBFTagged"
  };
  std::vector<TString> const channames{
    "2e2mu",
    "4e",
    "4mu"
  };

  for (auto const& strCategory:catnames){
    std::vector<std::vector<TH1F*>> hKDlist(3, std::vector<TH1F*>(2, nullptr));
    for (auto const& strChannel:channames){
      TFile* finput = TFile::Open(TString("combined_finalhists_") + strChannel + "_" + strCategory + ".root", "read");
      for (unsigned short ikd=0; ikd<3; ikd++){
        finput->cd();
        TH1F* hC = (TH1F*) finput->Get(Form("Sig_Correct_proj_KD%i", ikd+1));
        TH1F* hW = (TH1F*) finput->Get(Form("Sig_Wrong_proj_KD%i", ikd+1));

        hW->Add(hC);

        hC->SetLineColor(kViolet); hC->SetFillColor(kViolet);
        hW->SetLineColor(kGreen+2); hW->SetFillColor(kGreen+2);
        hC->SetLineWidth(2); hC->SetFillStyle(3001);
        hW->SetLineWidth(2); hW->SetFillStyle(3001);

        TString xtitle = hC->GetXaxis()->GetTitle();
        if (xtitle.Contains("[category]")){
          if (strCategory.Contains("VBF")) HelperFunctions::replaceString<TString, TString const>(xtitle, "[category]", "VBF");
          else if (strCategory.Contains("HadVH")) HelperFunctions::replaceString<TString, TString const>(xtitle, "[category]", "VH");
        }
        hC->GetXaxis()->SetTitle(xtitle);
        hW->GetXaxis()->SetTitle(xtitle);
        hC->GetYaxis()->SetTitle("Events / bin");
        hW->GetYaxis()->SetTitle("Events / bin");

        curdir->cd();
        if (!hKDlist.at(ikd).front()){
          hKDlist.at(ikd).front() = (TH1F*) hC->Clone(Form("hC_%s_KD%i", strCategory.Data(), ikd+1));
          hKDlist.at(ikd).back() = (TH1F*) hW->Clone(Form("hT_%s_KD%i", strCategory.Data(), ikd+1));
        }
        else{
          hKDlist.at(ikd).front()->Add(hC);
          hKDlist.at(ikd).back()->Add(hW);
        }
      }
      finput->Close();
    }
    for (unsigned short ikd=0; ikd<3; ikd++){
      makePlot(
        "./",
        1.f,
        Form("c_ZZCheck_%s_KD_%i",strCategory.Data(), ikd+1),
        hKDlist.at(ikd),
        std::vector<TString>{ "EW signal off-shell", "EW signal on-shell" },
        "", "hist",
        (ikd==0),
        false, -1, true
      );
    }
  }
}


void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts,
  bool useLogX,
  bool adjustYLow,
  float factorYHigh,
  bool addRatioPanel
){
  using namespace PlottingHelpers;

  size_t nplottables = hlist.size();
  if (hlabels.size()!=nplottables) return;
  for (auto const& hlabel:hlabels){ if (hlabel=="") nplottables--; }

  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  std::vector<bool> hHasErrors;

  int nbins = -1;
  double ymin = 0;
  if (adjustYLow) ymin=9e9;
  double ymax = -9e9;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    bool hasErrors=false;
    if (nbins<0) nbins = hist->GetNbinsX();
    for (int ix=1; ix<=nbins; ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      if (be!=0.f) hasErrors = true;
      ymax = std::max(ymax, bc+be);
      double bclow=bc; if (be<=bclow) bclow -= be;
      if (adjustYLow && !(bc==0.f && be==0.f)) ymin = std::min(ymin, bclow);
    }
    hHasErrors.push_back(hasErrors);
    //IVYout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : 1.5);
  else ymax /= (factorYHigh>0.f ? factorYHigh : 1.5);
  ymin *= (ymin>=0. ? 0.95 : 1.05);
  for (TH1F* const& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

  TString varlabel;
  TString quantlabel;
  TH1F* hdenom = nullptr;
  TH1F* hratio = nullptr;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();

    if (varlabel=="") varlabel = hist->GetXaxis()->GetTitle();
    if (quantlabel=="") quantlabel = hist->GetYaxis()->GetTitle();

    hist->GetXaxis()->SetTitle("");
    hist->GetYaxis()->SetTitle("");

    if (is==0) hratio = dynamic_cast<TH1F*>(hist->Clone("hratio"));
    else hdenom = dynamic_cast<TH1F*>(hist->Clone("hdenom"));
  }
  for (int ix=1; ix<=hratio->GetNbinsX(); ix++){
    double bt = hdenom->GetBinContent(ix);
    double bc = hratio->GetBinContent(ix);
    hdenom->SetBinContent(ix, 1.);
    hratio->SetBinContent(ix, bc/bt);
  }

  constexpr double npixels_pad_xy = 800;
  CMSLogoStep cmslogotype = kSimulation;
  PlotCanvas plot(canvasname, npixels_pad_xy, npixels_pad_xy, 1, (addRatioPanel ? 2 : 1), 0.2, 0.05, 0.15, 0.07, 0., 0.1, 0.2);
  plot.addCMSLogo(cmslogotype, 13, lumi, 0);

  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);

    hist->GetXaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hist->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hist->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) hist->GetXaxis()->SetLabelSize(0);
  }

  TH1F* hdummy_ratio = nullptr;
  if (addRatioPanel){
    // Iterate in reverse order to preserve the same order of plotting as in the main panel.
    double ymin_ratio = 0;
    double ymax_ratio = 1;

    hdummy_ratio = dynamic_cast<TH1F*>(hdenom->Clone("hdummy_ratio")); hdummy_ratio->Reset("ICESM");
    hdummy_ratio->GetYaxis()->SetRangeUser(ymin_ratio, ymax_ratio);

    hdummy_ratio->GetXaxis()->SetNdivisions(505);
    hdummy_ratio->GetXaxis()->SetLabelFont(43);
    hdummy_ratio->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hdummy_ratio->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hdummy_ratio->GetYaxis()->SetNdivisions(505);
    hdummy_ratio->GetYaxis()->SetLabelFont(43);
    hdummy_ratio->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hdummy_ratio->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
  }

  TPad* pad_hists = plot.getInsidePanels().front().back();
  TPad* pad_ratios = (addRatioPanel ? plot.getInsidePanels().front().front() : nullptr);
  if (useLogX){
    pad_hists->SetLogx();
    pad_ratios->SetLogx();
  }

  // Add x and y titles
  TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
  TLatex* xtitle = new TLatex(); plot.addText(xtitle);
  xtitle->SetTextAlign(22);
  xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  plot.addText(xtitle->DrawLatexNDC(0.5, 0.5, varlabel));

  TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
  TLatex* ytitle = new TLatex(); plot.addText(ytitle);
  ytitle->SetTextAlign(22);
  ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  ytitle->SetTextAngle(90);
  plot.addText(ytitle->DrawLatexNDC(0.5, (addRatioPanel ? 1.-0.5/1.4 : 0.5), quantlabel));
  if (addRatioPanel) plot.addText(ytitle->DrawLatexNDC(0.5, 0.13/1.4, "Comp."));

  pad_hists->cd();

  constexpr double legend_ymax = 0.90;
  double legend_pixelsize = plot.getStdPixelSize_XYTitle();
  double legend_reldy = legend_pixelsize/npixels_pad_xy*1.3;
  TLegend* legend = new TLegend(
    0.50,
    legend_ymax-legend_reldy*float(nplottables),
    0.90,
    legend_ymax,
    "", "NDC"
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  legend->SetTextSize(legend_pixelsize);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  plot.addLegend(legend);
  TText* text;

  pad_hists->cd();

  bool firstHist = true;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString const& hlabel = hlabels.at(is);

    TString stropt = drawopts;
    hist->SetTitle("");
    if (hlabel!=""){
      legend->AddEntry(hist, hlabel, "f");
    }

    if (firstHist){
      hist->Draw(stropt);
      firstHist = false;
    }
    else{
      hist->Draw(stropt+"same");
    }
  }

  pad_hists->cd();

  // Draw in reverse in order to make sure real data is drawn the last.
  for (int is=hlist.size()-1; is>=0; is--){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    TString stropt = drawopts;
    hist->Draw(stropt+"same");
  }

  pad_hists->cd();
  legend->Draw();

  pad_hists->cd();
  TLatex* selectionstitle = new TLatex(); plot.addText(selectionstitle);
  selectionstitle->SetTextAlign(12);
  selectionstitle->SetTextFont(43);
  selectionstitle->SetTextSize(legend_pixelsize);
  {
    double pt_ymax = legend_ymax;
    double pt_dy = legend_reldy;
    for (auto const& strSel:selectionList){
      plot.addText(selectionstitle->DrawLatexNDC(0.25/(1.+0.25+0.0625)+0.05, pt_ymax-pt_dy/2., strSel));
      pt_ymax -= pt_dy;
    }
  }

  if (pad_ratios){
    pad_ratios->cd();
    hdummy_ratio->SetTitle("");
    hdummy_ratio->Draw("hist");
    hdenom->Draw("histsame");
    hratio->Draw("histsame");
    hdenom->Draw("histsame");
    hratio->Draw("histsame");
  }

  plot.update();
  plot.save(coutput_main, "png");
  plot.save(coutput_main, "pdf");

  delete hdummy_ratio;
  delete hratio;
  delete hdenom;
}
