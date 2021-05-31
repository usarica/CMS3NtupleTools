#ifndef PRODUCESIMBKGESTIMATES_ZW3L1NU_H
#define PRODUCESIMBKGESTIMATES_ZW3L1NU_H


// Selection of recorded variables from produceDileptonEvents.cc
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringDn) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringUp) \
  BRANCH_COMMAND(float, event_wgt_PUDn) \
  BRANCH_COMMAND(float, event_wgt_PUUp) \
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
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaUp) \
  BRANCH_COMMAND(float, event_wgt_triggers_Trilepton) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_StatDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_StatUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_SystDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_SystUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_AltMCDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_AltMCUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_StatDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_StatUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_SystDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_SystUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_AltMCDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_AltMCUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons_EffUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId_EffUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging_EffUp) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_fakeableBase) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT) \
  BRANCH_COMMAND(float, event_m3l) \
  BRANCH_COMMAND(float, dPhi_Z_W) \
  BRANCH_COMMAND(float, dPhi_lepW_pTmiss) \
  BRANCH_COMMAND(float, event_mTl) \
  BRANCH_COMMAND(float, event_mWVis) \
  BRANCH_COMMAND(float, dPhi_pTleptonsjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_SF) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(cms3_listSize_t, dilepton_daughter_indices) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(bool, leptons_is_fakeableBase) \
  BRANCH_COMMAND(std::vector<bool>, leptons_pass_fakeable_ids) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(float, leptons_relIso) \
  BRANCH_COMMAND(float, leptons_eff) \
  BRANCH_COMMAND(float, leptons_eff_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_SystUp) \
  BRANCH_COMMAND(float, leptons_eff_DF) \
  BRANCH_COMMAND(float, leptons_eff_DF_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_DF_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_DF_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_DF_SystUp) \
  BRANCH_COMMAND(bool, electrons_conv_vtx_flag) \
  BRANCH_COMMAND(bool, electrons_pass_tightCharge) \
  BRANCH_COMMAND(cms3_electron_missinghits_t, electrons_n_missing_inner_hits) \
  BRANCH_COMMAND(cms3_electron_missinghits_t, electrons_n_all_missing_inner_hits) \
  BRANCH_COMMAND(float, electrons_full5x5_sigmaIEtaIEta) \
  BRANCH_COMMAND(float, electrons_full5x5_sigmaIPhiIPhi) \
  BRANCH_COMMAND(float, electrons_full5x5_r9) \
  BRANCH_COMMAND(float, electrons_seedTime) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched_fullCone) \
  BRANCH_COMMAND(unsigned char, ak4jets_btagWP_Bits) \
  BRANCH_COMMAND(float, ak4jets_NEMF) \
  BRANCH_COMMAND(float, ak4jets_CEMF) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(float, ak8jets_pt) \
  BRANCH_COMMAND(float, ak8jets_eta) \
  BRANCH_COMMAND(float, ak8jets_phi) \
  BRANCH_COMMAND(float, ak8jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


void getMCSampleSet_ZWTo3L1Nu(std::vector< std::pair< TString, std::vector<TString> > >& sampleSpecs){
  OffshellCutflow::FinalStateType const& fstype = OffshellCutflow::activeFinalState;
  if (fstype==OffshellCutflow::fs_ZW_3l1nu){
    switch (SampleHelpers::getDataYear()){
    case 2016:
      sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
        {
          "DY_2l",{ "DY_2l_M_50" }
        },
        {
          "qqZG",{ "ZGJets_ll_nlo_inclusive" }
        },
        {
          "qqZG",{ "ZGJets_ll_nlo_pTG_130-inf" }
        },
        {
          "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
        },
        {
          "qqZZ_4l",{ "qqZZ_4l", "qqZZ_4l_ext" }
        },
        {
          "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
        },
        {
          "qqWZ_3lnu_ext",{ "qqWZ_3lnu_POWHEG_mll_0p1-inf" }
        },
        {
          "TT_2l2nu",{ "TT_2l2nu" }
        },
        {
          "TTW_lnu",{ "TTW_lnu" }
        },
        {
          "TTZ_2l2nu",{ "TTZ_2l2nu_M_10" }
        },
        {
          "TZ_2l_4f",{ "TZ_2l_4f" }
        }
      };
      break;
    case 2017:
      sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
        {
          "DY_2l",{ "DY_2l_M_50", "DY_2l_M_50_ext" }
        },
        {
          "qqZG",{ "ZGJets_ll_pTG_40-130" }
        },
        {
          "qqZG",{ "ZGJets_ll_nlo_pTG_130-inf" }
        },
        {
          "qqZZ_2l2nu",{ "qqZZ_2l2nu_mZ_18-inf" }
        },
        {
          "qqZZ_2l2nu_ext",{ "qqZZ_2l2nu" }
        },
        {
          "qqZZ_4l",{ "qqZZ_4l", "qqZZ_4l_ext" }
        },
        {
          "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG_mll_0p1-inf" }
        },
        {
          "qqWZ_3lnu_ext",{ "qqWZ_3lnu_POWHEG" }
        },
        {
          "TT_2l2nu",{ "TT_2l2nu" }
        },
        {
          "TTW_lnu",{ "TTW_lnu" }
        },
        {
          "TTZ_2l2nu",{ "TTZ_2l2nu" }
        },
        {
          "TZ_2l_4f",{ "TZ_2l_4f" }
        }
      };
      break;
    case 2018:
      sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
        {
          "DY_2l",{ "DY_2l_M_50", "DY_2l_M_50_ext" }
        },
        {
          "qqZG",{ "ZGJets_ll_pTG_40-130" }
        },
        {
          "qqZG",{ "ZGJets_ll_nlo_pTG_130-inf" }
        },
        {
          "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
        },
        {
          "qqZZ_2l2nu_ext",{ "qqZZ_2l2nu_mZ_18-inf" }
        },
        {
          "qqZZ_4l",{ "qqZZ_4l" }
        },
        {
          "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
        },
        {
          "qqWZ_3lnu_ext",{ "qqWZ_3lnu_POWHEG_mll_0p1-inf" }
        },
        {
          "TT_2l2nu",{ "TT_2l2nu" }
        },
        {
          "TTW_lnu",{ "TTW_lnu" }
        },
        {
          "TTZ_2l2nu",{ "TTZ_2l2nu" }
        },
        {
          "TZ_2l_4f",{ "TZ_2l_4f" }
        }
      };
      break;
    }
  }
  else{
    MELAerr << "getMCSampleSet_ZWTo3L1Nu: Final state " << fstype << " is not of the correct type." << endl;
    exit(1);
  }
}


// period: The data period (i.e. "[year]")
// prodVersion: SkimTrees directory version (e.g. "201221_[year]")
// ntupleVersion: Version of trimmed DileptonEvents ntuples, which is separate from the SkimTrees version (e.g. "210107").
// strdate: Tag for the output
void produceSimBkgEstimates_ZWTo3L1Nu(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst, int iCRSF=0,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections

  using namespace OffshellCutflow;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZW_3l1nu);

  SampleHelpers::configure(period, Form("store_skims:%s", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TString cinput_main = "output/3LEvents/SkimTrees/" + ntupleVersion;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Input folder " << cinput_main << " does not exist locally and on the worker directory." << endl;
    exit(1);
  }
  TString const coutput_main = "output/SimBkgEstimates_ZWTo3L1Nu/" + strdate + "/FinalTrees/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  std::vector<TString> transfer_list;

  TriggerScaleFactorHandler triggerSFHandler;

  // Get list of samples
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC;
  getMCSampleDirs(sgroup_sname_sfname_pairs_MC, theGlobalSyst, _JETMETARGS_);
  // Only keep ZG if iCRSF!=0
  if (iCRSF!=0){
    std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC_new;
    sgroup_sname_sfname_pairs_MC_new.reserve(sgroup_sname_sfname_pairs_MC.size());
    for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
      auto const& sgroup = sgroup_sname_sfname_pair.first;
      if (sgroup.Contains("qqZG")) sgroup_sname_sfname_pairs_MC_new.push_back(sgroup_sname_sfname_pair);
    }
    std::swap(sgroup_sname_sfname_pairs_MC_new, sgroup_sname_sfname_pairs_MC);
  }
  std::vector<TString> sgroups;
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    if (!HelperFunctions::checkListVariable(sgroups, sgroup)) sgroups.push_back(sgroup);
  }

  // Get output files and trees
  std::unordered_map<TString, TFile*> sgroup_foutput_map;
  std::unordered_map<TString, BaseTree*> sgroup_tout_map;
  for (auto const& sgroup:sgroups){
    TString stroutput = coutput_main + Form("/finaltree_%s_%s", sgroup.Data(), strSyst.Data());
    if (iCRSF!=0){
      if (sgroup.Contains("qqZG")) HelperFunctions::replaceString<TString, TString const>(stroutput, strSyst, "ZGScaleFactor");
      stroutput += (iCRSF>0 ? "Up" : "Dn");
    }
    stroutput = stroutput + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    foutput->cd();

    BaseTree* tout = new BaseTree("FinalTree");

    tout->putBranch<float>("weight", 1.f);

    tout->putBranch<float>("mTWZ", 0.f);

    tout->putBranch<cms3_id_t>("dilepton_id", 0);
    tout->putBranch<float>("dilepton_mass", 0.f);
    tout->putBranch<float>("dilepton_pt", 0.f);
    tout->putBranch<float>("dilepton_eta", 0.f);

    tout->putBranch<cms3_id_t>("id_Z1", 0);
    tout->putBranch<float>("pt_Z1", 0.f);
    tout->putBranch<float>("eta_Z1", 0.f);
    tout->putBranch<cms3_id_t>("id_Z2", 0);
    tout->putBranch<float>("pt_Z2", 0.f);
    tout->putBranch<float>("eta_Z2", 0.f);
    tout->putBranch<cms3_id_t>("id_lW", 0);
    tout->putBranch<float>("pt_lW", 0.f);
    tout->putBranch<float>("eta_lW", 0.f);
    tout->putBranch<float>("mT_lW", 0.f);

    tout->putBranch<float>("pTmiss", 0.f);
    tout->putBranch<float>("phimiss", 0.f);

    tout->putBranch<unsigned int>("n_ak4jets_pt30", 0); // Number of ak4 jets with pT>=30 GeV
    tout->putBranch<unsigned int>("n_ak4jets_pt30_mass60", 0); // Number of ak4 jets with pT>=30 GeV AND mass>=60 GeV

    // Dijet variables are always calculated for the two leading-pT jets
    tout->putBranch<float>("dijet_mass", -1.f);
    tout->putBranch<float>("dijet_pt", -1.f);
    tout->putBranch<float>("dijet_dEta", -1.f);
    tout->putBranch<float>("dijet_dPhi", -1.f); // Signed dPhi = phi_forward - phi_backward (after re-ordering leading-pT and subleading-pT by pz)

    tout->putBranch<float>("ak4jet_leading_pt", -1.f);
    tout->putBranch<float>("ak4jet_leading_eta", 0.f);
    tout->putBranch<float>("ak4jet_leading_phi", 0.f);
    tout->putBranch<float>("ak4jet_leading_mass", -1.f);

    tout->putBranch<float>("ak4jet_subleading_pt", -1.f);
    tout->putBranch<float>("ak4jet_subleading_eta", 0.f);
    tout->putBranch<float>("ak4jet_subleading_phi", 0.f);
    tout->putBranch<float>("ak4jet_subleading_mass", -1.f);

    // ak8 jet variables, for potential future use
    tout->putBranch<unsigned int>("n_ak8jets_pt200", 0); // Number of ak8 jets with pT>=200 GeV
    tout->putBranch<unsigned int>("n_ak8jets_pt200_mass60to110", 0); // Number of ak8 jets with pT>=200 GeV AND mass within [60, 110) GeV (inclusive/exclusive range)
    tout->putBranch<unsigned int>("n_ak8jets_pt200_mass60to130", 0); // Number of ak8 jets with pT>=200 GeV AND mass within [60, 130) GeV (inclusive/exclusive range)
    tout->putBranch<unsigned int>("n_ak8jets_pt200_mass140", 0); // Number of ak8 jets with pT>=200 GeV AND mass>=140 GeV

    tout->putBranch<float>("ak8jet_leading_pt", -1.f);
    tout->putBranch<float>("ak8jet_leading_eta", 0.f);
    tout->putBranch<float>("ak8jet_leading_mass", -1.f);

    tout->putBranch<std::vector<float>*>("ak8jets_pt200_mass60to130_pt", nullptr);
    tout->putBranch<std::vector<float>*>("ak8jets_pt200_mass60to130_eta", nullptr);
    tout->putBranch<std::vector<float>*>("ak8jets_pt200_mass60to130_mass", nullptr);

    transfer_list.push_back(stroutput);
    sgroup_foutput_map[sgroup] = foutput;
    sgroup_tout_map[sgroup] = tout;
    curdir->cd();
  }

  // Get input trees
  std::unordered_map<TChain*, double> norm_map;
  std::unordered_map<TChain*, double> xsec_scale_map;
  std::vector<std::pair<TString, TChain*>> samples_all;
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    auto const& sname_sfname_pairs = sgroup_sname_sfname_pair.second;
    std::vector<TChain*> tins_collected;
    for (auto const& sname_sfname_pair:sname_sfname_pairs){
      auto const& sname = sname_sfname_pair.first;
      auto const& sfname = sname_sfname_pair.second;

      TString sid = SampleHelpers::getSampleIdentifier(sname);
      float xsec_scale = 1;
      SampleHelpers::hasXSecException(sid, SampleHelpers::getDataYear(), &xsec_scale);
      if (xsec_scale!=1.f) MELAout << "\t- Sample " << sname << " has a cross section exception with scale " << xsec_scale << "." << endl;

      TString cinput = cinput_main + "/" + sfname;
      curdir->cd();
      TChain* tin = new TChain("SkimTree");
      int nfiles = tin->Add(cinput);
      xsec_scale_map[tin] = xsec_scale;
      MELAout << "\t- Successfully added " << nfiles << " files for " << sname << " from " << cinput << "..." << endl;
      samples_all.emplace_back(sgroup, tin);
      tins_collected.push_back(tin);
      norm_map[tin] = 1;
      if (sname_sfname_pairs.size()>1){
        norm_map[tin] = 0;
        std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
        double sum_wgts = 0;
        bool hasCounters = true;
        {
          int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
          int bin_period = 1;
          for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
            if (validDataPeriods.at(iperiod)==SampleHelpers::getDataPeriod()){ bin_period += iperiod+1; break; }
          }
          for (auto const& fname:inputfilenames){
            TFile* ftmp = TFile::Open(fname, "read");
            TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
            if (!hCounters){
              hasCounters = false;
              sum_wgts = 0;
              break;
            }
            sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
            ftmp->Close();
            curdir->cd();
          }
          if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
        }
        norm_map[tin] += sum_wgts;
      }
      curdir->cd();
    }
    {
      double sum_wgts_MC = 0;
      for (auto const& tin:tins_collected) sum_wgts_MC += norm_map[tin];
      for (auto const& tin:tins_collected) norm_map[tin] /= sum_wgts_MC;
    }
  }
  for (auto const& sgroup_tin_pair:samples_all) MELAout
    << "Relative normalization for sample in group " << sgroup_tin_pair.first << " = " << norm_map[sgroup_tin_pair.second]
    << endl;

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float> ME_Kfactor_values;
  for (auto& pp:samples_all){
    auto const& tin = pp.second;

    std::vector<TString> allbranchnames;
    {
      // Then check all leaves
      const TList* llist = (const TList*) tin->GetListOfLeaves();
      if (llist){
        for (int ib=0; ib<llist->GetSize(); ib++){
          auto const& bmem = llist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
         }
      }
      // Then check all branches
      const TList* blist = (const TList*) tin->GetListOfBranches();
      if (blist){
        for (int ib=0; ib<blist->GetSize(); ib++){
          auto const& bmem = blist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
        }
      }

      for (auto const& bname:allbranchnames){
        if (
          (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
          ||
          bname.BeginsWith("KFactor")
          ) ME_Kfactor_values[bname] = -1;
      }
    }

    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    for (auto& it:ME_Kfactor_values){
      TString const& MEname = it.first;
      if (!HelperFunctions::checkListVariable(allbranchnames, MEname)) continue;
      float& MEval = it.second;
      tin->SetBranchStatus(MEname, 1); tin->SetBranchAddress(MEname, &MEval);
    }
  }

  std::vector<double> sfs_ZG_Njets; sfs_ZG_Njets.reserve(3);
  {
    TFile* finput_ZGSFs = TFile::Open(ANALYSISTREEPKGDATAPATH + Form("ScaleFactors/DataDriven/%i/ZGScaleFactors.root", SampleHelpers::getDataYear()), "read");
    TH1D* h_sfs_zg = dynamic_cast<TH1D*>(finput_ZGSFs->Get("sf_zgamma"));
    for (int ix=1; ix<=h_sfs_zg->GetNbinsX(); ix++) sfs_ZG_Njets.push_back(h_sfs_zg->GetBinContent(ix) + (iCRSF==0 ? 0. : (iCRSF>0 ? 1. : -1.))*h_sfs_zg->GetBinError(ix));
    finput_ZGSFs->Close();
  }


  // Keep track of sums of predicted number of events
  std::unordered_map<TString, std::pair<double, double> > sgroup_sumwgts_map;
  std::unordered_map<TString, std::unordered_map<TString, unsigned int> > sgroup_cutlabel_sumevts_map;
  std::vector<TString> strcutlabels{
    "Before cuts",
    "After dPhi_pTlljets_pTmiss",
    "After dPhi_pTll_pTmiss",
    "After mTlW vs pTmiss",
    "After pTmiss",
    "After min_abs_dPhi_pTj_pTmiss",
    "After id_ll",
    "After fakeable lepton veto",
    "After mass window",
    "After b veto",
    "After pT_Z1",
    "After pT_Z2",
    "After pT_lW",
    "After trigger",
    "After gen. match"
  };
  for (auto const& sgroup:sgroups){
    sgroup_sumwgts_map[sgroup] = std::pair<double, double>(0, 0);
    sgroup_cutlabel_sumevts_map[sgroup] = std::unordered_map<TString, unsigned int>();
    for (auto const& strcutlabel:strcutlabels) sgroup_cutlabel_sumevts_map[sgroup][strcutlabel] = 0;
  }

  // Loop over the samples
  for (auto const& spair:samples_all){
    auto const& sgroup = spair.first;
    auto const& tin = spair.second;
    double const& xsec_scale = xsec_scale_map.find(tin)->second;
    double const& norm_scale = norm_map.find(tin)->second;

    TFile* foutput = sgroup_foutput_map[sgroup];
    BaseTree* tout = sgroup_tout_map[sgroup];

    auto& sumwgts_map = sgroup_sumwgts_map[sgroup];
    auto& cutlabel_sumevts_map = sgroup_cutlabel_sumevts_map[sgroup];

    foutput->cd();

    // Reset ME and K factor values
    for (auto& it:ME_Kfactor_values) it.second = -1;
    bool const is_qqVV = sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW");
    bool const is_ggVV = sgroup.Contains("ggZZ") || sgroup.Contains("ggWW") || sgroup.Contains("GGH");
    bool const isData = (sgroup == "Data");

    bool const useNNPDF30 = !isData && sgroup.Contains("TZ_2l");
    bool const requireGenMatchedLeptons = (sgroup.Contains("TTZ_2l2nu") || sgroup.Contains("TZ_2l"));

    float* val_Kfactor_QCD = nullptr;
    float* val_Kfactor_EW = nullptr;
    if (is_qqVV){
      val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_qqVV_Bkg_Nominal")->second);
      switch (theGlobalSyst){
      case tEWDn:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_EWDn")->second);
        break;
      case tEWUp:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_EWUp")->second);
        break;
      default:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_Nominal")->second);
        break;
      }
    }
    if (is_ggVV){
      switch (theGlobalSyst){
      case tQCDScaleDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleDn")->second);
        break;
      case tQCDScaleUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleUp")->second);
        break;
      case tPDFScaleDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleDn")->second);
        break;
      case tPDFScaleUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleUp")->second);
        break;
      case tPDFReplicaDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaDn")->second);
        break;
      case tPDFReplicaUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaUp")->second);
        break;
      case tAsMZDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsDn")->second);
        break;
      case tAsMZUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsUp")->second);
        break;
      default:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second);
        break;
      }
    }

    float* ptr_event_wgt = &event_wgt;
    float* ptr_event_wgt_adjustment = (!useNNPDF30 ? nullptr : &event_wgt_adjustment_NNPDF30);
    float* ptr_event_wgt_syst_adjustment = nullptr;
    float* ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons;
    float* ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons;
    float* ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons;
    float* ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId;
    float* ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging;
    switch (theGlobalSyst){
      case tPDFScaleDn:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PDFScaleDn;
        break;
      case tPDFScaleUp:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PDFScaleUp;
        break;
      case tQCDScaleDn:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_QCDScaleDn;
        break;
      case tQCDScaleUp:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_QCDScaleUp;
        break;
      case tAsMZDn:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_AsMZDn : &event_wgt_adjustment_NNPDF30_AsMZDn);
        break;
      case tAsMZUp:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_AsMZUp : &event_wgt_adjustment_NNPDF30_AsMZUp);
        break;
      case tPDFReplicaDn:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_PDFReplicaDn : &event_wgt_adjustment_NNPDF30_PDFReplicaDn);
        break;
      case tPDFReplicaUp:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_PDFReplicaUp : &event_wgt_adjustment_NNPDF30_PDFReplicaUp);
        break;
      case tPythiaScaleDn:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PythiaScaleDn;
        break;
      case tPythiaScaleUp:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PythiaScaleUp;
        break;

      case eEleEffStatDn:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_StatDn;
        break;
      case eEleEffStatUp:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_StatUp;
        break;
      case eEleEffSystDn:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_SystDn;
        break;
      case eEleEffSystUp:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_SystUp;
        break;
      case eEleEffAltMCDn:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_AltMCDn;
        break;
      case eEleEffAltMCUp:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_AltMCUp;
        break;

      case eMuEffStatDn:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_StatDn;
        break;
      case eMuEffStatUp:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_StatUp;
        break;
      case eMuEffSystDn:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_SystDn;
        break;
      case eMuEffSystUp:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_SystUp;
        break;
      case eMuEffAltMCDn:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_AltMCDn;
        break;
      case eMuEffAltMCUp:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_AltMCUp;
        break;

      case ePhoEffDn:
        ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons_EffDn;
        break;
      case ePhoEffUp:
        ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons_EffUp;
        break;

      case ePUJetIdEffDn:
        ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId_EffDn;
        break;
      case ePUJetIdEffUp:
        ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId_EffUp;
        break;

      case eBTagSFDn:
        ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging_EffDn;
        break;
      case eBTagSFUp:
        ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging_EffUp;
        break;

      case eL1PrefiringDn:
        ptr_event_wgt = &event_wgt_L1PrefiringDn;
        break;
      case eL1PrefiringUp:
        ptr_event_wgt = &event_wgt_L1PrefiringUp;
        break;

      case ePUDn:
        ptr_event_wgt = &event_wgt_PUDn;
        break;
      case ePUUp:
        ptr_event_wgt = &event_wgt_PUUp;
        break;

      default:
        break;
    }

    int const nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      cutlabel_sumevts_map["Before cuts"]++;

      if (!check_dPhi_pTlljets_pTmiss(dPhi_pTleptonsjets_pTmiss)) continue;
      cutlabel_sumevts_map["After dPhi_pTlljets_pTmiss"]++;

      float dPhi_pTboson_pTmiss = TMath::Pi();
      HelperFunctions::deltaPhi(dilepton_phi, event_phimiss, dPhi_pTboson_pTmiss);
      if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      cutlabel_sumevts_map["After dPhi_pTll_pTmiss"]++;

      cms3_listSize_t idx_Z1 = dilepton_daughter_indices->front();
      cms3_listSize_t idx_Z2 = dilepton_daughter_indices->back();
      if (leptons_pt->at(idx_Z1)<leptons_pt->at(idx_Z2)) std::swap(idx_Z1, idx_Z2);
      cms3_listSize_t const idx_lW = 3 - idx_Z1 - idx_Z2;
      cms3_id_t const& idZ1 = leptons_id->at(idx_Z1);
      cms3_id_t const& idZ2 = leptons_id->at(idx_Z2);
      cms3_id_t const& idlW = leptons_id->at(idx_lW);

      if (!check_mTlW_pTmiss(idlW, event_mTl, event_pTmiss)) continue;
      cutlabel_sumevts_map["After mTlW vs pTmiss"]++;
      if (!check_pTmiss(event_pTmiss, event_n_ak4jets_pt30)) continue;
      cutlabel_sumevts_map["After pTmiss"]++;
      if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss, event_n_ak4jets_pt30)) continue;
      cutlabel_sumevts_map["After min_abs_dPhi_pTj_pTmiss"]++;
      if (dilepton_id!=-121 && dilepton_id!=-169) continue;
      cutlabel_sumevts_map["After id_ll"]++;
      if (event_n_leptons_fakeableBase!=0) continue;
      cutlabel_sumevts_map["After fakeable lepton veto"]++;
      // check_mll_QCDsuppression is already applied, so there is no need for it again.
      if (!check_mll(dilepton_mass, true)) continue;
      cutlabel_sumevts_map["After mass window"]++;
      if (!check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)) continue;
      cutlabel_sumevts_map["After b veto"]++;

      float const& pTZ1 = leptons_pt->at(idx_Z1);
      float const& pTZ2 = leptons_pt->at(idx_Z2);
      float const& pTlW = leptons_pt->at(idx_lW);
      // dilepton_daughter_indices is already sorted in pT, so there is no need to check for pT ordering. Just put a failure to the code for sanity.
      if (pTZ1<pTZ2){
        MELAerr << "pTZ1=" << pTZ1 << ", pTZ2=" << pTZ2 << endl;
        exit(1);
      }
      if (!check_pTZ1(pTZ1)) continue;
      cutlabel_sumevts_map["After pT_Z1"]++;
      if (!check_pTZ2(pTZ2)) continue;
      cutlabel_sumevts_map["After pT_Z2"]++;
      if (!check_pTlW(pTlW)) continue;
      cutlabel_sumevts_map["After pT_lW"]++;

      float const& etaZ1 = leptons_eta->at(idx_Z1);
      float const& etaZ2 = leptons_eta->at(idx_Z2);
      float const& etalW = leptons_eta->at(idx_lW);

      if (
        event_wgt_triggers_SingleLepton->at(idx_Z1)!=1.f
        &&
        event_wgt_triggers_SingleLepton->at(idx_Z2)!=1.f
        &&
        event_wgt_triggers_Dilepton_SF->at(idx_Z1+idx_Z2-1)!=1.f
        ) continue;
      cutlabel_sumevts_map["After trigger"]++;

      bool const hasGenMatchedPair = isData || (leptons_is_genMatched_prompt->front() && leptons_is_genMatched_prompt->back() && leptons_is_genMatched_prompt->at(1));
      if (requireGenMatchedLeptons && !hasGenMatchedPair) continue;
      cutlabel_sumevts_map["After gen. match"]++;

      *ptr_event_wgt_SFs_PUJetId = std::min(*ptr_event_wgt_SFs_PUJetId, 3.f);
      double wgt_adjustment = (ptr_event_wgt_adjustment ? *ptr_event_wgt_adjustment : 1.f) * (ptr_event_wgt_syst_adjustment ? *ptr_event_wgt_syst_adjustment : 1.f);
      if (std::abs(wgt_adjustment)>10.f) wgt_adjustment = 10./std::abs(wgt_adjustment);
      float wgt =
        (*ptr_event_wgt) * wgt_adjustment
        * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
        * (val_Kfactor_EW ? *val_Kfactor_EW : 1.f)
        * (*ptr_event_wgt_SFs_muons) * (*ptr_event_wgt_SFs_electrons) * (*ptr_event_wgt_SFs_photons) * (*ptr_event_wgt_SFs_PUJetId) * (*ptr_event_wgt_SFs_btagging)
        * norm_scale * xsec_scale;

      if (!isData){
        float SFself=1, effself=1;
        triggerSFHandler.getCombinedDileptonSFAndEff(
          theGlobalSyst,
          pTZ1, etaZ1, idZ1,
          pTZ2, etaZ2, idZ2,
          true,
          SFself, &effself
        );
        wgt *= SFself;
      }

      // NO MORE MODIFICATION TO wgt BEYOND THIS POINT!
      float const wgtsq = std::pow(wgt, 2);
      sumwgts_map.first += wgt; sumwgts_map.second += wgtsq;

      ROOT::Math::PtEtaPhiMVector p4_dilepton;
      p4_dilepton.SetCoordinates(dilepton_pt, dilepton_eta, dilepton_phi, dilepton_mass);

      ROOT::Math::PtEtaPhiMVector p4_miss;
      p4_miss.SetCoordinates(event_pTmiss, 0.f, event_phimiss, 0.f);

      ROOT::Math::PtEtaPhiMVector p4_lW;
      p4_miss.SetCoordinates(pTlW, etalW, leptons_phi->at(idx_lW), leptons_mass->at(idx_lW));

      ROOT::Math::PtEtaPhiMVector p4_W = p4_miss + p4_lW;

      float out_mTWZ = std::sqrt(
        std::pow(
        (
          std::sqrt(std::pow(dilepton_pt, 2) + std::pow(dilepton_mass, 2))
          + std::sqrt(std::pow(p4_W.Pt(), 2) + std::pow(PDGHelpers::Wmass, 2))
          ), 2
        )
        - std::pow((p4_dilepton + p4_W).Pt(), 2)
      );

      unsigned int out_n_ak4jets_pt30_mass60 = 0;
      ROOT::Math::PtEtaPhiMVector ak4jet_leadingpt, ak4jet_subleadingpt;
      for (unsigned int ijet=0; ijet<ak4jets_mass->size(); ijet++){
        if (ak4jets_mass->at(ijet)>=60.f) out_n_ak4jets_pt30_mass60++;
        if (ijet==0) ak4jet_leadingpt.SetCoordinates(ak4jets_pt->at(ijet), ak4jets_eta->at(ijet), ak4jets_phi->at(ijet), ak4jets_mass->at(ijet));
        else if (ijet==1) ak4jet_subleadingpt.SetCoordinates(ak4jets_pt->at(ijet), ak4jets_eta->at(ijet), ak4jets_phi->at(ijet), ak4jets_mass->at(ijet));
      }
      ROOT::Math::PtEtaPhiMVector p4_dijet = ak4jet_leadingpt + ak4jet_subleadingpt;
      float dijet_dEta, dijet_dPhi;
      HelperFunctions::deltaEta(float(ak4jet_leadingpt.Eta()), float(ak4jet_subleadingpt.Eta()), dijet_dEta);
      if (ak4jet_leadingpt.Pz()>ak4jet_subleadingpt.Pz()){
        HelperFunctions::deltaPhi(float(ak4jet_leadingpt.Phi()), float(ak4jet_subleadingpt.Phi()), dijet_dPhi);
      }
      else{
        HelperFunctions::deltaPhi(float(ak4jet_subleadingpt.Phi()), float(ak4jet_leadingpt.Phi()), dijet_dPhi);
      }

      unsigned int out_n_ak8jets_pt200(0), out_n_ak8jets_pt200_mass60to110(0), out_n_ak8jets_pt200_mass60to130(0), out_n_ak8jets_pt200_mass140(0);
      std::vector<float> out_ak8jets_pt200_mass60to130_pt, out_ak8jets_pt200_mass60to130_eta, out_ak8jets_pt200_mass60to130_mass;
      // Tight ak8 jet selection always ensures pT>=200 GeV, so we only need to look at mass.
      out_n_ak8jets_pt200 = ak8jets_mass->size();
      {
        unsigned int ijet=0;
        for (auto const& ak8jet_mass:(*ak8jets_mass)){
          if (ak8jet_mass>=60.f && ak8jet_mass<110.f){
            out_n_ak8jets_pt200_mass60to110++; out_n_ak8jets_pt200_mass60to130++;
            out_ak8jets_pt200_mass60to130_pt.push_back(ak8jets_pt->at(ijet));
            out_ak8jets_pt200_mass60to130_eta.push_back(ak8jets_eta->at(ijet));
            out_ak8jets_pt200_mass60to130_mass.push_back(ak8jet_mass);
          }
          else if (ak8jet_mass>=60.f && ak8jet_mass<130.f){
            out_n_ak8jets_pt200_mass60to130++;
            out_ak8jets_pt200_mass60to130_pt.push_back(ak8jets_pt->at(ijet));
            out_ak8jets_pt200_mass60to130_eta.push_back(ak8jets_eta->at(ijet));
            out_ak8jets_pt200_mass60to130_mass.push_back(ak8jet_mass);
          }
          else if (ak8jet_mass>=140.f) out_n_ak8jets_pt200_mass140++;
          ijet++;
        }
      }

      // Record the event to the output tree
      tout->setVal<float>("weight", wgt);

      tout->setVal<float>("mTWZ", out_mTWZ);

      tout->setVal<cms3_id_t>("dilepton_id", dilepton_id);
      tout->setVal<float>("dilepton_mass", dilepton_mass);
      tout->setVal<float>("dilepton_pt", dilepton_pt);
      tout->setVal<float>("dilepton_eta", dilepton_eta);

      tout->setVal<cms3_id_t>("id_Z1", idZ1);
      tout->setVal<float>("pt_Z1", pTZ1);
      tout->setVal<float>("eta_Z1", etaZ1);
      tout->setVal<cms3_id_t>("id_Z2", idZ2);
      tout->setVal<float>("pt_Z2", pTZ2);
      tout->setVal<float>("eta_Z2", etaZ2);
      tout->setVal<cms3_id_t>("id_lW", idlW);
      tout->setVal<float>("pt_lW", pTlW);
      tout->setVal<float>("eta_lW", etalW);
      tout->setVal<float>("mT_lW", event_mTl);

      tout->setVal<float>("pTmiss", event_pTmiss);
      tout->setVal<float>("phimiss", event_phimiss);

      tout->setVal<unsigned int>("n_ak4jets_pt30", event_n_ak4jets_pt30);
      tout->setVal<unsigned int>("n_ak4jets_pt30_mass60", out_n_ak4jets_pt30_mass60);
      // No need to set value if njets<2; this is why default values are provided and BaseTree::resetBranches() is called.
      if (event_n_ak4jets_pt30>=2){
        tout->setVal<float>("dijet_mass", p4_dijet.M());
        tout->setVal<float>("dijet_pt", p4_dijet.Pt());
        tout->setVal<float>("dijet_dEta", dijet_dEta);
        tout->setVal<float>("dijet_dPhi", dijet_dPhi);
      }

      if (event_n_ak4jets_pt30>0){
        tout->setVal<float>("ak4jet_leading_pt", ak4jet_leadingpt.Pt());
        tout->setVal<float>("ak4jet_leading_eta", ak4jet_leadingpt.Eta());
        tout->setVal<float>("ak4jet_leading_phi", ak4jet_leadingpt.Phi());
        tout->setVal<float>("ak4jet_leading_mass", ak4jet_leadingpt.M());
        if (event_n_ak4jets_pt30>1){
          tout->setVal<float>("ak4jet_subleading_pt", ak4jet_subleadingpt.Pt());
          tout->setVal<float>("ak4jet_subleading_eta", ak4jet_subleadingpt.Eta());
          tout->setVal<float>("ak4jet_subleading_phi", ak4jet_subleadingpt.Phi());
          tout->setVal<float>("ak4jet_subleading_mass", ak4jet_subleadingpt.M());
        }
      }

      tout->setVal<unsigned int>("n_ak8jets_pt200", out_n_ak8jets_pt200);
      tout->setVal<unsigned int>("n_ak8jets_pt200_mass60to110", out_n_ak8jets_pt200_mass60to110);
      tout->setVal<unsigned int>("n_ak8jets_pt200_mass60to130", out_n_ak8jets_pt200_mass60to130);
      tout->setVal<unsigned int>("n_ak8jets_pt200_mass140", out_n_ak8jets_pt200_mass140);
      tout->setVal("ak8jets_pt200_mass60to130_pt", &out_ak8jets_pt200_mass60to130_pt);
      tout->setVal("ak8jets_pt200_mass60to130_eta", &out_ak8jets_pt200_mass60to130_eta);
      tout->setVal("ak8jets_pt200_mass60to130_mass", &out_ak8jets_pt200_mass60to130_mass);
      if (out_n_ak8jets_pt200>0){
        tout->setVal<float>("ak8jet_leading_pt", ak8jets_pt->front());
        tout->setVal<float>("ak8jet_leading_eta", ak8jets_eta->front());
        tout->setVal<float>("ak8jet_leading_mass", ak8jets_mass->front());
      }

      tout->fill();
      tout->resetBranches();
    }

    curdir->cd();
  }

  for (auto const& sgroup:sgroups){
    MELAout << "Finalizing " << sgroup << ":" << endl;
    TFile* foutput = sgroup_foutput_map[sgroup];
    BaseTree* tout = sgroup_tout_map[sgroup];

    auto const& sumwgts_map = sgroup_sumwgts_map[sgroup];
    MELAout << "\t- Sum of weights: " << sumwgts_map.first << " +- " << std::sqrt(sumwgts_map.second) << endl;

    auto const& cutlabel_sumevts_map = sgroup_cutlabel_sumevts_map.find(sgroup)->second;
    for (auto const& strcutlabel:strcutlabels) MELAout << "\t- " << strcutlabel << ": " << ((double) cutlabel_sumevts_map.find(strcutlabel)->second)/((double) cutlabel_sumevts_map.find(strcutlabels.front())->second) << endl;

    foutput->cd();
    tout->writeToFile(foutput);

    delete tout;
    foutput->Close();
    curdir->cd();
  }

  for (auto& pp:samples_all) delete pp.second;

  for (auto const& fname:transfer_list) SampleHelpers::addToCondorTransferList(fname);
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


#endif
