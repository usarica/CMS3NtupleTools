#include <thread>
#include <chrono>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


constexpr bool useJetOverlapStripping=false;

#define _USE_GENJETS_SYST_ 1
#if _USE_GENJETS_SYST_==1
typedef TH3F hsyst_genjets_t;
typedef TH2F hsyst_genjets_slice_t;
#else
typedef TH2F hsyst_genjets_t;
typedef TH1F hsyst_genjets_slice_t;
#endif


// Dummy for now, but can be extended to veto certain samples for certain systematics
bool checkSystematicForSampleGroup(TString const& sgroup, SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;
  switch (syst){
  case tEWDn:
  case tEWUp:
    return (sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW"));
  default:
    return true;
  }
}

void getMCSampleDirs(
  TString strSampleSet,
  std::vector< std::pair<TString, TString> >& strsamples, SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  bool applyPUIdToAK4Jets, bool applyTightLeptonVetoIdToAK4Jets,
  bool use_MET_Puppi,
  bool use_MET_XYCorr, bool use_MET_JERCorr, bool use_MET_ParticleMomCorr, bool use_MET_p4Preservation, bool use_MET_corrections
){
  using namespace SystematicsHelpers;
  using namespace ACHypothesisHelpers;

  std::vector<SystematicsHelpers::SystematicVariationTypes> const disallowedSysts{
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tPythiaTuneDn, tPythiaTuneUp,
    tHardJetsDn, tHardJetsUp,
    tEWDn, tEWUp,

    eEleEffDn, eEleEffUp,
    eEleEffStatDn, eEleEffStatUp,
    eEleEffSystDn, eEleEffSystUp,
    eEleEffAltMCDn, eEleEffAltMCUp,

    eMuEffDn, eMuEffUp,
    eMuEffStatDn, eMuEffStatUp,
    eMuEffSystDn, eMuEffSystUp,
    eMuEffAltMCDn, eMuEffAltMCUp,

    ePhoEffDn, ePhoEffUp,

    ePUJetIdEffDn, ePUJetIdEffUp,
    eBTagSFDn, eBTagSFUp,

    ePUDn, ePUUp,
    eL1PrefiringDn, eL1PrefiringUp,

    eTriggerEffDn, eTriggerEffUp
  };
  if (HelperFunctions::checkListVariable(disallowedSysts, theGlobalSyst)) theGlobalSyst = sNominal;

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString period = SampleHelpers::getDataPeriod();

  TString cinput_main =
    TString("AK4Jets")
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
  cinput_main = cinput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
  cinput_main = cinput_main
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  cinput_main = cinput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  cinput_main = cinput_main + "/" + period;

  if (!checkSystematicForSampleGroup(strSampleSet, theGlobalSyst)) return;

  std::vector<TString> sdirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sdirs);
  for (auto const& sname:sdirs){
    TString cinput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(cinput, "_MINIAOD", "");
    cinput = cinput + "_" + strSyst + "*.root";
    cinput = cinput_main + "/" + cinput;
    strsamples.emplace_back(sname, cinput);
  }
}


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
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(float, event_wgt_triggers_PFHT_Control) \
  BRANCH_COMMAND(float, event_wgt_triggers_PFMET_MHT_Control) \
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
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(float, event_mZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_fakeableBase) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(float, ak8jets_pt) \
  BRANCH_COMMAND(float, ak8jets_eta) \
  BRANCH_COMMAND(float, ak8jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


using namespace SystematicsHelpers;
using namespace PhysicsProcessHelpers;


PhysicsProcessHandler* getPhysicsProcessHandler(TString strSampleSet, ACHypothesisHelpers::DecayType dktype){
  PhysicsProcessHandler* res = nullptr;
  if (strSampleSet.Contains("GGH") || strSampleSet.Contains("GluGluH")) res = new GGProcessHandler(dktype);
  else if (strSampleSet.Contains("VBF")) res = new VVProcessHandler(dktype, kProcess_VBF);
  else if (strSampleSet.Contains("ZH")) res = new VVProcessHandler(dktype, kProcess_ZH);
  else if (strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH")) res = new VVProcessHandler(dktype, kProcess_WH);
  else{
    MELAerr << "getPhysicsProcessHandler: Cannot identify process " << strSampleSet;
    assert(0);
  }
  return res;
}

void produceReweightingRecords(
  TString strSampleSet,
  TString period, TString prodVersion, TString strdate
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strSampleSet.Contains("/MINIAOD")){
    MELAerr << "Processing single samples is not the design goal of produceReweightingRecords." << endl;
    return;
  }

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  //bool const isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  bool const isPowheg = strSampleSet.Contains("POWHEG");

  ACHypothesisHelpers::DecayType dktype = ACHypothesisHelpers::kZZ2l2nu_offshell;

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampledirs);
  if (sampledirs.empty()){
    MELAerr << "No samples are found for the set " << strSampleSet << "." << endl;
    return;
  }

  // Set output directory
  TString coutput_main = "output/ReweightingRecords/" + strdate + "/" + period;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Declare handlers
  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireGenParticles(false);
  // We will set up the LHE information acquisition options per sample, so no need to do them here.

  /***** REWEIGHTING SETUP *****/
  // Acquire the binning
  ExtendedBinning binning_rewgt;
  binning_rewgt.addBinBoundary(70);
  binning_rewgt.addBinBoundary(13000);
  {
    std::vector<double> sample_masses;
    for (auto const& sname:sampledirs){
      TString sid = SampleHelpers::getSampleIdentifier(sname);
      sample_masses.push_back(SampleHelpers::findPoleMass(sid));
    }
    if (sample_masses.size()>1){
      for (unsigned int im=0; im<sample_masses.size()-1; im++){
        double const& thisMass = sample_masses.at(im);
        double const& nextMass = sample_masses.at(im+1);

        // Special treatment for mH=125 GeV
        // Need to account for the relatively large gap with the next mass point
        if (std::abs(thisMass-125.)<0.8 && std::abs(nextMass-125.)>0.8){
          binning_rewgt.addBinBoundary(thisMass + (nextMass - thisMass)/3.);
          binning_rewgt.addBinBoundary(thisMass + (nextMass - thisMass)*2./3.);
        }
        else binning_rewgt.addBinBoundary((thisMass + nextMass)/2.);
      }
    }
  }
  MELAout << "Reweighting bin boundaries: " << binning_rewgt.getBinningVector() << endl;

  // Acquire the MEs
  std::vector<TString> strMEs;
  {
    PhysicsProcessHandler* proc_handler = getPhysicsProcessHandler(strSampleSet, dktype);
    for (unsigned int iac=0; iac<(unsigned int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis hypo = (ACHypothesisHelpers::ACHypothesis) iac;
      std::vector<TString> strMEs_hypo = proc_handler->getMELAHypothesisWeights(hypo, false);
      HelperFunctions::appendVector(strMEs, strMEs_hypo);
    }
    delete proc_handler;
  }

  // Construct the bulk reweighting builder
  BulkReweightingBuilder rewgtBuilder(
    binning_rewgt,
    { "LHECandMass" },
    { "genHEPMCweight_default" },
    { "xsec" },
    ReweightingFunctions::getSimpleVariableBin,
    ReweightingFunctions::getSimpleWeight,
    ReweightingFunctions::getSimpleWeight
  );

  // Specify the reweighting threshold conditions
  constexpr double tol_wgt = 5;
  float const thr_frac_Neff = (isGG ? 0.005 : 0.01);
  for (auto const& strME:strMEs){
    if (isGG) rewgtBuilder.addReweightingWeights(
      { strME, "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      0.9995, tol_wgt
    );
    else if (strME == "p_Gen_JJEW_SIG_ghv1_1_MCFM") rewgtBuilder.addReweightingWeights(
      { strME, "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      0.9995, tol_wgt
    );
    else rewgtBuilder.addReweightingWeights(
      { strME, "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (isVBF ? 0.999 : 0.9995), tol_wgt
    );
  }

  curdir->cd();

  /***** INPUT SETUP *****/
  // Acquire the samples and register them to the reweighting builder
  // Keep track of the H125 and H200 (or closest higher mass point) samples in order to anchor the nornalization of other samples
  bool allTreesValid = true;
  BaseTree* tree_MH125 = nullptr;
  BaseTree* tree_MHLowestOffshell = nullptr;
  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    TString strinput = SampleHelpers::getDatasetFileName(sname);
    MELAout << "Acquiring " << sname << " from input file(s) " << strinput << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strinput, "cms3ntuple/Events", "", ""); sample_trees.push_back(sample_tree);
    if (!sample_tree->isValid()){
      MELAerr << "\t- Tree is invalid. Aborting..." << endl;
      delete sample_tree;
      for (auto& ss:sample_trees) delete ss;
      allTreesValid = false;
      break;
    }
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    float const sampleMH = SampleHelpers::findPoleMass(sample_tree->sampleIdentifier);
    if (std::abs(sampleMH-125.f)<0.8f) tree_MH125 = sample_tree;
    else if (!tree_MHLowestOffshell && sampleMH>=200.f) tree_MHLowestOffshell = sample_tree;

    bool has_lheMEweights = false;
    bool has_lheparticles = false;
    {
      std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);
      for (auto const& bname:allbranchnames){
        if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights = true;
        else if (bname.Contains(GenInfoHandler::colName_lheparticles)) has_lheparticles = true;
      }
    }

    const int nEntries = sample_tree->getSelectedNEvents();
    bool hasTaus = false;
    double frac_zero_genwgts = 0;
    double sum_wgts_raw_withveto = 0;
    double sum_wgts_raw_withveto_defaultMemberZero = 0;
    double xsec_scale = 1;
    double BR_scale = 1;
    float xsec = 1;

    // Get cross section
    sample_tree->bookBranch<float>("xsec", 0.f);
    sample_tree->getSelectedEvent(0);
    sample_tree->getVal("xsec", xsec);
    xsec *= 1000.;

    // Book branches
    genInfoHandler.setAcquireLHEMEWeights(false);
    // LHE particles are needed to determine the hasTaus flag. This is only needed for the high-mass POWHEG samples.
    genInfoHandler.setAcquireLHEParticles(isPowheg && has_lheparticles);
    genInfoHandler.bookBranches(sample_tree);
    genInfoHandler.wrapTree(sample_tree);

    sample_tree->silenceUnused();

    // Get sums of weights in the sample
    {
      MELAout << "Initiating loop over " << nEntries << " events to determine the sample normalization:" << endl;

      unsigned int n_zero_genwgts=0;
      for (int ev=0; ev<nEntries; ev++){
        HelperFunctions::progressbar(ev, nEntries);
        sample_tree->getEvent(ev);

        genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
        auto const& genInfo = genInfoHandler.getGenInfo();

        auto const& lheparticles = genInfoHandler.getLHEParticles();
        if (!hasTaus){
          for (auto const& part:lheparticles){
            if (part->status()==1 && std::abs(part->pdgId())==15){
              hasTaus = true;
              genInfoHandler.setAcquireLHEParticles(false);
              break;
            }
          }
        }

        double genwgt = genInfo->getGenWeight(true);
        double genwgt_defaultMemberZero = genInfo->extras.LHEweight_defaultMemberZero;
        if (genwgt==0.){
          n_zero_genwgts++;
          continue;
        }

        sum_wgts_raw_withveto_defaultMemberZero += genwgt_defaultMemberZero;
        sum_wgts_raw_withveto += genwgt;
      }
      if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
    }
    xsec_scale = (1. - frac_zero_genwgts);
    // Determine the BR * filter efficiency for the high-mass samples
    if (sampleMH>0.f) BR_scale = SampleHelpers::calculateAdjustedHiggsBREff(sname, sum_wgts_raw_withveto_defaultMemberZero, sum_wgts_raw_withveto, hasTaus);

    // Reset gen. and LHE particle settings, and book those branches as well
    genInfoHandler.setAcquireLHEMEWeights(has_lheMEweights);
    genInfoHandler.setAcquireLHEParticles(false);
    genInfoHandler.bookBranches(sample_tree);

    // Register tree
    MELAout << "\t- Registering the sample for reweighting..." << endl;
    rewgtBuilder.registerTree(sample_tree, xsec_scale * BR_scale / sum_wgts_raw_withveto);

    sample_tree->silenceUnused();
  }
  if (!allTreesValid) return;

  // Setup the pairwise normalization
  std::vector< std::pair<BaseTree*, BaseTree*> > tree_normTree_pairs; tree_normTree_pairs.reserve(sample_trees.size()-1);
  for (unsigned int itree=0; itree<sample_trees.size(); itree++){
    BaseTree* const& sample_tree = sample_trees.at(itree);
    if (sample_tree==tree_MH125 || sample_tree==tree_MHLowestOffshell) continue;
    float const sampleMH = SampleHelpers::findPoleMass(sample_tree->sampleIdentifier);
    if ((tree_MH125 && sampleMH>SampleHelpers::findPoleMass(tree_MH125->sampleIdentifier) && sampleMH<160.f) || (tree_MHLowestOffshell && sampleMH>SampleHelpers::findPoleMass(tree_MHLowestOffshell->sampleIdentifier))){
      tree_normTree_pairs.emplace_back(sample_trees.at(itree), sample_trees.at(itree-1));
      MELAout << "Normalizing mass " << sampleMH << " to mass " << SampleHelpers::findPoleMass(sample_trees.at(itree-1)->sampleIdentifier) << endl;
    }
  }

  /***** REWEIGHTING OUTPUT *****/
  // Determine the output file names
  TString stroutput = Form("%s/%s%s", coutput_main.Data(), strSampleSet.Data(), ".root");
  TString stroutput_txt = stroutput; HelperFunctions::replaceString<TString, TString const>(stroutput_txt, ".root", ".txt");

  // Allow the reweighting record output to be directed to a file so that cross-checks can be made outside this run
  MELAout.open(stroutput_txt.Data());
  // Determine the weight thresholds
  rewgtBuilder.setup(0, &tree_normTree_pairs, thr_frac_Neff);
  rewgtBuilder.print();
  MELAout.close();

  // Record the reweighting information to the output file
  TFile* foutput = TFile::Open(stroutput, "recreate");
  rewgtBuilder.writeToFile(foutput);
  foutput->Close();
  curdir->cd();

  // No need for the inputs anymore
  for (auto& ss:sample_trees) delete ss;

  curdir->cd();

  // Add the output files to condor transfer list
  SampleHelpers::addToCondorTransferList(stroutput);
  SampleHelpers::addToCondorTransferList(stroutput_txt);
}

using namespace SystematicsHelpers;

void produceSystematicsReweighting_MINLO_Pythia(
  TString strSampleSet,
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  bool const isWH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH");
  bool const isZH = strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  bool const usePSWeightsforPythiaScale = (
    (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp)
    &&
    (period=="2017" || period=="2018" || (period=="2016" && (strSampleSet.Contains("WWTo2L2Nu") || strSampleSet.Contains("WW_2LOSFilter"))))
    );
  //unsigned int n_genak4jets_expected_mult = (isWH || isZH || strSampleSet.Contains("2L2Q") || strSampleSet.Contains("LNuQQ"));

  if (
    !(
      theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp
      ||
      (period.Contains("2016") && (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp))
      ||
      (isGG && (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp))
      )
    ) return;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  if (usePSWeightsforPythiaScale) MELAout << "PS weights will be used for Pythia reweighting..." << endl;

  // Acquire the nominal / syst tree pairs
  std::vector<TString> sampledirs_nominal;
  std::vector<TString> sampledirs_syst;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampledirs_nominal);
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs_syst);
  if (sampledirs_syst.empty() || sampledirs_nominal.empty()) return;

  std::vector<float> sample_common_mHvals;
  std::vector< std::pair<TString, TString> > sampledirs_common;
  for (auto const& sname_nominal:sampledirs_nominal){
    TString sid_nominal = SampleHelpers::getSampleIdentifier(sname_nominal);
    float const sampleMH_nominal = SampleHelpers::findPoleMass(sid_nominal);
    bool hasCommon = false;
    TString sname_match;
    for (auto const& sname_syst:sampledirs_syst){
      TString sid_syst = SampleHelpers::getSampleIdentifier(sname_syst);
      float const sampleMH_syst = SampleHelpers::findPoleMass(sid_syst);
      if (sampleMH_nominal == sampleMH_syst){ hasCommon = true; sname_match = sname_syst; break; }
    }
    if (hasCommon){
      sample_common_mHvals.push_back(sampleMH_nominal);
      sampledirs_common.emplace_back(sname_nominal, sname_match);
    }
  }

  // Create the output file
  TString const strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString const coutput_main = "output/SystematicsReweighting/" + strdate + "/" + period;

  // Set output directory
  gSystem->mkdir(coutput_main, true);

  TString stroutput = Form("%s/%s_%s.root", coutput_main.Data(), strSampleSet.Data(), strSyst.Data());
  TFile* foutput = TFile::Open(stroutput, "recreate");

  foutput->cd();

  ExtendedBinning binning_mass({ 0, 150, 600 }, "lheHiggs_mass");
  ExtendedBinning binning_n_genak4jets;
  /*if (n_genak4jets_expected_mult) binning_n_genak4jets = ExtendedBinning({ 0, 2, 4, 6 }, "n_genak4jets");
  else if (isGG && !(theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp)) binning_n_genak4jets = ExtendedBinning({ 0, 1, 2 }, "n_genak4jets");
  else */binning_n_genak4jets = ExtendedBinning({ 0, 1, 2, 3 }, "n_genak4jets"); // This is actually the count without V->qq decays in the different exclusive final state trees.
  ExtendedBinning binning_ptRatio({ 0, 0.4, 0.8, 1.2, 1.6, 2.0 }, "lheHiggs_pt_over_lheHiggs_mass");
  ExtendedBinning binning_genpromptparticles_sump4_pt({ 0, 0.15, 0.3, 0.5, 0.8, 1.2, 1.6, 2.0 }, "genpromptparticles_sump4_pt_over_lheHiggs_mass");
  ExtendedBinning const* zbinning = (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp ? &binning_ptRatio : &binning_genpromptparticles_sump4_pt);

  hsyst_genjets_t h_nominal(
    "h_nominal", "",
    binning_mass.getNbins(), binning_mass.getBinning(),
#if _USE_GENJETS_SYST_==1
    binning_n_genak4jets.getNbins(), binning_n_genak4jets.getBinning(),
    zbinning->getNbins(), zbinning->getBinning()
#else
    zbinning->getNbins(), zbinning->getBinning()
#endif
  );
  hsyst_genjets_t h_syst(
    "h_syst", "",
    binning_mass.getNbins(), binning_mass.getBinning(),
#if _USE_GENJETS_SYST_==1
    binning_n_genak4jets.getNbins(), binning_n_genak4jets.getBinning(),
    zbinning->getNbins(), zbinning->getBinning()
#else
    zbinning->getNbins(), zbinning->getBinning()
#endif
  );
  h_nominal.GetXaxis()->SetTitle(binning_mass.getLabel());
  h_syst.GetXaxis()->SetTitle(binning_mass.getLabel());
#if _USE_GENJETS_SYST_==1
  h_nominal.GetYaxis()->SetTitle(binning_n_genak4jets.getLabel());
  h_syst.GetYaxis()->SetTitle(binning_n_genak4jets.getLabel());
  h_nominal.GetZaxis()->SetTitle(zbinning->getLabel()); h_nominal.Sumw2();
  h_syst.GetZaxis()->SetTitle(zbinning->getLabel()); h_syst.Sumw2();
#else
  h_nominal.GetYaxis()->SetTitle(zbinning->getLabel()); h_nominal.Sumw2();
  h_syst.GetYaxis()->SetTitle(zbinning->getLabel()); h_syst.Sumw2();
#endif

  curdir->cd();

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireGenParticles(true);
  genInfoHandler.setAcquireGenAK4Jets(true);
  // Specify this flag to omit V->qq->jj
  genInfoHandler.setDoGenJetsVDecayCleaning(true);

  for (auto const& sname_pair:sampledirs_common){
    TString const& sname_nominal = sname_pair.first;
    TString const& sname_syst = sname_pair.second;
    TString sid_nominal = SampleHelpers::getSampleIdentifier(sname_nominal);
    TString sid_syst = SampleHelpers::getSampleIdentifier(sname_syst);
    float const sampleMH = SampleHelpers::findPoleMass(sid_nominal);

    MELAout << "Looping over the pair ( " << sname_nominal << ", " << sname_syst << " )..." << endl;

    BaseTree* sample_tree_nominal = new BaseTree(SampleHelpers::getDatasetFileName(sname_nominal), "cms3ntuple/Events", "", "");
    sample_tree_nominal->sampleIdentifier = sid_nominal;
    BaseTree* sample_tree_syst = new BaseTree(SampleHelpers::getDatasetFileName(sname_syst), "cms3ntuple/Events", "", "");
    sample_tree_syst->sampleIdentifier = sid_syst;

    genInfoHandler.bookBranches(sample_tree_nominal);
    genInfoHandler.bookBranches(sample_tree_syst);

    int nEntries = 0;
    double sum_wgts = 0;
    hsyst_genjets_t* htmp = nullptr;

    htmp = (hsyst_genjets_t*) h_nominal.Clone("htmp"); htmp->Reset("ICESM");
    sum_wgts = 0;
    nEntries = sample_tree_nominal->getNEvents();
    genInfoHandler.wrapTree(sample_tree_nominal);
    MELAout << "\t- Looping over " << nEntries << " entries:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree_nominal->getEvent(ev);
      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal);
      auto const& genInfo = genInfoHandler.getGenInfo();
      double genwgt = genInfo->getGenWeight(true);
      double genwgt_adjustment_NNPDF30 = (genwgt==0. ? 1. : genInfo->getGenWeight(false)/genwgt);
      if (std::abs(genwgt_adjustment_NNPDF30)>10.) genwgt_adjustment_NNPDF30 = 10. / std::abs(genwgt_adjustment_NNPDF30);
      genwgt *= genwgt_adjustment_NNPDF30;

      sum_wgts += genwgt;

      auto const& genak4jets = genInfoHandler.getGenAK4Jets();
      unsigned int n_genak4jets = genak4jets.size();

      double lheHiggs_pt=0, lheHiggs_mass=0;
      auto const& lheparticles = genInfoHandler.getLHEParticles();
      for (auto const& part:lheparticles){
        if (PDGHelpers::isAHiggs(part->pdgId())){
          lheHiggs_pt = part->pt();
          lheHiggs_mass = part->mass();
          break;
        }
      }

      ParticleObject::LorentzVector_t genpromptparticles_sump4;
      auto const& genparticles = genInfoHandler.getGenParticles();
      for (auto const& part:genparticles){
        if (part->extras.isPromptFinalState && (PDGHelpers::isAPhoton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()) || PDGHelpers::isALepton(part->pdgId()))){
          genpromptparticles_sump4 += part->p4();
        }
      }

      htmp->Fill(
        lheHiggs_mass,
#if _USE_GENJETS_SYST_==1
        n_genak4jets,
#endif
        (zbinning == &binning_ptRatio ? lheHiggs_pt : static_cast<float>(genpromptparticles_sump4.Pt()))/lheHiggs_mass,
        genwgt
      );
    }
    htmp->Scale(double(nEntries)/sum_wgts);
    HelperFunctions::wipeOverUnderFlows(htmp, false, true);
    h_nominal.Add(htmp, 1.);
    delete htmp;

    htmp = (hsyst_genjets_t*) h_syst.Clone("htmp"); htmp->Reset("ICESM");
    sum_wgts = 0;
    nEntries = sample_tree_syst->getNEvents();
    genInfoHandler.wrapTree(sample_tree_syst);
    MELAout << "\t- Looping over " << nEntries << " entries:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree_syst->getEvent(ev);
      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal);
      auto const& genInfo = genInfoHandler.getGenInfo();
      double genwgt = genInfo->getGenWeight(true);
      double genwgt_adjustment_NNPDF30 = (genwgt==0. ? 1. : genInfo->getGenWeight(false)/genwgt);
      if (usePSWeightsforPythiaScale){
        if (theGlobalSyst==tPythiaScaleDn) genwgt_adjustment_NNPDF30 *= genInfo->extras.PythiaWeight_isr_muR0p25 * genInfo->extras.PythiaWeight_fsr_muR0p25;
        else if (theGlobalSyst==tPythiaScaleUp) genwgt_adjustment_NNPDF30 *= genInfo->extras.PythiaWeight_isr_muR4 * genInfo->extras.PythiaWeight_fsr_muR4;
      }
      if (std::abs(genwgt_adjustment_NNPDF30)>10.) genwgt_adjustment_NNPDF30 = 10. / std::abs(genwgt_adjustment_NNPDF30);
      genwgt *= genwgt_adjustment_NNPDF30;

      sum_wgts += genwgt;

      double genwgt_sf = 1;
      if ((theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp) && sid_syst.Contains("WWTo2L2Nu")){
        auto const& genInfoExtras = genInfo->extras;

        if (theGlobalSyst==tPythiaScaleDn) genwgt_sf = genInfoExtras.PythiaWeight_isr_muR0p25 * genInfoExtras.PythiaWeight_fsr_muR0p25;
        else genwgt_sf = genInfoExtras.PythiaWeight_isr_muR4 * genInfoExtras.PythiaWeight_fsr_muR4;
        if (std::abs(genwgt_sf)>10.) genwgt_sf *= 10./std::abs(genwgt_sf);
      }

      auto const& genak4jets = genInfoHandler.getGenAK4Jets();
      unsigned int n_genak4jets = genak4jets.size();

      double lheHiggs_pt=0, lheHiggs_mass=0;
      auto const& lheparticles = genInfoHandler.getLHEParticles();
      for (auto const& part:lheparticles){
        if (PDGHelpers::isAHiggs(part->pdgId())){
          lheHiggs_pt = part->pt();
          lheHiggs_mass = part->mass();
          break;
        }
      }

      ParticleObject::LorentzVector_t genpromptparticles_sump4;
      auto const& genparticles = genInfoHandler.getGenParticles();
      for (auto const& part:genparticles){
        if (part->extras.isPromptFinalState && (PDGHelpers::isAPhoton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()) || PDGHelpers::isALepton(part->pdgId()))){
          genpromptparticles_sump4 += part->p4();
        }
      }

      htmp->Fill(
        lheHiggs_mass,
#if _USE_GENJETS_SYST_==1
        n_genak4jets,
#endif
        (zbinning == &binning_ptRatio ? lheHiggs_pt : static_cast<float>(genpromptparticles_sump4.Pt()))/lheHiggs_mass,
        genwgt*genwgt_sf
      );
    }
    htmp->Scale(double(nEntries)/sum_wgts);
    HelperFunctions::wipeOverUnderFlows(htmp, false, true);
    h_syst.Add(htmp, 1.);
    delete htmp;

    delete sample_tree_nominal;
    delete sample_tree_syst;
  }

  foutput->WriteTObject(&h_nominal);
  foutput->WriteTObject(&h_syst);

  HelperFunctions::conditionalizeHistogram<hsyst_genjets_t>(&h_nominal, 0, nullptr, false, false);
  HelperFunctions::conditionalizeHistogram<hsyst_genjets_t>(&h_syst, 0, nullptr, false, false);

  hsyst_genjets_t* h_ratio = (hsyst_genjets_t*) h_nominal.Clone("h_ratio");
  if (theGlobalSyst==tHardJetsDn) HelperFunctions::divideHistograms(&h_nominal, &h_syst, h_ratio, false);
  else HelperFunctions::divideHistograms(&h_syst, &h_nominal, h_ratio, false);
  foutput->WriteTObject(h_ratio);

  for (unsigned int isl=0; isl<binning_mass.getNbins(); isl++){
#if _USE_GENJETS_SYST_==1
    hsyst_genjets_slice_t* h_nominal_slice = HelperFunctions::getHistogramSlice(&h_nominal, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_nominal.GetName(), isl));
    hsyst_genjets_slice_t* h_syst_slice = HelperFunctions::getHistogramSlice(&h_syst, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_syst.GetName(), isl));
    hsyst_genjets_slice_t* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));
#else
    hsyst_genjets_slice_t* h_nominal_slice = HelperFunctions::getHistogramSlice(&h_nominal, 1, isl+1, isl+1, Form("%s_slice%u", h_nominal.GetName(), isl));
    hsyst_genjets_slice_t* h_syst_slice = HelperFunctions::getHistogramSlice(&h_syst, 1, isl+1, isl+1, Form("%s_slice%u", h_syst.GetName(), isl));
    hsyst_genjets_slice_t* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));
#endif

    foutput->WriteTObject(h_nominal_slice); delete h_nominal_slice;
    foutput->WriteTObject(h_syst_slice); delete h_syst_slice;
    foutput->WriteTObject(h_ratio_slice); delete h_ratio_slice;
  }

  delete h_ratio;

  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void produceSystematicsReweighting_LHEWeights(
  TString strSampleSet,
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  std::vector<SystematicVariationTypes> allowedSysts;
  if (isGG) allowedSysts = std::vector<SystematicVariationTypes>{
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp
  };
  else if (isVBF) allowedSysts = std::vector<SystematicVariationTypes>{
    tAsMZDn, tAsMZUp
  };
  if (!HelperFunctions::checkListVariable(allowedSysts, theGlobalSyst)) return;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  // Acquire the nominal / syst tree pairs
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampledirs);
  if (sampledirs.empty()) return;

  // Create the output file
  TString const strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString const coutput_main = "output/SystematicsReweighting/" + strdate + "/" + period;

  // Set output directory
  gSystem->mkdir(coutput_main, true);

  TString stroutput = Form("%s/%s_%s.root", coutput_main.Data(), strSampleSet.Data(), strSyst.Data());
  TFile* foutput = TFile::Open(stroutput, "recreate");

  foutput->cd();

  ExtendedBinning binning_mass(18, 100., 1000., "lheHiggs_mass"); binning_mass.addBinBoundary(1250); binning_mass.addBinBoundary(1750); binning_mass.addBinBoundary(2250); binning_mass.addBinBoundary(2750); binning_mass.addBinBoundary(6000);
  //ExtendedBinning binning_lheHiggs_pt(20, 0., 1000., "lheHiggs_pt"); binning_lheHiggs_pt.addBinBoundary(1250); binning_lheHiggs_pt.addBinBoundary(1750); binning_lheHiggs_pt.addBinBoundary(2250); binning_lheHiggs_pt.addBinBoundary(2750); binning_lheHiggs_pt.addBinBoundary(6000);

  TH1F/*TH2F*/ h_nominal_mass_pt(
    "h_nominal", "",
    binning_mass.getNbins(), binning_mass.getBinning()/*,
    binning_lheHiggs_pt.getNbins(), binning_lheHiggs_pt.getBinning()*/
  );
  TH1F/*TH2F*/ h_syst_mass_pt(
    "h_syst", "",
    binning_mass.getNbins(), binning_mass.getBinning()/*,
    binning_lheHiggs_pt.getNbins(), binning_lheHiggs_pt.getBinning()*/
  );
  h_nominal_mass_pt.GetXaxis()->SetTitle(binning_mass.getLabel());
  h_syst_mass_pt.GetXaxis()->SetTitle(binning_mass.getLabel());
  //h_nominal_mass_pt.GetYaxis()->SetTitle(binning_lheHiggs_pt.getLabel());
  //h_syst_mass_pt.GetYaxis()->SetTitle(binning_lheHiggs_pt.getLabel());

  curdir->cd();

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireGenParticles(false);
  genInfoHandler.setAcquireGenAK4Jets(false);

  for (auto const& sname:sampledirs){
    TString sid = SampleHelpers::getSampleIdentifier(sname);
    float const sampleMH = SampleHelpers::findPoleMass(sid);

    MELAout << "Looping over " << sname << "..." << endl;

    BaseTree* sample_tree = new BaseTree(SampleHelpers::getDatasetFileName(sname), "cms3ntuple/Events", "", "");
    sample_tree->sampleIdentifier = sid;

    genInfoHandler.bookBranches(sample_tree);

    TH1F* htmp_nominal_mass_pt = (TH1F*) h_nominal_mass_pt.Clone("htmp_nominal_mass_pt"); htmp_nominal_mass_pt->Reset("ICESM");
    TH1F* htmp_syst_mass_pt = (TH1F*) h_syst_mass_pt.Clone("htmp_syst_mass_pt"); htmp_syst_mass_pt->Reset("ICESM");

    double sum_wgts = 0;
    int nEntries = sample_tree->getNEvents();
    genInfoHandler.wrapTree(sample_tree);
    MELAout << "\t- Looping over " << nEntries << " entries:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree->getEvent(ev);

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal);
      auto const& genInfo = genInfoHandler.getGenInfo();
      auto const& genInfoExtras = genInfo->extras;

      double genwgt = genInfo->getGenWeight(true);
      double genwgt_adjustment_NNPDF30 = (genwgt==0. ? 1. : genInfo->getGenWeight(false)/genwgt);
      if (std::abs(genwgt_adjustment_NNPDF30)>10.) genwgt_adjustment_NNPDF30 = 10. / std::abs(genwgt_adjustment_NNPDF30);
      genwgt *= genwgt_adjustment_NNPDF30;
      sum_wgts += genwgt;

      double lheHiggs_pt=0, lheHiggs_mass=0;
      auto const& lheparticles = genInfoHandler.getLHEParticles();
      for (auto const& part:lheparticles){
        if (PDGHelpers::isAHiggs(part->pdgId())){
          lheHiggs_pt = part->pt();
          lheHiggs_mass = part->mass();
          break;
        }
      }

      double genwgt_syst_sf = 1;
      switch (theGlobalSyst){
      case tQCDScaleDn:
        genwgt_syst_sf = genInfoExtras.LHEweight_QCDscale_muR0p5_muF1;
        break;
      case tQCDScaleUp:
        genwgt_syst_sf = genInfoExtras.LHEweight_QCDscale_muR2_muF1;
        break;
      case tPDFScaleDn:
        genwgt_syst_sf = genInfoExtras.LHEweight_QCDscale_muR1_muF0p5;
        break;
      case tPDFScaleUp:
        genwgt_syst_sf = genInfoExtras.LHEweight_QCDscale_muR1_muF2;
        break;
      case tPDFReplicaDn:
        //genwgt_syst_sf = genInfoExtras.LHEweight_PDFVariation_Dn_default;
        genwgt_syst_sf = genInfoExtras.LHEweight_PDFVariation_Dn_NNPDF30;
        break;
      case tPDFReplicaUp:
        //genwgt_syst_sf = genInfoExtras.LHEweight_PDFVariation_Up_default;
        genwgt_syst_sf = genInfoExtras.LHEweight_PDFVariation_Up_NNPDF30;
        break;
      case tAsMZDn:
        //genwgt_syst_sf = (period=="2016" ? 1.f : genInfoExtras.LHEweight_AsMZ_Dn_default);
        genwgt_syst_sf = (period=="2016" ? 1.f : genInfoExtras.LHEweight_AsMZ_Dn_NNPDF30);
        break;
      case tAsMZUp:
        //genwgt_syst_sf = (period=="2016" ? 1.f : genInfoExtras.LHEweight_AsMZ_Up_default);
        genwgt_syst_sf = (period=="2016" ? 1.f : genInfoExtras.LHEweight_AsMZ_Up_NNPDF30);
        break;
      default:
        break;
      }
      if (std::abs(genwgt_syst_sf)>10.) genwgt_syst_sf = 10. / std::abs(genwgt_syst_sf);

      htmp_nominal_mass_pt->Fill(
        lheHiggs_mass,
        //lheHiggs_pt,
        genwgt
      );
      htmp_syst_mass_pt->Fill(
        lheHiggs_mass,
        //lheHiggs_pt,
        genwgt*genwgt_syst_sf
      );
    }
    htmp_nominal_mass_pt->Scale(double(nEntries)/sum_wgts);
    HelperFunctions::wipeOverUnderFlows(htmp_nominal_mass_pt, false, true);
    h_nominal_mass_pt.Add(htmp_nominal_mass_pt, 1.);
    delete htmp_nominal_mass_pt;
    htmp_syst_mass_pt->Scale(double(nEntries)/sum_wgts);
    HelperFunctions::wipeOverUnderFlows(htmp_syst_mass_pt, false, true);
    h_syst_mass_pt.Add(htmp_syst_mass_pt, 1.);
    delete htmp_syst_mass_pt;

    delete sample_tree;
  }

  foutput->WriteTObject(&h_nominal_mass_pt);
  foutput->WriteTObject(&h_syst_mass_pt);

  if (isVBF){
    TH1F* h_ratio = (TH1F*) h_nominal_mass_pt.Clone("h_ratio"); h_ratio->Reset("ICESM");
    HelperFunctions::divideHistograms(&h_syst_mass_pt, &h_nominal_mass_pt, h_ratio, false);
    foutput->WriteTObject(h_ratio);
    delete h_ratio;
  }
  else if (isGG){
    //TH1F* h_nominal_mass = HelperFunctions::getHistogramSlice(&h_nominal_mass_pt, 0, 1, binning_mass.getNbins(), "h_nominal_mass");
    //TH1F* h_syst_mass = HelperFunctions::getHistogramSlice(&h_syst_mass_pt, 0, 1, binning_mass.getNbins(), "h_syst_mass");

    // In ggH, the ratio along mass goes the other way around to make sure we renormalize the mass shape to actual K factors.
    //TH1F* h_ratio = (TH1F*) h_nominal_mass->Clone("h_ratio"); h_ratio->Reset("ICESM");
    //HelperFunctions::divideHistograms(h_nominal_mass, h_syst_mass, h_ratio, false);
    TH1F* h_ratio = (TH1F*) h_nominal_mass_pt.Clone("h_ratio"); h_ratio->Reset("ICESM");
    HelperFunctions::divideHistograms(&h_nominal_mass_pt, &h_syst_mass_pt, h_ratio, false);
    foutput->WriteTObject(h_ratio);

    /*
    if (theGlobalSyst==tAsMZDn || theGlobalSyst==tAsMZUp){
      // This is the normalized systematics replacement for AsMZ
      TH2F* h_ratio_corr = (TH2F*) h_nominal_mass_pt.Clone("h_ratio_corr"); h_ratio_corr->Reset("ICESM");
      HelperFunctions::divideHistograms(&h_syst_mass_pt, &h_nominal_mass_pt, h_ratio_corr, false);
      for (int ix=1; ix<=h_syst_mass_pt.GetNbinsX(); ix++){
        double r = h_ratio->GetBinContent(ix);
        for (int iy=1; iy<=h_syst_mass_pt.GetNbinsY(); iy++){
          h_ratio_corr->SetBinContent(ix, iy, h_ratio_corr->GetBinContent(ix, iy)*r);
          h_ratio_corr->SetBinError(ix, iy, h_ratio_corr->GetBinError(ix, iy)*r);
        }
      }
      foutput->WriteTObject(h_ratio_corr); delete h_ratio_corr;
    }
    */

    delete h_ratio;
    //delete h_syst_mass;
    //delete h_nominal_mass;
  }

  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void produceSystematicsReweighting(
  TString strSampleSet,
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");

  if (
    theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp
    ||
    (period.Contains("2016") && (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp))
    ||
    (isGG && (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp))
    ) produceSystematicsReweighting_MINLO_Pythia(strSampleSet, period, prodVersion, strdate, theGlobalSyst);
  else produceSystematicsReweighting_LHEWeights(strSampleSet, period, prodVersion, strdate, theGlobalSyst);
}

void combineSystematicsReweightings(
  TString strSampleSet,
  TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  bool const isWH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH");
  bool const isZH = strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");

  std::vector<SystematicVariationTypes> allowedSysts_LHE;
  if (isGG) allowedSysts_LHE = std::vector<SystematicVariationTypes>{
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp
  };
  else if (isVBF) allowedSysts_LHE = std::vector<SystematicVariationTypes>{
    tAsMZDn, tAsMZUp
  };
  bool const isLHESyst = HelperFunctions::checkListVariable(allowedSysts_LHE, theGlobalSyst);
  if (
    !(
      theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp
      ||
      theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp
      ||
      (isGG && (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp))
      ||
      isLHESyst
      )
    ) return;

  TString const strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString strOutPeriod;
  if (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp || (isVBF && isLHESyst)) strOutPeriod = "Combined/2016";
  else if (isGG && (theGlobalSyst==tAsMZDn || theGlobalSyst==tAsMZUp)) strOutPeriod = "Combined/2017_2018";
  else strOutPeriod = "Combined/2016_2017_2018";
  TString const coutput_main = "output/SystematicsReweighting/" + strdate + "/" + strOutPeriod;

  // Set output directory
  gSystem->mkdir(coutput_main, true);

  TString stroutput = Form("%s/%s_%s.root", coutput_main.Data(), strSampleSet.Data(), strSyst.Data());
  TFile* foutput = TFile::Open(stroutput, "recreate");

  std::vector< std::pair<TString, TString> > period_sname_pairs;
  if (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp || theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp){
    std::vector<TString> periods{ "2016", "2017", "2018" };
    for (auto const& period:periods) period_sname_pairs.emplace_back(period, strSampleSet);
    if (strSampleSet == "GGH_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "GGH_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_ZZ2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_ZZ2L2Nu_POWHEG");
    }
  }
  else if (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp){
    std::vector<TString> periods{ "2016" };
    for (auto const& period:periods) period_sname_pairs.emplace_back(period, strSampleSet);
    if (strSampleSet == "GGH_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "GGH_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_ZZ2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_ZZ2L2Nu_POWHEG");
    }
  }
  else if (theGlobalSyst==tAsMZDn || theGlobalSyst==tAsMZUp){
    std::vector<TString> periods{ "2017", "2018" };
    for (auto const& period:periods) period_sname_pairs.emplace_back(period, strSampleSet);
    if (strSampleSet == "GGH_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "GGH_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_ZZ2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_ZZ2L2Nu_POWHEG");
    }
  }
  else if (isLHESyst){
    std::vector<TString> periods{ "2016", "2017", "2018" };
    for (auto const& period:periods) period_sname_pairs.emplace_back(period, strSampleSet);
    if (strSampleSet == "GGH_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "GGH_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "GGH_ZZ2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_ZZ2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_WW2L2Nu_POWHEG");
    }
    else if (strSampleSet == "VBF_WW2L2Nu_POWHEG"){
      for (auto const& period:periods) period_sname_pairs.emplace_back(period, "VBF_ZZ2L2Nu_POWHEG");
    }
  }

  foutput->cd();

  hsyst_genjets_t* h_nominal=nullptr;
  hsyst_genjets_t* h_syst=nullptr;

  hsyst_genjets_slice_t* h_nominal_reduced=nullptr;
  hsyst_genjets_slice_t* h_syst_reduced=nullptr;

  TH1F* h_nominal_1D=nullptr;
  TH1F* h_syst_1D=nullptr;

  for (auto const& period_sname_pair:period_sname_pairs){
    TString const cinput_main = "output/SystematicsReweighting/" + strdate + "/" + period_sname_pair.first;

    TString strinput = Form("%s/%s_%s.root", cinput_main.Data(), period_sname_pair.second.Data(), strSyst.Data());
    if (!SampleHelpers::checkFileOnWorker(strinput)) continue;

    TFile* finput = TFile::Open(strinput, "read");

    finput->cd();

    if (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp || theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp){
      hsyst_genjets_t* htmp_nominal = (hsyst_genjets_t*) finput->Get("h_nominal");
      hsyst_genjets_t* htmp_syst = (hsyst_genjets_t*) finput->Get("h_syst");

      foutput->cd();
      if (!h_nominal) h_nominal = (hsyst_genjets_t*) htmp_nominal->Clone(htmp_nominal->GetName());
      else h_nominal->Add(htmp_nominal, 1.);
      if (!h_syst) h_syst = (hsyst_genjets_t*) htmp_syst->Clone(htmp_syst->GetName());
      else h_syst->Add(htmp_syst, 1.);
    }
    else if (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp){
      hsyst_genjets_t* htmp_nominal = (hsyst_genjets_t*) finput->Get("h_nominal");
      hsyst_genjets_t* htmp_syst = (hsyst_genjets_t*) finput->Get("h_syst");

      foutput->cd();

#if _USE_GENJETS_SYST_==1
      hsyst_genjets_slice_t* htmp_nominal_reduced = HelperFunctions::getHistogramSlice(htmp_nominal, 1, 2, 1, htmp_nominal->GetNbinsX(), Form("%s_reduced", htmp_nominal->GetName()));
      hsyst_genjets_slice_t* htmp_syst_reduced = HelperFunctions::getHistogramSlice(htmp_syst, 1, 2, 1, htmp_syst->GetNbinsX(), Form("%s_reduced", htmp_syst->GetName()));
#else
      hsyst_genjets_slice_t* htmp_nominal_reduced = HelperFunctions::getHistogramSlice(htmp_nominal, 1, 1, htmp_nominal->GetNbinsX(), Form("%s_reduced", htmp_nominal->GetName()));
      hsyst_genjets_slice_t* htmp_syst_reduced = HelperFunctions::getHistogramSlice(htmp_syst, 1, 1, htmp_syst->GetNbinsX(), Form("%s_reduced", htmp_syst->GetName()));
#endif

      if (!h_nominal_reduced) h_nominal_reduced = (hsyst_genjets_slice_t*) htmp_nominal_reduced->Clone(htmp_nominal->GetName());
      else h_nominal_reduced->Add(htmp_nominal_reduced, 1.);
      if (!h_syst_reduced) h_syst_reduced = (hsyst_genjets_slice_t*) htmp_syst_reduced->Clone(htmp_syst->GetName());
      else h_syst_reduced->Add(htmp_syst_reduced, 1.);

      delete htmp_nominal_reduced;
      delete htmp_syst_reduced;
    }
    else if (isLHESyst){
      TH1F* htmp_nominal = (TH1F*) finput->Get("h_nominal");
      TH1F* htmp_syst = (TH1F*) finput->Get("h_syst");

      foutput->cd();
      if (!h_nominal_1D) h_nominal_1D = (TH1F*) htmp_nominal->Clone(htmp_nominal->GetName());
      else h_nominal_1D->Add(htmp_nominal, 1.);
      if (!h_syst_1D) h_syst_1D = (TH1F*) htmp_syst->Clone(htmp_syst->GetName());
      else h_syst_1D->Add(htmp_syst, 1.);
    }

    finput->Close();

    foutput->cd();
  }

  if (h_nominal){
    foutput->WriteTObject(h_nominal);
    foutput->WriteTObject(h_syst);

    HelperFunctions::conditionalizeHistogram<hsyst_genjets_t>(h_nominal, 0, nullptr, false, false);
    HelperFunctions::conditionalizeHistogram<hsyst_genjets_t>(h_syst, 0, nullptr, false, false);

    hsyst_genjets_t* h_ratio = (hsyst_genjets_t*) h_nominal->Clone("h_ratio");
    if (theGlobalSyst==tHardJetsDn) HelperFunctions::divideHistograms(h_nominal, h_syst, h_ratio, false);
    else HelperFunctions::divideHistograms(h_syst, h_nominal, h_ratio, false);
    foutput->WriteTObject(h_ratio);

    for (int isl=0; isl<h_nominal->GetNbinsX(); isl++){
#if _USE_GENJETS_SYST_==1
      hsyst_genjets_slice_t* h_nominal_slice = HelperFunctions::getHistogramSlice(h_nominal, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_nominal->GetName(), isl));
      hsyst_genjets_slice_t* h_syst_slice = HelperFunctions::getHistogramSlice(h_syst, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_syst->GetName(), isl));
      hsyst_genjets_slice_t* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));
#else
      hsyst_genjets_slice_t* h_nominal_slice = HelperFunctions::getHistogramSlice(h_nominal, 1, isl+1, isl+1, Form("%s_slice%u", h_nominal->GetName(), isl));
      hsyst_genjets_slice_t* h_syst_slice = HelperFunctions::getHistogramSlice(h_syst, 1, isl+1, isl+1, Form("%s_slice%u", h_syst->GetName(), isl));
      hsyst_genjets_slice_t* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));
#endif

      foutput->WriteTObject(h_nominal_slice); delete h_nominal_slice;
      foutput->WriteTObject(h_syst_slice); delete h_syst_slice;
      foutput->WriteTObject(h_ratio_slice); delete h_ratio_slice;
    }

    delete h_ratio;
  }
  else if (h_nominal_reduced){
    foutput->WriteTObject(h_nominal_reduced);
    foutput->WriteTObject(h_syst_reduced);

    hsyst_genjets_slice_t* h_ratio = (hsyst_genjets_slice_t*) h_nominal_reduced->Clone("h_ratio");
    HelperFunctions::divideHistograms(h_syst_reduced, h_nominal_reduced, h_ratio, false);
    foutput->WriteTObject(h_ratio);

    delete h_ratio;
  }
  else if (h_nominal_1D){
    foutput->WriteTObject(h_nominal_1D);
    foutput->WriteTObject(h_syst_1D);

    TH1F* h_ratio = (TH1F*) h_nominal_1D->Clone("h_ratio");
    HelperFunctions::divideHistograms(h_syst_1D, h_nominal_1D, h_ratio, false);

    if (isVBF) HelperFunctions::divideHistograms(h_syst_1D, h_nominal_1D, h_ratio, false);
    else if (isGG) HelperFunctions::divideHistograms(h_nominal_1D, h_syst_1D, h_ratio, false);
    foutput->WriteTObject(h_ratio);

    delete h_ratio;
  }

  delete h_syst_1D;
  delete h_nominal_1D;
  delete h_syst_reduced;
  delete h_nominal_reduced;
  delete h_syst;
  delete h_nominal;

  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void runSystematicsReweightingsChain(TString strdate){
  std::vector<SystematicVariationTypes> const allowedSysts{
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaTuneDn, tPythiaTuneUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tHardJetsDn, tHardJetsUp
  };
  std::vector<TString> const strprocesses{
    "GGH_ZZTo2L2Nu_POWHEG",
    "GGH_WWTo2L2Nu_POWHEG",

    "VBF_ZZTo2L2Nu_POWHEG",
    "VBF_WWTo2L2Nu_POWHEG",

    "WminusH_ZZTo2L2Nu_POWHEG",
    //"WminusH_ZZTo2L2Q_POWHEG",
    "WminusH_HToWW_2LOSFilter_POWHEG",

    "WplusH_ZZTo2L2Nu_POWHEG",
    //"WplusH_ZZTo2L2Q_POWHEG",
    "WplusH_HToWW_2LOSFilter_POWHEG",

    "ZH_HTo2Nu2X_2LFilter_POWHEG",
    "ZH_HTo2L2Q_2LFilter_POWHEG",
    "ZH_HTo4Q_2LFilter_POWHEG",
    "ZH_WWTo2L2Nu_POWHEG"/*,
    "ZH_HToLNuQQ_2LFilter_POWHEG"*/
  };
  for (auto const& syst:allowedSysts){ for (auto const& strprocess:strprocesses) combineSystematicsReweightings(strprocess, strdate, syst); }
}


double EvalSystHistogram(TH1F* hist, float const& xvar){
  if (!hist) return 1;

  const int nbinsx = hist->GetNbinsX();
  int ibin = hist->GetXaxis()->FindBin(xvar);
  if (ibin<1) ibin = 1;
  else if (ibin>nbinsx) ibin = nbinsx;

  return hist->GetBinContent(ibin);
}
double EvalSystHistogram(TH2F* hist, float const& xvar, float const& yvar){
  if (!hist) return 1;

  const int nbinsx = hist->GetNbinsX();
  int ibin = hist->GetXaxis()->FindBin(xvar);
  if (ibin<1) ibin = 1;
  else if (ibin>nbinsx) ibin = nbinsx;

  const int nbinsy = hist->GetNbinsY();
  int jbin = hist->GetYaxis()->FindBin(yvar);
  if (jbin<1) jbin = 1;
  else if (jbin>nbinsy) jbin = nbinsy;

  return hist->GetBinContent(ibin, jbin);
}
double EvalSystHistogram(TH3F* hist, float const& xvar, float const& yvar, float const& zvar){
  if (!hist) return 1;

  const int nbinsx = hist->GetNbinsX();
  int ibin = hist->GetXaxis()->FindBin(xvar);
  if (ibin<1) ibin = 1;
  else if (ibin>nbinsx) ibin = nbinsx;

  const int nbinsy = hist->GetNbinsY();
  int jbin = hist->GetYaxis()->FindBin(yvar);
  if (jbin<1) jbin = 1;
  else if (jbin>nbinsy) jbin = nbinsy;

  const int nbinz = hist->GetNbinsZ();
  int kbin = hist->GetZaxis()->FindBin(zvar);
  if (kbin<1) kbin = 1;
  else if (kbin>nbinz) kbin = nbinz;

  return hist->GetBinContent(ibin, jbin, kbin);
}

void getTrees_ZZTo2L2Nu(
  TString strSampleSet,
  TString period, TString prodVersion, TString ntupleVersion, TString rewgtRcdVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
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

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);
  ACHypothesisHelpers::DecayType dktype = ACHypothesisHelpers::kZZ2l2nu_offshell;

  SampleHelpers::configure(period, Form("store_skims:%s", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TString const coutput_main = "output/SignalEstimates_ZZTo2L2Nu/" + strdate + "/FinalTrees/" + period;

  // Acquire the weight record
  TString cinput_rewgtRcd = "output/ReweightingRecords/" + rewgtRcdVersion + "/" + period;
  TString strinput_weights = Form("%s/%s%s", cinput_rewgtRcd.Data(), strSampleSet.Data(), ".root");
  if (!SampleHelpers::checkFileOnWorker(strinput_weights)){
    MELAerr << "Reweighting file " << strinput_weights << " does not exist locally and on the worker directory." << endl;
    assert(0);
  }

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  bool const isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  bool const hasPartialAsMZ = (isGG || isVBF) && SampleHelpers::getDataYear()==2016;

  PhysicsProcessHandler* proc_handler = getPhysicsProcessHandler(strSampleSet, dktype);
  double const proc_norm_scale = proc_handler->getProcessScale();

  std::vector<TString> transfer_list;

  TriggerScaleFactorHandler triggerSFHandler;

  // Get list of samples
  std::vector< std::pair<TString, TString> > sname_sfname_pairs_MC;
  getMCSampleDirs(strSampleSet, sname_sfname_pairs_MC, theGlobalSyst, _JETMETARGS_);

  // Build discriminants
  std::vector<DiscriminantClasses::Type> KDtypes;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kSM, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kL1, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kA2, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kA3, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kL1ZGs, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  // Construct empty KD specs
  std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(KDtypes.size());
  for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);
  // Construct the discriminants
  DiscriminantClasses::constructDiscriminants(KDlist, 0, "JJVBFTagged");

  // Get output file and tree
  TString stroutput = coutput_main + Form("/finaltrees_%s_%s", strSampleSet.Data(), strSyst.Data());
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  // Acquire the systematics histogram as necessary
  TH1F* h_ratio_syst_1D = nullptr;
  TH2F* h_ratio_syst_2D = nullptr;
  TH3F* h_ratio_syst_3D = nullptr;
  {
    std::vector<SystematicVariationTypes> const allowedSysts_1D{
      tPDFScaleDn, tPDFScaleUp,
      tQCDScaleDn, tQCDScaleUp,
      tAsMZDn, tAsMZUp,
      tPDFReplicaDn, tPDFReplicaUp
    };
#if _USE_GENJETS_SYST_==1
    std::vector<SystematicVariationTypes> const allowedSysts_2D{
      tPythiaScaleDn, tPythiaScaleUp
    };
    std::vector<SystematicVariationTypes> const allowedSysts_3D{
#else
    std::vector<SystematicVariationTypes> const allowedSysts_3D;
    std::vector<SystematicVariationTypes> const allowedSysts_2D{
      tPythiaScaleDn, tPythiaScaleUp,
#endif
      tPythiaTuneDn, tPythiaTuneUp,
      tHardJetsDn, tHardJetsUp
    };

    TString strinput_customSyst_main = ANALYSISTREEPKGDATAPATH + "ScaleFactors/SystematicsCustomReweighting/";
    HostHelpers::ExpandEnvironmentVariables(strinput_customSyst_main);
    TString strinput_customSyst;
    if (SampleHelpers::getDataYear()==2016){
      strinput_customSyst = strinput_customSyst_main + "2016/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = strinput_customSyst_main + "2016_2017_2018/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = "";
    }
    else if (SampleHelpers::getDataYear()==2017 || SampleHelpers::getDataYear()==2018){
      strinput_customSyst = strinput_customSyst_main + "2017_2018/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = strinput_customSyst_main + "2016_2017_2018/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = "";
    }
    if (strinput_customSyst!=""){
      MELAout << "Acquiring the systematics file " << strinput_customSyst << ":" << endl;
      TFile* finput_syst = TFile::Open(strinput_customSyst, "read");
      if (HelperFunctions::checkListVariable(allowedSysts_1D, theGlobalSyst)){
        TH1F* htmp = (TH1F*) finput_syst->Get("h_ratio");
        foutput->cd();
        h_ratio_syst_1D = (TH1F*) htmp->Clone(htmp->GetName());
        if (h_ratio_syst_1D) MELAout << "\t- A 1D reweighting histogram " << h_ratio_syst_1D->GetName() << " is found." << endl;
      }
      else if (HelperFunctions::checkListVariable(allowedSysts_2D, theGlobalSyst)){
        TH2F* htmp = (TH2F*) finput_syst->Get("h_ratio");
        foutput->cd();
        h_ratio_syst_2D = (TH2F*) htmp->Clone(htmp->GetName());
        if (h_ratio_syst_2D) MELAout << "\t- A 2D reweighting histogram " << h_ratio_syst_2D->GetName() << " is found." << endl;
      }
      else if (HelperFunctions::checkListVariable(allowedSysts_3D, theGlobalSyst)){
        TH3F* htmp = (TH3F*) finput_syst->Get("h_ratio");
        foutput->cd();
        h_ratio_syst_3D = (TH3F*) htmp->Clone(htmp->GetName());
        if (h_ratio_syst_3D) MELAout << "\t- A 3D reweighting histogram " << h_ratio_syst_3D->GetName() << " is found." << endl;
      }
      finput_syst->Close();
      if (!h_ratio_syst_1D && !h_ratio_syst_2D && !h_ratio_syst_3D){
        MELAerr << "\t- No 1, 2, or 3D reweighting could be acquired. Aborting..." << endl;
        assert(0);
      }
    }
  }

  foutput->cd();

  std::vector<TString> strMEs_allhypos;
  std::unordered_map<TString, TString> strME_strweight_map;
  BaseTree* tout = new BaseTree("FinalTree");
  {
    // Branch weights
    for (unsigned int iac=0; iac<(unsigned int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis hypo = (ACHypothesisHelpers::ACHypothesis) iac;
      TString strHypoName = ACHypothesisHelpers::getACHypothesisName(hypo);

      std::vector<TString> strMEs = proc_handler->getMELAHypothesisWeights(hypo, false);
      HelperFunctions::appendVector(strMEs_allhypos, strMEs);
      std::vector<TString> strOutTreeNames = proc_handler->getOutputTreeNames(hypo, false);
      for (unsigned int itree=0; itree<strMEs.size(); itree++){
        auto const& strME = strMEs.at(itree);
        TString strweight = Form("weight_%s_%s", strHypoName.Data(), strOutTreeNames.at(itree).Data());
        tout->putBranch<float>(strweight, 1.f);
        strME_strweight_map[strME] = strweight;
      }
    }

    tout->putBranch<float>("LHECandMass", -1.f);

    tout->putBranch<float>("mTZZ", 0.f);

    tout->putBranch<cms3_id_t>("dilepton_id", 0);
    tout->putBranch<float>("dilepton_mass", 0.f);
    tout->putBranch<float>("dilepton_pt", 0.f);
    tout->putBranch<float>("dilepton_eta", 0.f);

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

    // Values of various discriminants
    for (auto const& KDspec:KDlist) tout->putBranch<float>(KDspec.KDname, -1.f);
  }

  transfer_list.push_back(stroutput);
  curdir->cd();

  BulkReweightingBuilder rewgtBuilder(
    ExtendedBinning(),
    { "LHECandMass" },
    { "genHEPMCweight_default" },
    { "xsec" },
    ReweightingFunctions::getSimpleVariableBin,
    ReweightingFunctions::getSimpleWeight,
    ReweightingFunctions::getSimpleWeight
  );
  for (auto const& strME:strMEs_allhypos) rewgtBuilder.addReweightingWeights(
    { strME, "p_Gen_CPStoBWPropRewgt" },
    ReweightingFunctions::getSimpleWeight,
    0.9995, 5.
  );

  // Get input trees
  TString const cinput_main = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + ntupleVersion;
  std::unordered_map<BaseTree*, double> norm_map;
  std::unordered_map<BaseTree*, double> xsec_scale_map;
  std::vector< std::pair<TString, BaseTree*> > samples_all;
  for (auto const& sname_sfname_pair:sname_sfname_pairs_MC){
    auto const& sname = sname_sfname_pair.first;
    auto const& sfname = sname_sfname_pair.second;

    TString sid = SampleHelpers::getSampleIdentifier(sname);
    float xsec_scale = 1;
    SampleHelpers::hasXSecException(sid, SampleHelpers::getDataYear(), &xsec_scale);
    if (xsec_scale!=1.f) MELAout << "\t- Sample " << sname << " has a cross section exception with scale " << xsec_scale << "." << endl;

    curdir->cd();

    TString cinput = cinput_main + "/" + sfname;
    BaseTree* tin = new BaseTree(cinput, "SkimTree", "", "");
    tin->sampleIdentifier = sid;
    xsec_scale_map[tin] = xsec_scale;
    MELAout << "\t- Successfully added the input files for " << sname << " from " << cinput << "..." << endl;
    samples_all.emplace_back(sname, tin);
    norm_map[tin] = 1;

    rewgtBuilder.registerTree(tin, 1.);

    curdir->cd();
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE* inptr_##NAME = nullptr;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>** inptr_##NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  // Extra variables for systematics
  std::vector<float>** inptr_genak4jets_pt = nullptr;
  float* inptr_lheHiggs_mass = nullptr;
  float* inptr_lheHiggs_pt = nullptr;
  float* inptr_genpromptparticles_sump4_pt = nullptr;

  std::unordered_map<TString, float*> ME_Kfactor_values;
  for (auto& pp:samples_all){
    auto const& tin = pp.second;

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
      }
    }

#define BRANCH_COMMAND(TYPE, NAME) tin->bookBranch<TYPE>(#NAME, 0);
    BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) tin->bookBranch<std::vector<TYPE>*>(#NAME, nullptr);
    BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

    if (HelperFunctions::checkListVariable<TString>(allbranchnames, "genak4jets_pt")) tin->bookBranch<std::vector<float>*>("genak4jets_pt", nullptr);
    if (HelperFunctions::checkListVariable<TString>(allbranchnames, "lheHiggs_mass")) tin->bookBranch<float>("lheHiggs_mass", -1);
    if (HelperFunctions::checkListVariable<TString>(allbranchnames, "lheHiggs_pt")) tin->bookBranch<float>("lheHiggs_pt", -1);
    if (HelperFunctions::checkListVariable<TString>(allbranchnames, "genpromptparticles_sump4_pt")) tin->bookBranch<float>("genpromptparticles_sump4_pt", -1);

    tin->silenceUnused();
  }

  // Build the reweighting handler from the records
  rewgtBuilder.setupFromFile(strinput_weights);
  rewgtBuilder.print();

  // Loop over the samples
  for (auto const& spair:samples_all){
    auto const& sname = spair.first;
    auto const& tin = spair.second;
    double const& xsec_scale = xsec_scale_map.find(tin)->second;
    double const& norm_scale = norm_map.find(tin)->second;

    MELAout << "Setting up " << tin->sampleIdentifier << "..." << endl;

    constexpr bool useNNPDF30 = true;
    constexpr bool requireGenMatchedLeptons = false;

    MELAout << "\t- Setting up references to standard variables..." << endl;
#define BRANCH_COMMAND(TYPE, NAME) tin->getValRef(#NAME, inptr_##NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    MELAout << "\t- Setting up references to ME and K factor variables..." << endl;
    for (auto& it:ME_Kfactor_values) tin->getValRef(it.first, it.second);

    {
      MELAout << "\t- Setting up references to custom. systematics evaluation variables..." << endl;

      std::vector<TString> allbranchnames_booked;
      tin->getValidBranchNamesWithoutAlias(allbranchnames_booked, true);

      if (HelperFunctions::checkListVariable<TString>(allbranchnames_booked, "genak4jets_pt")) tin->getValRef("genak4jets_pt", inptr_genak4jets_pt);
      if (HelperFunctions::checkListVariable<TString>(allbranchnames_booked, "lheHiggs_mass")) tin->getValRef("lheHiggs_mass", inptr_lheHiggs_mass);
      if (HelperFunctions::checkListVariable<TString>(allbranchnames_booked, "lheHiggs_pt")) tin->getValRef("lheHiggs_pt", inptr_lheHiggs_pt);
      if (HelperFunctions::checkListVariable<TString>(allbranchnames_booked, "genpromptparticles_sump4_pt")) tin->getValRef("genpromptparticles_sump4_pt", inptr_genpromptparticles_sump4_pt);
    }

    foutput->cd();

    MELAout << "\t- Assigning pointers to references for standard variables..." << endl;
#define BRANCH_COMMAND(TYPE, NAME) auto& NAME = *inptr_##NAME;
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND

    // Reset ME and K factor values
    MELAout << "\t- Assigning the K factor pointer..." << endl;
    float* val_Kfactor_QCD = nullptr;
    if (isGG){
      switch (theGlobalSyst){
      case tQCDScaleDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleDn")->second;
        break;
      case tQCDScaleUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleUp")->second;
        break;
      case tPDFScaleDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleDn")->second;
        break;
      case tPDFScaleUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleUp")->second;
        break;
      case tPDFReplicaDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaDn")->second;
        break;
      case tPDFReplicaUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaUp")->second;
        break;
      case tAsMZDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsDn")->second;
        break;
      case tAsMZUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsUp")->second;
        break;
      default:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second;
        break;
      }
    }

    MELAout << "\t- Assigning the CPS reweighting pointer..." << endl;
    float* val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
    if (!val_ME_CPS) MELAerr << "\t\t- Cannot find 'p_Gen_CPStoBWPropRewgt' in ME_Kfactor_values." << endl;

    MELAout << "\t- Assigning the LHECandMass pointer..." << endl;
    float* val_LHECandMass = ME_Kfactor_values.find("LHECandMass")->second;
    if (!val_LHECandMass) MELAerr << "\t\t- Cannot find 'LHECandMass' in ME_Kfactor_values." << endl;


    MELAout << "\t- Assigning pointers to various event weight factors..." << endl;
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
        if (!hasPartialAsMZ) ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_AsMZDn : &event_wgt_adjustment_NNPDF30_AsMZDn);
        break;
      case tAsMZUp:
        if (!hasPartialAsMZ) ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_AsMZUp : &event_wgt_adjustment_NNPDF30_AsMZUp);
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

    int const nEntries = tin->getNEvents();
    MELAout << "\t- Begin looping over " << nEntries << " events:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->getEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (!check_pTmiss(event_pTmiss, event_n_ak4jets_pt30)) continue;
      if (!check_pTboson(dilepton_pt)) continue;
      if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss, event_n_ak4jets_pt30)) continue;
      if (dilepton_id==-143) continue;
      if (event_n_leptons_fakeableBase!=0) continue;
      if (event_wgt_triggers_SingleLepton!=1.f && event_wgt_triggers_Dilepton!=1.f) continue;
      if (!check_mll(dilepton_mass, true)) continue;
      if (!check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)) continue;

      float const pTl1 = std::max(leptons_pt->front(), leptons_pt->back());
      float const pTl2 = std::min(leptons_pt->front(), leptons_pt->back());
      if (!check_pTl1(pTl1)) continue;
      if (!check_pTl2(pTl2)) continue;

      bool const hasGenMatchedPair = (leptons_is_genMatched_prompt->front() && leptons_is_genMatched_prompt->back());
      if (requireGenMatchedLeptons && !hasGenMatchedPair) continue;

      if (!rewgtBuilder.checkWeightsBelowThreshold(tin)) continue;
      double const sample_rewgtnorm_wgt = rewgtBuilder.getOverallReweightingNormalization(tin);

      double syst_corr = 1;
      if (
        inptr_lheHiggs_mass
        &&
        (
          theGlobalSyst==tPDFScaleDn || theGlobalSyst==tPDFScaleUp
          || theGlobalSyst==tQCDScaleDn || theGlobalSyst==tQCDScaleUp
          || theGlobalSyst==tAsMZDn || theGlobalSyst==tAsMZUp
          || theGlobalSyst==tPDFReplicaDn || theGlobalSyst==tPDFReplicaUp
          )
        ) syst_corr = EvalSystHistogram(h_ratio_syst_1D, *inptr_lheHiggs_mass);
      else if (
        inptr_lheHiggs_mass && inptr_genpromptparticles_sump4_pt
        && inptr_genak4jets_pt
        &&
        (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp)
        ) syst_corr = EvalSystHistogram(h_ratio_syst_2D, (*inptr_genak4jets_pt)->size(), *inptr_genpromptparticles_sump4_pt / *inptr_lheHiggs_mass);
      else if (
        inptr_lheHiggs_mass && inptr_genpromptparticles_sump4_pt
#if _USE_GENJETS_SYST_==1
        && inptr_genak4jets_pt
#endif
        &&
        (theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp)
#if _USE_GENJETS_SYST_==1
        ) syst_corr = EvalSystHistogram(h_ratio_syst_3D, *inptr_lheHiggs_mass, (*inptr_genak4jets_pt)->size(), *inptr_genpromptparticles_sump4_pt / *inptr_lheHiggs_mass);
#else
        ) syst_corr = EvalSystHistogram(h_ratio_syst_2D, *inptr_lheHiggs_mass, *inptr_genpromptparticles_sump4_pt / *inptr_lheHiggs_mass);
#endif
      else if (
        inptr_lheHiggs_mass && inptr_lheHiggs_pt
#if _USE_GENJETS_SYST_==1
        && inptr_genak4jets_pt
#endif
        &&
        (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp)
#if _USE_GENJETS_SYST_==1
        ) syst_corr = EvalSystHistogram(h_ratio_syst_3D, *inptr_lheHiggs_mass, (*inptr_genak4jets_pt)->size(), *inptr_lheHiggs_pt / *inptr_lheHiggs_mass);
#else
        ) syst_corr = EvalSystHistogram(h_ratio_syst_2D, *inptr_lheHiggs_mass, *inptr_lheHiggs_pt / *inptr_lheHiggs_mass);
#endif

      *ptr_event_wgt_SFs_PUJetId = std::min(*ptr_event_wgt_SFs_PUJetId, 3.f);
      float wgt =
        (*ptr_event_wgt) * (ptr_event_wgt_adjustment ? *ptr_event_wgt_adjustment : 1.f) * (ptr_event_wgt_syst_adjustment ? *ptr_event_wgt_syst_adjustment : 1.f)
        * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
        * (*ptr_event_wgt_SFs_muons) * (*ptr_event_wgt_SFs_electrons) * (*ptr_event_wgt_SFs_photons) * (*ptr_event_wgt_SFs_PUJetId) * (*ptr_event_wgt_SFs_btagging)
        * norm_scale * xsec_scale * sample_rewgtnorm_wgt * proc_norm_scale
        * syst_corr;

      {
        float SF_trigger=1;
        triggerSFHandler.getCombinedDileptonSFAndEff(
          theGlobalSyst,
          leptons_pt->front(), leptons_eta->front(), leptons_id->front(),
          leptons_pt->back(), leptons_eta->back(), leptons_id->back(),
          true,
          SF_trigger, nullptr
        );
        wgt *= SF_trigger;
      }

      // NO MORE MODIFICATION TO wgt BEYOND THIS POINT!

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
      // Tight ak8 jet selection always ensures pT>=200 GeV, so we only need to look at mass.
      out_n_ak8jets_pt200 = ak8jets_mass->size();
      for (auto const& ak8jet_mass:(*ak8jets_mass)){
        if (ak8jet_mass>=60.f && ak8jet_mass<110.f){ out_n_ak8jets_pt200_mass60to110++; out_n_ak8jets_pt200_mass60to130++; }
        else if (ak8jet_mass>=60.f && ak8jet_mass<130.f) out_n_ak8jets_pt200_mass60to130++;
        else if (ak8jet_mass>=140.f) out_n_ak8jets_pt200_mass140++;
      }

      // Update discriminants
      for (auto& KDspec:KDlist){
        std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
        for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(*(ME_Kfactor_values[strKDvar]));
        KDspec.KD->update(KDvars, event_mZZ); // Use mZZ!
      }

      // Record the event to the output trees
      {
        for (auto const& strME:strMEs_allhypos){
          TString const& strweight = strME_strweight_map.find(strME)->second;

          float* val_ME = ME_Kfactor_values.find(strME)->second;
          float wgt_ME_comb = wgt * (*val_ME) * (*val_ME_CPS);

          tout->setVal<float>(strweight, wgt_ME_comb);
        }

        tout->setVal<float>("LHECandMass", *val_LHECandMass);

        tout->setVal<float>("mTZZ", event_mTZZ);

        tout->setVal<cms3_id_t>("dilepton_id", dilepton_id);
        tout->setVal<float>("dilepton_mass", dilepton_mass);
        tout->setVal<float>("dilepton_pt", dilepton_pt);
        tout->setVal<float>("dilepton_eta", dilepton_eta);

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
        if (out_n_ak8jets_pt200>0){
          tout->setVal<float>("ak8jet_leading_pt", ak8jets_pt->front());
          tout->setVal<float>("ak8jet_leading_eta", ak8jets_eta->front());
          tout->setVal<float>("ak8jet_leading_mass", ak8jets_mass->front());
        }

        if (event_n_ak4jets_pt30>=2){ for (auto const& KDspec:KDlist) tout->setVal<float>(KDspec.KDname, *(KDspec.KD)); }

        tout->fill();
        tout->resetBranches();
      }
    }
    MELAout << "Accumulated " << tout->getNEvents() << " events." << endl;

    curdir->cd();
  }

  // Delete KDs
  for (auto& KDspec:KDlist) KDspec.resetKD();

  // Write the output tree
  foutput->cd();
  BaseTree::setRobustSaveWrite(true);
  tout->writeToFile(foutput);
  BaseTree::setRobustSaveWrite(false);
  delete tout;

  delete h_ratio_syst_3D;
  delete h_ratio_syst_2D;
  delete h_ratio_syst_1D;

  foutput->Close();

  curdir->cd();

  for (auto& pp:samples_all) delete pp.second;

  delete proc_handler;

  for (auto const& fname:transfer_list) SampleHelpers::addToCondorTransferList(fname);
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS

// period: The data period (i.e. "[year]")
// prodVersion: SkimTrees directory version (e.g. "201221_[year]")
// ntupleVersion: Version of trimmed DileptonEvents ntuples, which is separate from the SkimTrees version (e.g. "210107").
// strdate: Tag for the output
void runDistributionsChain(
  TString strSampleSet,
  TString period, TString prodVersion, TString ntupleVersion, TString rewgtRcdVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _VERSIONARGS_ period, prodVersion, ntupleVersion, rewgtRcdVersion, strdate
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections
  getTrees_ZZTo2L2Nu(strSampleSet, _VERSIONARGS_, theGlobalSyst, _JETMETARGS_);
#undef _JETMETARGS_
#undef _VERSIONARGS_
}
