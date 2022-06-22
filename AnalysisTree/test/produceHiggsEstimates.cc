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
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>


constexpr bool useJetOverlapStripping=false;


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


using namespace SystematicsHelpers;
using namespace PhysicsProcessHelpers;


PhysicsProcessHandler* getPhysicsProcessHandler(TString strSampleSet, ACHypothesisHelpers::DecayType dktype){
  PhysicsProcessHandler* res = nullptr;
  if (strSampleSet.Contains("GGH") || strSampleSet.Contains("GluGluH")) res = new GGProcessHandler(dktype);
  else if (strSampleSet.Contains("VBF")) res = new VVProcessHandler(dktype, kProcess_VBF);
  else if (strSampleSet.Contains("ZH")) res = new VVProcessHandler(dktype, kProcess_ZH);
  else if (strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH")) res = new VVProcessHandler(dktype, kProcess_WH);
  else{
    IVYerr << "getPhysicsProcessHandler: Cannot identify process " << strSampleSet;
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
    IVYerr << "Processing single samples is not the design goal of produceReweightingRecords." << endl;
    return;
  }

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  bool const isWH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH");
  //bool const isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  SampleHelpers::HiggsSampleDecayMode const hdecaymode = SampleHelpers::getHiggsSampleDecayMode(strSampleSet);
  bool const hasHWWDecay = SampleHelpers::isHiggsToWWDecay(hdecaymode);
  bool const isPowheg = strSampleSet.Contains("POWHEG");
  bool const hasDirectHWW = isWH || hasHWWDecay;
  bool const isWHWW = isWH && hasDirectHWW;

  ACHypothesisHelpers::DecayType dktype = ACHypothesisHelpers::kZZ2l2nu_offshell;

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampledirs);
  if (sampledirs.empty()){
    IVYerr << "No samples are found for the set " << strSampleSet << "." << endl;
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
  IVYout << "Reweighting bin boundaries: " << binning_rewgt.getBinningVector() << endl;

  // Acquire the MEs
  std::vector<std::pair<TString, bool>> strMEs;
  {
    PhysicsProcessHandler* proc_handler = getPhysicsProcessHandler(strSampleSet, (hdecaymode==SampleHelpers::kZZTo4L ? ACHypothesisHelpers::kZZ4l_offshell : ACHypothesisHelpers::kZZ2l2nu_offshell));
    for (unsigned int iac=0; iac<(unsigned int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis hypo = (ACHypothesisHelpers::ACHypothesis) iac;
      std::vector<TString> strMEs_hypo = proc_handler->getMELAHypothesisWeights(hypo, false);
      std::vector<bool> excludeFromNormRewgt_list(strMEs_hypo.size(), (hasDirectHWW && hypo==ACHypothesisHelpers::kL1ZGs));
      // WW-only interaction can also happen in VBF, but that concerns only one of the L1ZGs ME types.
      if (isVBF && hypo==ACHypothesisHelpers::kL1ZGs){
        for (unsigned int ibsi=0; ibsi<strMEs_hypo.size(); ibsi++){
          if (ibsi==(static_cast<unsigned int>(VVProcessHandler::VVTplSigBSM)-static_cast<unsigned int>(VVProcessHandler::nVVTplSMTypes))) excludeFromNormRewgt_list.at(ibsi)=true;
        }
      }
      std::vector<std::pair<TString, bool>> strMEs_tmp;
      HelperFunctions::zipVectors(strMEs_hypo, excludeFromNormRewgt_list, strMEs_tmp);
      HelperFunctions::appendVector(strMEs, strMEs_tmp);
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
  int iRefHypo = 0;
  constexpr double tol_wgt = 5;
  float const thr_frac_Neff = (isGG ? 0.005 : 0.01);
  for (auto const& strME_pp:strMEs){
    auto const& strME = strME_pp.first;
    auto const& excludeFromNormRewgt = strME_pp.second;
    double thr_wgt = 0.9995;
    if (isGG){ /* Do nothing, 0.9995 is good enough. */ }
    else if (strME == "p_Gen_JJEW_SIG_ghv1_1_MCFM"){ /* Do nothing, 0.9995 is good enough. */ }
    else if (isVBF) thr_wgt = 0.999;
    else if (isWHWW) thr_wgt = 0.995;

    std::vector<TString> wgtcoll{ strME };
    if (isPowheg) wgtcoll.push_back("p_Gen_CPStoBWPropRewgt");
    rewgtBuilder.addReweightingWeights(
      wgtcoll,
      ReweightingFunctions::getSimpleWeight,
      thr_wgt, tol_wgt,
      excludeFromNormRewgt
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
    IVYout << "Acquiring " << sname << " from input file(s) " << strinput << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strinput, "cms3ntuple/Events", "", ""); sample_trees.push_back(sample_tree);
    if (!sample_tree->isValid()){
      IVYerr << "\t- Tree is invalid. Aborting..." << endl;
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
      IVYout << "Initiating loop over " << nEntries << " events to determine the sample normalization:" << endl;

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
    IVYout << "\t- Registering the sample for reweighting..." << endl;
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
      IVYout << "Normalizing mass " << sampleMH << " to mass " << SampleHelpers::findPoleMass(sample_trees.at(itree-1)->sampleIdentifier) << endl;
    }
  }

  /***** REWEIGHTING OUTPUT *****/
  // Determine the output file names
  TString stroutput = Form("%s/%s%s", coutput_main.Data(), strSampleSet.Data(), ".root");
  TString stroutput_txt = stroutput; HelperFunctions::replaceString<TString, TString const>(stroutput_txt, ".root", ".txt");

  // Allow the reweighting record output to be directed to a file so that cross-checks can be made outside this run
  IVYout.open(stroutput_txt.Data());
  // Determine the weight thresholds
  rewgtBuilder.setup(iRefHypo, &tree_normTree_pairs, thr_frac_Neff);
  rewgtBuilder.print();
  IVYout.close();

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

  if (usePSWeightsforPythiaScale) IVYout << "PS weights will be used for Pythia reweighting..." << endl;

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

  TH3F h_nominal(
    "h_nominal", "",
    binning_mass.getNbins(), binning_mass.getBinning(),
    binning_n_genak4jets.getNbins(), binning_n_genak4jets.getBinning(),
    zbinning->getNbins(), zbinning->getBinning()
  );
  TH3F h_syst(
    "h_syst", "",
    binning_mass.getNbins(), binning_mass.getBinning(),
    binning_n_genak4jets.getNbins(), binning_n_genak4jets.getBinning(),
    zbinning->getNbins(), zbinning->getBinning()
  );
  h_nominal.GetXaxis()->SetTitle(binning_mass.getLabel());
  h_syst.GetXaxis()->SetTitle(binning_mass.getLabel());
  h_nominal.GetYaxis()->SetTitle(binning_n_genak4jets.getLabel());
  h_syst.GetYaxis()->SetTitle(binning_n_genak4jets.getLabel());
  h_nominal.GetZaxis()->SetTitle(zbinning->getLabel()); h_nominal.Sumw2();
  h_syst.GetZaxis()->SetTitle(zbinning->getLabel()); h_syst.Sumw2();

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

    IVYout << "Looping over the pair ( " << sname_nominal << ", " << sname_syst << " )..." << endl;

    BaseTree* sample_tree_nominal = new BaseTree(SampleHelpers::getDatasetFileName(sname_nominal), "cms3ntuple/Events", "", "");
    sample_tree_nominal->sampleIdentifier = sid_nominal;
    BaseTree* sample_tree_syst = new BaseTree(SampleHelpers::getDatasetFileName(sname_syst), "cms3ntuple/Events", "", "");
    sample_tree_syst->sampleIdentifier = sid_syst;

    genInfoHandler.bookBranches(sample_tree_nominal);
    genInfoHandler.bookBranches(sample_tree_syst);

    int nEntries = 0;
    double sum_wgts = 0;
    TH3F* htmp = nullptr;

    htmp = (TH3F*) h_nominal.Clone("htmp"); htmp->Reset("ICESM");
    sum_wgts = 0;
    nEntries = sample_tree_nominal->getNEvents();
    genInfoHandler.wrapTree(sample_tree_nominal);
    IVYout << "\t- Looping over " << nEntries << " entries:" << endl;
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
        if (
          part->testSelectionBit(GenParticleSelectionHelpers::kHardPromptFinalVisibleParticle)
          ||
          (part->extras.isPromptFinalState && PDGHelpers::isANeutrino(part->pdgId()))
          ) genpromptparticles_sump4 += part->p4();
      }

      htmp->Fill(
        lheHiggs_mass,
        n_genak4jets,
        (zbinning == &binning_ptRatio ? lheHiggs_pt : static_cast<float>(genpromptparticles_sump4.Pt()))/lheHiggs_mass,
        genwgt
      );
    }
    htmp->Scale(double(nEntries)/sum_wgts);
    HelperFunctions::wipeOverUnderFlows(htmp, false, true);
    h_nominal.Add(htmp, 1.);
    delete htmp;

    htmp = (TH3F*) h_syst.Clone("htmp"); htmp->Reset("ICESM");
    sum_wgts = 0;
    nEntries = sample_tree_syst->getNEvents();
    genInfoHandler.wrapTree(sample_tree_syst);
    IVYout << "\t- Looping over " << nEntries << " entries:" << endl;
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
        if (
          part->testSelectionBit(GenParticleSelectionHelpers::kHardPromptFinalVisibleParticle)
          ||
          (part->extras.isPromptFinalState && PDGHelpers::isANeutrino(part->pdgId()))
          ) genpromptparticles_sump4 += part->p4();
      }

      htmp->Fill(
        lheHiggs_mass,
        n_genak4jets,
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

  HelperFunctions::conditionalizeHistogram<TH3F>(&h_nominal, 0, nullptr, false, false);
  HelperFunctions::conditionalizeHistogram<TH3F>(&h_syst, 0, nullptr, false, false);

  TH3F* h_ratio = (TH3F*) h_nominal.Clone("h_ratio");
  if (theGlobalSyst==tHardJetsDn) HelperFunctions::divideHistograms(&h_nominal, &h_syst, h_ratio, false);
  else HelperFunctions::divideHistograms(&h_syst, &h_nominal, h_ratio, false);
  foutput->WriteTObject(h_ratio);

  for (unsigned int isl=0; isl<binning_mass.getNbins(); isl++){
    TH2F* h_nominal_slice = HelperFunctions::getHistogramSlice(&h_nominal, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_nominal.GetName(), isl));
    TH2F* h_syst_slice = HelperFunctions::getHistogramSlice(&h_syst, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_syst.GetName(), isl));
    TH2F* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));

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

    IVYout << "Looping over " << sname << "..." << endl;

    BaseTree* sample_tree = new BaseTree(SampleHelpers::getDatasetFileName(sname), "cms3ntuple/Events", "", "");
    sample_tree->sampleIdentifier = sid;

    genInfoHandler.bookBranches(sample_tree);

    TH1F* htmp_nominal_mass_pt = (TH1F*) h_nominal_mass_pt.Clone("htmp_nominal_mass_pt"); htmp_nominal_mass_pt->Reset("ICESM");
    TH1F* htmp_syst_mass_pt = (TH1F*) h_syst_mass_pt.Clone("htmp_syst_mass_pt"); htmp_syst_mass_pt->Reset("ICESM");

    double sum_wgts = 0;
    int nEntries = sample_tree->getNEvents();
    genInfoHandler.wrapTree(sample_tree);
    IVYout << "\t- Looping over " << nEntries << " entries:" << endl;
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

  TH3F* h_nominal=nullptr;
  TH3F* h_syst=nullptr;

  TH2F* h_nominal_reduced=nullptr;
  TH2F* h_syst_reduced=nullptr;

  TH1F* h_nominal_1D=nullptr;
  TH1F* h_syst_1D=nullptr;

  for (auto const& period_sname_pair:period_sname_pairs){
    TString const cinput_main = "output/SystematicsReweighting/" + strdate + "/" + period_sname_pair.first;

    TString strinput = Form("%s/%s_%s.root", cinput_main.Data(), period_sname_pair.second.Data(), strSyst.Data());
    if (!SampleHelpers::checkFileOnWorker(strinput)) continue;

    TFile* finput = TFile::Open(strinput, "read");

    finput->cd();

    if (theGlobalSyst==tHardJetsDn || theGlobalSyst==tHardJetsUp || theGlobalSyst==tPythiaTuneDn || theGlobalSyst==tPythiaTuneUp || theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp){
      TH3F* htmp_nominal = (TH3F*) finput->Get("h_nominal");
      TH3F* htmp_syst = (TH3F*) finput->Get("h_syst");

      foutput->cd();
      if (!h_nominal) h_nominal = (TH3F*) htmp_nominal->Clone(htmp_nominal->GetName());
      else h_nominal->Add(htmp_nominal, 1.);
      if (!h_syst) h_syst = (TH3F*) htmp_syst->Clone(htmp_syst->GetName());
      else h_syst->Add(htmp_syst, 1.);
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

    HelperFunctions::conditionalizeHistogram<TH3F>(h_nominal, 0, nullptr, false, false);
    HelperFunctions::conditionalizeHistogram<TH3F>(h_syst, 0, nullptr, false, false);

    TH3F* h_ratio = (TH3F*) h_nominal->Clone("h_ratio");
    if (theGlobalSyst==tHardJetsDn) HelperFunctions::divideHistograms(h_nominal, h_syst, h_ratio, false);
    else HelperFunctions::divideHistograms(h_syst, h_nominal, h_ratio, false);
    foutput->WriteTObject(h_ratio);

    for (int isl=0; isl<h_nominal->GetNbinsX(); isl++){
      TH2F* h_nominal_slice = HelperFunctions::getHistogramSlice(h_nominal, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_nominal->GetName(), isl));
      TH2F* h_syst_slice = HelperFunctions::getHistogramSlice(h_syst, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_syst->GetName(), isl));
      TH2F* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));

      foutput->WriteTObject(h_nominal_slice); delete h_nominal_slice;
      foutput->WriteTObject(h_syst_slice); delete h_syst_slice;
      foutput->WriteTObject(h_ratio_slice); delete h_ratio_slice;
    }

    delete h_ratio;
  }
  else if (h_nominal_reduced){
    foutput->WriteTObject(h_nominal_reduced);
    foutput->WriteTObject(h_syst_reduced);

    h_nominal_reduced->Scale(1./HelperFunctions::getHistogramIntegralAndError(h_nominal_reduced, 1, h_nominal_reduced->GetNbinsX(), 1, h_nominal_reduced->GetNbinsY(), false));
    h_syst_reduced->Scale(1./HelperFunctions::getHistogramIntegralAndError(h_syst_reduced, 1, h_syst_reduced->GetNbinsX(), 1, h_syst_reduced->GetNbinsY(), false));

    TH2F* h_ratio = (TH2F*) h_nominal_reduced->Clone("h_ratio");
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
    "WminusH_ZZTo2L2Q_POWHEG",
    "WminusH_HToWW_2LOSFilter_POWHEG",

    "WplusH_ZZTo2L2Nu_POWHEG",
    "WplusH_ZZTo2L2Q_POWHEG",
    "WplusH_HToWW_2LOSFilter_POWHEG",

    "ZH_HTo2Nu2X_2LFilter_POWHEG",
    "ZH_HTo2L2Q_2LFilter_POWHEG",
    "ZH_HTo4Q_2LFilter_POWHEG",
    "ZH_WWTo2L2Nu_POWHEG",
    "ZH_HToLNuQQ_2LFilter_POWHEG"
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


#include "produceHiggsEstimates_ZZTo2L2Nu.h"
#include "produceHiggsEstimates_ZWTo3L1Nu.h"
