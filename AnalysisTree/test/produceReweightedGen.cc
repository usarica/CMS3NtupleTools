#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include "TStyle.h"


// Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(bool, invalidReweightingWgts) \
  BRANCH_COMMAND(float, sample_wgt) \
  BRANCH_COMMAND(float, lheHiggs_pt) \
  BRANCH_COMMAND(float, lheLeptonicDecay_pt) \
  BRANCH_COMMAND(float, lheLeptonicDecay_mass) \
  BRANCH_COMMAND(float, genLeptonicDecay_pt) \
  BRANCH_COMMAND(float, genLeptonicDecay_mass)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, lhemothers_id) \
  BRANCH_COMMAND(cms3_id_t, genparticles_id) \
  BRANCH_COMMAND(float, lhemothers_pz) \
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


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace MELAStreamHelpers;
  using namespace OffshellCutflow;

  std::vector<TString> selectedMEs;

  bool looperRule(BaseTreeLooper*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const&, SimpleEntry&);

}
bool LooperFunctionHelpers::looperRule(BaseTreeLooper* theLooper, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& extWgt, SimpleEntry& commonEntry){
  // Define handlers
#define OBJECT_HANDLER_DIRECTIVES \
  HANDLER_DIRECTIVE(GenInfoHandler, genInfoHandler)

  // Get the current tree
  BaseTree* currentTree = theLooper->getWrappedTree();
  if (!currentTree) return false;

  // Acquire global variables
  SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst = theLooper->getSystematic();

  // Acquire all handlers
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* NAME = nullptr;
  OBJECT_HANDLER_DIRECTIVES;
#undef HANDLER_DIRECTIVE
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* tmp_##NAME = dynamic_cast<TYPE*>(handler); if (tmp_##NAME){ NAME = tmp_##NAME; continue; }
  for (auto const& handler:theLooper->getObjectHandlers()){
    OBJECT_HANDLER_DIRECTIVES;
  }
#undef HANDLER_DIRECTIVE

  auto* rewgtBuilder = theLooper->getRegisteredRewgtBuilders().find("MERewgt")->second;


  /************************/
  /* EVENT INTERPRETATION */
  /************************/
#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE> NAME;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  // Always assign the external weight first
  auto it_extWgt = extWgt.find(theGlobalSyst);
  if (it_extWgt==extWgt.cend()) it_extWgt = extWgt.find(SystematicsHelpers::nSystematicVariations);
  if (it_extWgt==extWgt.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: External normalization map does not have a proper weight assigned!" << endl;
    assert(0);
  }
  double const& extWgt_central = it_extWgt->second;

  event_wgt = extWgt_central;
  // Set NNPDF 3.0 adjustment to 1
  event_wgt_adjustment_NNPDF30 = 1;

  genInfoHandler->constructGenInfo(theGlobalSyst);
  auto const& genInfo = genInfoHandler->getGenInfo();
  double genwgt_NNPDF30 = genInfo->getGenWeight(false);
  double genwgt_default = genInfo->getGenWeight(true);
  event_wgt_adjustment_NNPDF30 = (genwgt_default!=0. ? genwgt_NNPDF30 / genwgt_default : 0.);
  event_wgt *= genwgt_default;
  if (event_wgt==0.f) return false;

  bool hasTaus = false;
  unsigned int n_leps_nus=0;
  ParticleObject::LorentzVector_t p4_lheHiggs;
  ParticleObject::LorentzVector_t p4_lheLeptonicDecay;
  auto const& lheparticles = genInfoHandler->getLHEParticles();
  for (auto const& part:lheparticles){
    if (part->status()==2 && PDGHelpers::isAHiggs(part->pdgId())) p4_lheHiggs = part->p4();
    else if (part->status()==1 && (PDGHelpers::isALepton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()))){
      p4_lheLeptonicDecay = p4_lheLeptonicDecay + part->p4();
      if (std::abs(part->pdgId())==15) hasTaus = true;
      n_leps_nus++;
    }
    else if (part->status()==-1){
      lhemothers_id.push_back(part->pdgId());
      lhemothers_pz.push_back(part->pz());
    }
  }
  if (n_leps_nus!=4) return false;
  if (hasTaus) return false;
  lheHiggs_pt = p4_lheHiggs.Pt();
  lheLeptonicDecay_pt = p4_lheLeptonicDecay.Pt();
  lheLeptonicDecay_mass = p4_lheLeptonicDecay.M();

  ParticleObject::LorentzVector_t p4_genLeptonicDecay;
  auto const& genparticles = genInfoHandler->getGenParticles();
  for (auto const& part:genparticles){
    if (part->extras.isPromptFinalState && (PDGHelpers::isAPhoton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()) || PDGHelpers::isALepton(part->pdgId()))){
      p4_genLeptonicDecay += part->p4();
      genparticles_id.push_back(part->pdgId());
      genparticles_pt.push_back(part->pt());
      genparticles_eta.push_back(part->eta());
      genparticles_phi.push_back(part->phi());
      genparticles_mass.push_back(part->mass());
    }
  }
  genLeptonicDecay_pt = p4_genLeptonicDecay.Pt();
  genLeptonicDecay_mass = p4_genLeptonicDecay.M();

  auto const& genak4jets = genInfoHandler->getGenAK4Jets();
  for (auto const& jet:genak4jets){
    if (jet->pt()<20.f) continue;
    genak4jets_pt.push_back(jet->pt());
    genak4jets_eta.push_back(jet->eta());
    genak4jets_phi.push_back(jet->phi());
    genak4jets_mass.push_back(jet->mass());
  }

  sample_wgt = rewgtBuilder->getOverallReweightingNormalization(currentTree);
  invalidReweightingWgts = !rewgtBuilder->checkWeightsBelowThreshold(currentTree);

  // Record LHE MEs and K factors
  for (auto const& it:genInfo->extras.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
  for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);

  /*********************/
  /* RECORD THE OUTPUT */
  /*********************/
#define BRANCH_COMMAND(TYPE, NAME) commonEntry.setNamedVal(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  return true;

#undef SCALEFACTOR_HANDLER_DIRECTIVES
#undef OBJECT_HANDLER_DIRECTIVES
}


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


void produceReweightedGen(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate
){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strSampleSet.Contains("/MINIAOD")){
    MELAerr << "Processing single samples is not the design goal of produceReweightingRecords." << endl;
    return;
  }

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", "store", prodVersion.Data()));

  constexpr float lumi = 1;

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  bool isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool isVBF = strSampleSet.Contains("VBF");
  bool const hasDirectHWW = (
    strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH")
    ||
    SampleHelpers::isHiggsToWWDecay(SampleHelpers::getHiggsSampleDecayMode(strSampleSet))
    );

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);

  if (sampledirs.empty()) return;
  bool isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  if (isData) return;

  // Set output directory
  TString coutput_main = "output/ReweightedGenTrees/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  stroutput += ".root";
  TString stroutput_txt = stroutput;
  HelperFunctions::replaceString<TString, TString const>(stroutput_txt, ".root", ".txt");
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
  MELAout << "Created output file " << stroutput << "..." << endl;
  curdir->cd();

  // Declare handlers
  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(true);
  genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireGenParticles(true);

  curdir->cd();

  // Configure the looper
  BaseTreeLooper theLooper;
  // Set systematic
  theLooper.setSystematic(theGlobalSyst);
  // Set looper function
  theLooper.setLooperFunction(LooperFunctionHelpers::looperRule);
  // Set object handlers
  theLooper.addObjectHandler(&genInfoHandler);

  // Set output tree
  theLooper.addOutputTree(tout);

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
  BulkReweightingBuilder rewgtBuilder(
    binning_rewgt,
    { "LHECandMass" },
    { "genHEPMCweight_default" },
    { "xsec" },
    ReweightingFunctions::getSimpleVariableBin,
    ReweightingFunctions::getSimpleWeight,
    ReweightingFunctions::getSimpleWeight
  );
  theLooper.addReweightingBuilder("MERewgt", &rewgtBuilder);

  // Acquire the MEs
  std::vector<TString> strMEs;
  {
    PhysicsProcessHandler* proc_handler = getPhysicsProcessHandler(strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);
    for (unsigned int iac=0; iac<(unsigned int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis hypo = (ACHypothesisHelpers::ACHypothesis) iac;
      if (hasDirectHWW && hypo==ACHypothesisHelpers::kL1ZGs) continue;
      std::vector<TString> strMEs_hypo = proc_handler->getMELAHypothesisWeights(hypo, false);
      HelperFunctions::appendVector(strMEs, strMEs_hypo);
    }
    delete proc_handler;
  }

  constexpr double tol_wgt = 5;
  float const thr_frac_Neff = (isGG ? 0.005 : 0.01);
  for (auto const& strME:strMEs){
    MELAout << "Registering ME " << strME << " for reweighting..." << endl;
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

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getSelectedNEvents();
    bool hasTaus = false;
    double sum_wgts_raw_noveto = 0;
    double sum_wgts_raw_withveto = 0;
    double sum_wgts_raw_withveto_defaultMemberZero = 0;
    float xsec = 1;
    float xsec_scale = 1;
    float BR_scale = 1;
    if (!isData){
      // Get cross section
      sample_tree->bookBranch<float>("xsec", 0.f);
      sample_tree->getSelectedEvent(0);
      sample_tree->getVal("xsec", xsec);
      xsec *= 1000.;

      // Book branches
      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(true);
      genInfoHandler.setAcquireGenParticles(false);
      genInfoHandler.bookBranches(sample_tree);

      sample_tree->silenceUnused();

      // Get sum of weights
      {
        MELAout << "No counters histograms are found. Initiation loop over " << nEntries << " events to determine the sample normalization:" << endl;

        genInfoHandler.wrapTree(sample_tree);

        unsigned int n_zero_genwgts=0;
        double frac_zero_genwgts=0;
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
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
      if (sampleMH>0.f) BR_scale = SampleHelpers::calculateAdjustedHiggsBREff(sname, sum_wgts_raw_withveto_defaultMemberZero, sum_wgts_raw_withveto, hasTaus);
    }
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> globalWeights;
    double globalWeight = xsec * xsec_scale * BR_scale / sum_wgts_raw_withveto; globalWeights[theGlobalSyst] = globalWeight;
    MELAout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts_raw_withveto << "." << endl;
    MELAout << "\t- Raw xsec = " << xsec << endl;
    MELAout << "\t- xsec scale * BR scale = " << xsec_scale * BR_scale << endl;
    MELAout << "\t- xsec * BR = " << xsec * xsec_scale * BR_scale << endl;
    MELAout << "\t- Global weight = " << globalWeight << endl;

    // Reset gen. and LHE particle settings, and book those branches as well
    {
      genInfoHandler.setAcquireLHEMEWeights(true);
      genInfoHandler.setAcquireLHEParticles(true);
      genInfoHandler.setAcquireGenParticles(true);
      genInfoHandler.setAcquireGenAK4Jets(true);
      // Specify this flag to omit V->qq->jj
      genInfoHandler.setDoGenJetsVDecayCleaning(false);
      // Do not clean jets, allow user to do this on their own
      genInfoHandler.setDoGenJetsCleaning(false);
      genInfoHandler.bookBranches(sample_tree);
    }

    // Register tree
    MELAout << "\t- Registering the sample for reweighting..." << endl;
    rewgtBuilder.registerTree(sample_tree, xsec_scale * BR_scale / sum_wgts_raw_withveto);

    sample_tree->silenceUnused();

    MELAout << "\t- Registering the sample to the looper..." << endl;
    // Add the input tree to the looper
    theLooper.addTree(sample_tree, globalWeights);
  }

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

  rewgtBuilder.setup(0, &tree_normTree_pairs, thr_frac_Neff);

  // Loop over all events
  theLooper.loop(true);

  // No need for the inputs
  for (auto& ss:sample_trees) delete ss;

  // Write output
  foutput->cd();
  tout->writeToFile(foutput);
  delete tout;
  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void makePlots(TString strSampleSet, TString period, TString strdate){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  //bool isVVVV = strSampleSet.Contains("VBF") || strSampleSet.Contains("ZH") || strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH");
  MELAout << "Sample is " << (isGG ? "a gg" : "an EW") << " sample." << endl;

  // Set output directory
  TString cinput_main = "output/ReweightedGenTrees/" + strdate + "/" + period;
  TString coutput_main = "output/ReweightedGenTrees/" + strdate + "/" + period + "/Plots";

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  MELAout << "Created output file " << stroutput << "..." << endl;
  foutput->cd();
  
  TString strinput = cinput_main + "/" + coutput + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
  MELAout << "Acquiring input " << strinput << "..." << endl;
  BaseTree* tin = new BaseTree(strinput, "SkimTree", "", "");
  if (!tin->isValid()){
    MELAerr << "\t- Failed to acquire." << endl;
    delete tin;
    foutput->Close();
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
  float* val_ME_SIG = nullptr;
  if (isGG){
    val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second;
    val_ME_SIG = ME_Kfactor_values.find("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM")->second;
  }
  else{
    val_ME_SIG = ME_Kfactor_values.find("p_Gen_JJEW_SIG_ghv1_1_MCFM")->second;
  }
  float* val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
  float* LHECandMass = ME_Kfactor_values.find("LHECandMass")->second;

  bool hasError = false;
  if (!val_ME_SIG){
    MELAerr << "val_ME_SIG is null!" << endl;
    hasError = true;
  }
  if (!val_ME_CPS){
    MELAerr << "val_ME_CPS is null!" << endl;
    hasError = true;
  }
  if (isGG && !val_Kfactor_QCD){
    MELAerr << "val_Kfactor_QCD is null!" << endl;
    hasError = true;
  }
  if (hasError){
    delete tin;
    foutput->Close();
    exit(1);
  }

  std::vector<TString> const strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  std::vector<TString> const strCatLabels{ "N_{j}=0", "N_{j}=1", "N_{j} #geq 2" };
  size_t const nCats = strCatNames.size();

  ExtendedBinning binning_mass("genmass", "m_{ZZ} (GeV)");
  constexpr bool useConstBinning = false;
  if (!useConstBinning){
    binning_mass.addBinBoundary(100);
    binning_mass.addBinBoundary(124);
    binning_mass.addBinBoundary(126);
    binning_mass.addBinBoundary(150);
    {
      double bin_low_edge = 180;
      while (bin_low_edge<=500.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 20.;
      }
      while (bin_low_edge<=1000.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 50.;
      }
      while (bin_low_edge<=2000.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 100.;
      }
      while (bin_low_edge<=3000.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 500.;
      }
    }
  }
  else{
    double bin_low_edge = 100;
    while (bin_low_edge<=3000.){
      binning_mass.addBinBoundary(bin_low_edge);
      bin_low_edge += 25.;
    }
  }

  std::vector<TH1F> hlist_genmass; hlist_genmass.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetpt; hlist_genak4jetpt.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetselectedpt; hlist_genak4jetselectedpt.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetselected_leadingpt; hlist_genak4jetselected_leadingpt.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetselected_subleadingpt; hlist_genak4jetselected_subleadingpt.reserve(strCatNames.size());
  for (unsigned short icat=0; icat<nCats; icat++){
    hlist_genmass.emplace_back(
      Form("h_genmass_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      binning_mass.getNbins(), binning_mass.getBinning()
    );
    hlist_genak4jetpt.emplace_back(
      Form("h_genak4jetpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
    hlist_genak4jetselectedpt.emplace_back(
      Form("h_genak4jetselectedpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
    hlist_genak4jetselected_leadingpt.emplace_back(
      Form("h_genak4jetselected_leadingpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
    hlist_genak4jetselected_subleadingpt.emplace_back(
      Form("h_genak4jetselected_subleadingpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
  }

  float sum_wgts = 0;
  float sum_wgts_accepted = 0;
  float sum_wgts_rejected_Nleptons_lt_4 = 0;
  float sum_wgts_rejected_Nleptons_gt_4 = 0;
  float sum_wgts_rejected_noZZ = 0;
  float sum_wgts_hasHardPhotons = 0;
  int n_printouts = -1;
  int nEntries = tin->getNEvents();
  MELAout << "Looping over " << nEntries << " events:" << endl;
  for (int ev=0; ev<nEntries; ev++){
    tin->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt = (*event_wgt) * (*sample_wgt) * (*invalidReweightingWgts ? 0.f : 1.f) * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f) * (*val_ME_SIG) * (*val_ME_CPS);
    sum_wgts += wgt;

    unsigned int n_genparticles = (*genparticles_id)->size();
    std::vector<MELAParticle> genparticles; genparticles.reserve(n_genparticles);
    for (unsigned int ipart=0; ipart<n_genparticles; ipart++){
      TLorentzVector tmp_p4;
      tmp_p4.SetPtEtaPhiM((*genparticles_pt)->at(ipart), (*genparticles_eta)->at(ipart), (*genparticles_phi)->at(ipart), (*genparticles_mass)->at(ipart));
      genparticles.emplace_back((*genparticles_id)->at(ipart), tmp_p4);
    }

    int sumid_genleptons_selected = 0;
    TLorentzVector sump4_genleptons_selected;
    std::vector<MELAParticle const*> genleptons_selected; genleptons_selected.reserve(n_genparticles);
    std::vector<MELAParticle const*> genphotons_selected; genphotons_selected.reserve(n_genparticles);
    for (auto const& part:genparticles){
      int pid = part.id;
      if (PDGHelpers::isALepton(pid) && part.pt()>=5. && std::abs(part.eta())<(std::abs(pid)==11 ? 2.5 : 2.4)){
        sumid_genleptons_selected += part.id;
        sump4_genleptons_selected += part.p4;
        genleptons_selected.push_back(&part);
      }
      else if (PDGHelpers::isAPhoton(pid) && part.pt()>=20. && std::abs(part.eta())<2.5) genphotons_selected.push_back(&part);
    }

    if (genleptons_selected.size()!=4){
      if (n_printouts>=0 && n_printouts<100){
        MELAout << "Event " << ev << " is rejected because it has " << genleptons_selected.size() << " leptons." << endl;
        n_printouts++;
      }
      if (genleptons_selected.size()<4) sum_wgts_rejected_Nleptons_lt_4 += wgt;
      if (genleptons_selected.size()>4) sum_wgts_rejected_Nleptons_gt_4 += wgt;
      continue;
    }
    if (sumid_genleptons_selected!=0){
      if (n_printouts>=0 && n_printouts<100){
        MELAout << "Event " << ev << " is rejected because it doesn't have a proper ZZ candidate." << endl;
        n_printouts++;
      }
      sum_wgts_rejected_noZZ += wgt;
      continue;
    }
    sum_wgts_accepted += wgt;

    if (!genphotons_selected.empty()){
      bool hasHardPhotons = false;
      for (auto const& photon:genphotons_selected){
        bool isSeparated = true;
        for (auto const& lepton:genleptons_selected){
          if (photon->deltaR(lepton)<0.4){
            isSeparated = false;
            break;
          }
        }
        if (isSeparated){
          hasHardPhotons = true;
          break;
        }
      }
      if (hasHardPhotons) sum_wgts_hasHardPhotons += wgt;
    }

    unsigned int n_genak4jets = (*genak4jets_pt)->size();
    std::vector<MELAParticle> genak4jets; genak4jets.reserve(n_genak4jets);
    for (unsigned int ipart=0; ipart<n_genak4jets; ipart++){
      TLorentzVector tmp_p4;
      tmp_p4.SetPtEtaPhiM((*genak4jets_pt)->at(ipart), (*genak4jets_eta)->at(ipart), (*genak4jets_phi)->at(ipart), (*genak4jets_mass)->at(ipart));
      genak4jets.emplace_back(0, tmp_p4);
    }

    std::vector<MELAParticle const*> genak4jets_selected; genak4jets_selected.reserve(genak4jets.size());
    for (auto const& jet:genak4jets){
      if (jet.pt()<30. || std::abs(jet.eta())>=4.7) continue;
      bool doSkip = false;
      for (auto const& part:genleptons_selected){
        if (jet.deltaR(part)<0.4){ doSkip = true; break; }
      }
      if (!doSkip) genak4jets_selected.push_back(&jet);
    }

    unsigned int icat = std::min(genak4jets_selected.size(), nCats-1);

    hlist_genmass.at(icat).Fill(
      //*LHECandMass,
      sump4_genleptons_selected.M(),
      wgt
    );
    if (/*(*LHECandMass)*/sump4_genleptons_selected.M()>200.){
      for (auto const& jet:genak4jets) hlist_genak4jetpt.at(icat).Fill(jet.pt(), wgt);
      for (auto const& jet:genak4jets_selected) hlist_genak4jetselectedpt.at(icat).Fill(jet->pt(), wgt);
      if (genak4jets_selected.size()>=1){
        hlist_genak4jetselected_leadingpt.at(icat).Fill(genak4jets_selected.at(0)->pt(), wgt);
        if (genak4jets_selected.size()>=2){
          hlist_genak4jetselected_subleadingpt.at(icat).Fill(genak4jets_selected.at(1)->pt(), wgt);
        }
      }
    }
  }

  MELAout << "Sum of weights: " << sum_wgts << endl;
  MELAout << "\t- Fraction of events accepted: " << sum_wgts_accepted / sum_wgts << endl;
  MELAout << "\t- Fraction of events rejected because selected number of leptons < 4: " << sum_wgts_rejected_Nleptons_lt_4 / sum_wgts << endl;
  MELAout << "\t- Fraction of events rejected because selected number of leptons > 4: " << sum_wgts_rejected_Nleptons_gt_4 / sum_wgts << endl;
  MELAout << "\t- Fraction of events rejected because no ZZ candidate can be constructed: " << sum_wgts_rejected_noZZ / sum_wgts << endl;
  MELAout << "\t- Fraction of accepted events with at least one hard photon: " << sum_wgts_hasHardPhotons / sum_wgts_accepted << endl;

  for (auto& hh:hlist_genmass){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetselectedpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetselected_leadingpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetselected_subleadingpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }

  delete tin;
  foutput->Close();
}
