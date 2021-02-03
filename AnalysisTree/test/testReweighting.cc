#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include "TStyle.h"


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace MELAStreamHelpers;
  using namespace OffshellCutflow;

  std::vector<TString> selectedMEs;

  bool looperRule(BaseTreeLooper*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const&, SimpleEntry&);

}
bool LooperFunctionHelpers::looperRule(BaseTreeLooper* theLooper, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& extWgt, SimpleEntry& commonEntry){
  // Define handlers
#define OBJECT_HANDLER_COMMON_DIRECTIVES
#define OBJECT_HANDLER_SIM_DIRECTIVES \
  HANDLER_DIRECTIVE(SimEventHandler, simEventHandler) \
  HANDLER_DIRECTIVE(GenInfoHandler, genInfoHandler)
#define OBJECT_HANDLER_DIRECTIVES \
  OBJECT_HANDLER_COMMON_DIRECTIVES \
  OBJECT_HANDLER_SIM_DIRECTIVES
#define SCALEFACTOR_HANDLER_COMMON_DIRECTIVES
#define SCALEFACTOR_HANDLER_SIM_DIRECTIVES
#define SCALEFACTOR_HANDLER_DIRECTIVES \
  SCALEFACTOR_HANDLER_COMMON_DIRECTIVES \
  SCALEFACTOR_HANDLER_SIM_DIRECTIVES

  // Get the current tree
  BaseTree* currentTree = theLooper->getWrappedTree();
  if (!currentTree) return false;

  // Acquire global variables
  SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst = theLooper->getSystematic();

  // Acquire sample flags
  bool const& isData = theLooper->getCurrentTreeFlag_IsData();
  bool const& isQCD = theLooper->getCurrentTreeFlag_QCDException();
  bool const& isGJets_HT = theLooper->getCurrentTreeFlag_GJetsHTException();
  float pTG_true_exception_range[2]={ -1, -1 };
  bool hasPTGExceptionRange = theLooper->getPTGExceptionRange(pTG_true_exception_range[0], pTG_true_exception_range[1]);
  bool needGenParticleChecks = isQCD || isGJets_HT || hasPTGExceptionRange;

  // Acquire all handlers
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* NAME = nullptr;
  OBJECT_HANDLER_DIRECTIVES;
  SCALEFACTOR_HANDLER_DIRECTIVES;
#undef HANDLER_DIRECTIVE
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* tmp_##NAME = dynamic_cast<TYPE*>(handler); if (tmp_##NAME){ NAME = tmp_##NAME; continue; }
  for (auto const& handler:theLooper->getObjectHandlers()){
    OBJECT_HANDLER_COMMON_DIRECTIVES;
    if (!isData){
      OBJECT_HANDLER_SIM_DIRECTIVES;
    }
  }
  for (auto const& handler:theLooper->getSFHandlers()){
    SCALEFACTOR_HANDLER_COMMON_DIRECTIVES;
    if (!isData){
      SCALEFACTOR_HANDLER_SIM_DIRECTIVES;
    }
  }
#undef HANDLER_DIRECTIVE
#define HANDLER_DIRECTIVE(TYPE, NAME) \
  if (!NAME){ \
    MELAerr << "LooperFunctionHelpers::looperRule: " << #TYPE << " " << #NAME << " is not registered. Please register and re-run." << endl; \
    assert(0); \
  }
  OBJECT_HANDLER_COMMON_DIRECTIVES;
  SCALEFACTOR_HANDLER_COMMON_DIRECTIVES;
  if (!isData){
    OBJECT_HANDLER_SIM_DIRECTIVES;
    SCALEFACTOR_HANDLER_SIM_DIRECTIVES;
  }
#undef HANDLER_DIRECTIVE

  auto* rewgtBuilder = theLooper->getRegisteredRewgtBuilders().find("MERewgt")->second;


  /************************/
  /* EVENT INTERPRETATION */
  /************************/
  // Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(bool, invalidReweightingWgts) \
  BRANCH_COMMAND(float, sample_wgt) \
  BRANCH_COMMAND(float, sample_pairwiseNorm_wgt) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
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

  if (!isData){
    genInfoHandler->constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler->getGenInfo();
    double genwgt_NNPDF30 = genInfo->getGenWeight(false);
    double genwgt_default = genInfo->getGenWeight(true);
    event_wgt_adjustment_NNPDF30 = (genwgt_default!=0. ? genwgt_NNPDF30 / genwgt_default : 0.);
    event_wgt *= genwgt_default;
    genmet_pTmiss = genInfo->extras.genmet_met;
    genmet_phimiss = genInfo->extras.genmet_metPhi;
    auto const& genparticles = genInfoHandler->getGenParticles();

    if (needGenParticleChecks){
      if (isQCD){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState && part->pt()>=25.f){
            event_wgt = 0;
            break;
          }
        }
      }
      if (isGJets_HT){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState){
            event_wgt *= std::max(1., 1.71691-0.001221*part->pt());
            break;
          }
        }
      }
      if (hasPTGExceptionRange){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
            if ((pTG_true_exception_range[0]>=0.f && part->pt()<pTG_true_exception_range[0]) || (pTG_true_exception_range[1]>=0.f && part->pt()>=pTG_true_exception_range[1])) event_wgt = 0;
            break;
          }
        }
      }
    }

    simEventHandler->constructSimEvent();
    event_wgt *= simEventHandler->getPileUpWeight(theGlobalSyst)*simEventHandler->getL1PrefiringWeight(theGlobalSyst);

    if (event_wgt==0.f) return false;

    sample_wgt = rewgtBuilder->getOverallReweightingNormalization(currentTree);
    sample_pairwiseNorm_wgt = rewgtBuilder->getSamplePairwiseNormalization(currentTree);
    invalidReweightingWgts = !rewgtBuilder->checkWeightsBelowThreshold(currentTree);

    // Record LHE MEs and K factors
    for (auto const& it:genInfo->extras.LHE_ME_weights){
      if (HelperFunctions::checkListVariable(selectedMEs, it.first)) commonEntry.setNamedVal(it.first, it.second);
    }
    for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);
  }

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


void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  double thr_frac=-1
){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", "store", prodVersion.Data()));

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  bool isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool isVBF = strSampleSet.Contains("VBF");
  bool isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);

  /*
  if (isGG){
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M125_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M160_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M170_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M180_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M190_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M210_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M230_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M250_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M300_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M350_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M400_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M450_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M500_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M550_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M600_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M700_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M800_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M900_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M1000_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M1500_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M2000_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M2500_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("GGH_ZZTo2L2Nu_M3000_POWHEG", theGlobalSyst, sampledirs);
  }
  if (isVBF){
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M125_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M160_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M170_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M180_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M190_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M210_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M230_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M250_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M300_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M350_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M400_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M450_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M500_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M550_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M600_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M700_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M800_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M900_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M1000_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M1500_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M2000_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M2500_POWHEG", theGlobalSyst, sampledirs);
    SampleHelpers::constructSamplesList("VBF_ZZTo2L2Nu_M3000_POWHEG", theGlobalSyst, sampledirs);
  }
  */

  if (sampledirs.empty()) return;
  bool isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  if (isData) return;

  // Set output directory
  TString coutput_main = "output/ReweightingTest/SkimTrees/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  if (thr_frac>0.f) stroutput += Form("_thr_%.5f", thr_frac);
  stroutput += ".root";
  TString stroutput_txt = stroutput;
  HelperFunctions::replaceString<TString, TString const>(stroutput_txt, ".root", ".txt");
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
  MELAout << "Created output file " << stroutput << "..." << endl;
  curdir->cd();

  // Declare handlers
  SimEventHandler simEventHandler;

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(true);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(true);

  curdir->cd();

  // Configure the looper
  BaseTreeLooper theLooper;
  // Set systematic
  theLooper.setSystematic(theGlobalSyst);
  // Set looper function
  theLooper.setLooperFunction(LooperFunctionHelpers::looperRule);
  // Set object handlers
  theLooper.addObjectHandler(&simEventHandler);
  theLooper.addObjectHandler(&genInfoHandler);

  // Set output tree
  theLooper.addOutputTree(tout);

  double tol_wgt=5;
  float thr_frac_Neff=0.005;
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
  BulkReweightingBuilder rewgtBuilder_readTest(
    binning_rewgt,
    { "LHECandMass" },
    { "genHEPMCweight_default" },
    { "xsec" },
    ReweightingFunctions::getSimpleVariableBin,
    ReweightingFunctions::getSimpleWeight,
    ReweightingFunctions::getSimpleWeight
  );
  theLooper.addReweightingBuilder("MERewgt", &rewgtBuilder);
  if (isGG){
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
    /*
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
    */
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_GG_BKG_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );

    // Fill also the reading test, but in swapped order
    rewgtBuilder_readTest.addReweightingWeights(
      { "p_Gen_GG_BKG_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
    rewgtBuilder_readTest.addReweightingWeights(
      { "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
  }
  else if (isVBF){
    thr_frac_Neff = 0.01;
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_JJEW_SIG_ghv1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
    /*
    rewgtBuilder.addReweightingWeights(
    { "p_Gen_JJEW_BSI_ghv1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
    ReweightingFunctions::getSimpleWeight,
    (thr_frac>0. ? thr_frac : 0.995), tol_wgt
    );
    */
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_JJEW_BKG_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.999), tol_wgt
    );

    // Fill also the reading test, but in swapped order
    rewgtBuilder_readTest.addReweightingWeights(
      { "p_Gen_JJEW_BKG_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.999), tol_wgt
    );
    rewgtBuilder_readTest.addReweightingWeights(
      { "p_Gen_JJEW_SIG_ghv1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
  }
  else{
    thr_frac_Neff = 0.01;
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_JJEW_SIG_ghv1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
    /*
    rewgtBuilder.addReweightingWeights(
    { "p_Gen_JJEW_BSI_ghv1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
    ReweightingFunctions::getSimpleWeight,
    (thr_frac>0. ? thr_frac : 0.995), tol_wgt
    );
    */
    rewgtBuilder.addReweightingWeights(
      { "p_Gen_JJEW_BKG_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );

    // Fill also the reading test, but in swapped order
    rewgtBuilder_readTest.addReweightingWeights(
      { "p_Gen_JJEW_BKG_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
    rewgtBuilder_readTest.addReweightingWeights(
      { "p_Gen_JJEW_SIG_ghv1_1_MCFM", "p_Gen_CPStoBWPropRewgt" },
      ReweightingFunctions::getSimpleWeight,
      (thr_frac>0. ? thr_frac : 0.9995), tol_wgt
    );
  }

  // Add the list of ME weights and LHECandMass to LooperFunctionHelpers::selectedMEs so that we record them.
  for (auto const& strvar:rewgtBuilder.getBinningVars()){
    if (!HelperFunctions::checkListVariable(LooperFunctionHelpers::selectedMEs, strvar)) LooperFunctionHelpers::selectedMEs.push_back(strvar);
  }
  for (auto const& strReweightingWeightVars:rewgtBuilder.getReweightingWeightVarList()){
    for (auto const& strvar:strReweightingWeightVars){
      if (!HelperFunctions::checkListVariable(LooperFunctionHelpers::selectedMEs, strvar)) LooperFunctionHelpers::selectedMEs.push_back(strvar);
    }
  }

  curdir->cd();

  BaseTree* tree_MH125 = nullptr;
  BaseTree* tree_MHLowestOffshell = nullptr;
  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    TString strinput = SampleHelpers::getDatasetFileName(sname);
    MELAout << "Acquiring " << sname << " from input file(s) " << strinput << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strinput, "cms3ntuple/Events", "", ""); sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    float const sampleMH = SampleHelpers::findPoleMass(sample_tree->sampleIdentifier);
    if (std::abs(sampleMH-125.f)<0.8f) tree_MH125 = sample_tree;
    else if (!tree_MHLowestOffshell && sampleMH>=200.f) tree_MHLowestOffshell = sample_tree;

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getSelectedNEvents();
    bool hasTaus = false;
    double sum_wgts = (isData ? 1.f : 0.f);
    double sum_wgts_PUDn = sum_wgts;
    double sum_wgts_PUUp = sum_wgts;
    double sum_wgts_raw_noveto = sum_wgts;
    double sum_wgts_raw_withveto = sum_wgts;
    double sum_wgts_raw_withveto_defaultMemberZero = sum_wgts;
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
      simEventHandler.bookBranches(sample_tree);

      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(sampleMH>0.f);
      genInfoHandler.setAcquireGenParticles(false);
      genInfoHandler.bookBranches(sample_tree);

      sample_tree->silenceUnused();

      // Get sum of weights
      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      bool hasCounters = true;
      {
        int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
        int bin_period = 1;
        for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
          if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
        }
        MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
        for (auto const& fname:inputfilenames){
          TFile* ftmp = TFile::Open(fname, "read");
          TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
          if (!hCounters){
            hasCounters = false;
            sum_wgts = sum_wgts_PUDn = sum_wgts_PUUp = 0;
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_PUDn += hCounters->GetBinContent(2, bin_period);
          sum_wgts_PUUp += hCounters->GetBinContent(3, bin_period);
          sum_wgts_raw_withveto += hCounters->GetBinContent(0, 0);
          sum_wgts_raw_noveto += hCounters->GetBinContent(0, 0) / (1. - hCounters->GetBinContent(0, 1));
          ftmp->Close();
        }
        if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters){
        MELAout << "No counters histograms are found. Initiation loop over " << nEntries << " events to determine the sample normalization:" << endl;

        simEventHandler.wrapTree(sample_tree);
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

          simEventHandler.constructSimEvent();

          sum_wgts_raw_withveto_defaultMemberZero += genwgt_defaultMemberZero;
          sum_wgts_raw_withveto += genwgt;
          sum_wgts += genwgt * simEventHandler.getPileUpWeight(theGlobalSyst);
          sum_wgts_PUDn += genwgt * simEventHandler.getPileUpWeight(SystematicsHelpers::ePUDn);
          sum_wgts_PUUp += genwgt * simEventHandler.getPileUpWeight(SystematicsHelpers::ePUUp);
        }
        if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
      if (sampleMH>0.f) BR_scale = SampleHelpers::calculateAdjustedHiggsBREff(sname, sum_wgts_raw_withveto_defaultMemberZero, sum_wgts_raw_withveto, hasTaus);
    }
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> globalWeights;
    double globalWeight = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts; globalWeights[theGlobalSyst] = globalWeight;
    double globalWeight_PUDn = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts_PUDn; globalWeights[SystematicsHelpers::ePUDn] = globalWeight_PUDn;
    double globalWeight_PUUp = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts_PUUp; globalWeights[SystematicsHelpers::ePUUp] = globalWeight_PUUp;
    MELAout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts << " (PU dn: " << sum_wgts_PUDn << ", PU up: " << sum_wgts_PUUp << ")." << endl;
    MELAout << "\t- Raw xsec = " << xsec << endl;
    MELAout << "\t- xsec scale * BR scale = " << xsec_scale * BR_scale << endl;
    MELAout << "\t- xsec * BR * lumi = " << xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) << endl;
    MELAout << "\t- Global weight = " << globalWeight << endl;
    MELAout << "\t- Global weight (PU dn) = " << globalWeight_PUDn << endl;
    MELAout << "\t- Global weight (PU up) = " << globalWeight_PUUp << endl;

    // Reset gen. and LHE particle settings, and book those branches as well
    {
      bool has_lheMEweights = false;
      bool has_lheparticles = false;
      bool has_genparticles = false;
      for (auto const& bname:allbranchnames){
        if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights=true;
        else if (bname.Contains(GenInfoHandler::colName_lheparticles)) has_lheparticles = true;
        else if (bname.Contains(GenInfoHandler::colName_genparticles)) has_genparticles = true;
      }
      genInfoHandler.setAcquireLHEMEWeights(has_lheMEweights);
      genInfoHandler.setAcquireLHEParticles(false);
      //genInfoHandler.setAcquireLHEParticles(has_lheparticles);
      //genInfoHandler.setAcquireGenParticles(has_genparticles);
      genInfoHandler.bookBranches(sample_tree);
    }

    // Register tree
    MELAout << "\t- Registering the sample for reweighting..." << endl;
    rewgtBuilder.registerTree(sample_tree, BR_scale/sum_wgts_raw_noveto);
    rewgtBuilder_readTest.registerTree(sample_tree, BR_scale/sum_wgts_raw_noveto);

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
    else{
      //tree_normTree_pairs.emplace_back(sample_trees.at(itree), sample_trees.at(itree+1));
      //MELAout << "Normalizing mass " << sampleMH << " to mass " << SampleHelpers::findPoleMass(sample_trees.at(itree+1)->sampleIdentifier) << endl;
      continue;
    }
  }

  MELAout.open(stroutput_txt.Data());
  rewgtBuilder.setup(0, &tree_normTree_pairs, thr_frac_Neff);
  rewgtBuilder.print();
  MELAout.close();

  // Loop over all events
  theLooper.loop(true);

  TString stroutput_weights = stroutput;
  HelperFunctions::replaceString<TString, TString const>(stroutput_weights, ".root", "_ReweightingRecord.root");
  TFile* foutput_rewgtrcd = TFile::Open(stroutput_weights, "recreate");
  rewgtBuilder.writeToFile(foutput_rewgtrcd);
  foutput_rewgtrcd->Close();

  rewgtBuilder_readTest.setupFromFile(stroutput_weights);
  rewgtBuilder_readTest.print();

  // No need for the inputs
  for (auto& ss:sample_trees) delete ss;

  // Write output
  foutput->cd();
  tout->writeToFile(foutput);
  delete tout;
  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
  SampleHelpers::addToCondorTransferList(stroutput_txt);
  SampleHelpers::addToCondorTransferList(stroutput_weights);

  LooperFunctionHelpers::selectedMEs.clear();
}
