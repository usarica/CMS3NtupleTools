#include <cassert>
#include <chrono>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
#include "TStyle.h"


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace IvyStreamHelpers;

  bool looperRule(BaseTreeLooper*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const&, SimpleEntry&);

  // Helpers to keep track of elapsed time
  std::vector<std::pair<TString, std::chrono::microseconds>> type_accTime_pairs;
  void addTimeDuration(TString const& strname, std::chrono::microseconds const& dur);

}
bool LooperFunctionHelpers::looperRule(BaseTreeLooper* theLooper, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& extWgt, SimpleEntry& commonEntry){
  // Define handlers
#define OBJECT_HANDLER_COMMON_DIRECTIVES
#define OBJECT_HANDLER_SIM_DIRECTIVES \
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
  IvyMELAHelpers::GMECBlock& MEblock = theLooper->getMEblock();

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
    IVYerr << "LooperFunctionHelpers::looperRule: " << #TYPE << " " << #NAME << " is not registered. Please register and re-run." << endl; \
    assert(0); \
  }
  OBJECT_HANDLER_COMMON_DIRECTIVES;
  SCALEFACTOR_HANDLER_COMMON_DIRECTIVES;
  if (!isData){
    OBJECT_HANDLER_SIM_DIRECTIVES;
    SCALEFACTOR_HANDLER_SIM_DIRECTIVES;
  }
#undef HANDLER_DIRECTIVE

  /************************/
  /* EVENT INTERPRETATION */
  /************************/
  // Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(bool, passSMPSelection) \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, m4l_true)
#define BRANCH_VECTOR_COMMANDS
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE> NAME;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  std::vector<TString> bnames_exclude;
#define BRANCH_COMMAND(TYPE, NAME) if (theGlobalSyst!=SystematicsHelpers::sNominal){ TString strtmp = #NAME; if (strtmp.EndsWith("Up") || strtmp.EndsWith("Dn")) bnames_exclude.emplace_back(strtmp); }
  BRANCH_SCALAR_COMMANDS;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND


  // Always assign the external weight first
  auto it_extWgt = extWgt.find(theGlobalSyst);
  if (it_extWgt==extWgt.cend()) it_extWgt = extWgt.find(SystematicsHelpers::nSystematicVariations);
  if (it_extWgt==extWgt.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: External normalization map does not have a proper weight assigned!" << endl;
    assert(0);
  }
  double const& extWgt_central = it_extWgt->second;

  event_wgt = extWgt_central;
  // Set NNPDF 3.0 adjustment to 1
  event_wgt_adjustment_NNPDF30 = 1;

  if (isData) return false;

  genInfoHandler->constructGenInfo(theGlobalSyst);
  auto const& genInfo = genInfoHandler->getGenInfo();
  auto const& genparticles = genInfoHandler->getGenParticles();

  {
    double genwgt_NNPDF30 = genInfo->getGenWeight(false);
    double genwgt_default = genInfo->getGenWeight(true);
    event_wgt_adjustment_NNPDF30 = (genwgt_default!=0. ? genwgt_NNPDF30 / genwgt_default : 0.);
    event_wgt *= genwgt_default;

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
            double wgt_gjets = std::max(1., 1.71691-0.001221*part->pt());;
            event_wgt *= wgt_gjets;
            break;
          }
        }
      }
      if (hasPTGExceptionRange){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
            if ((pTG_true_exception_range[0]>=0.f && part->pt()<pTG_true_exception_range[0]) || (pTG_true_exception_range[1]>=0.f && part->pt()>=pTG_true_exception_range[1])){
              event_wgt = 0;
            }
            break;
          }
        }
      }
    }

    if (event_wgt==0.f) return false;

    // Record LHE MEs and K factors
    for (auto const& it:genInfo->extras.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
    for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);
  }
  theLooper->incrementSelection("Valid gen. weights");

  m4l_true = -1;
  commonEntry.getNamedVal("KFactor_QCD_NNLO_qqVV_Bkg_arg_mass", m4l_true);
  if (m4l_true<=75. || m4l_true>=105.) return false;

  // Compute MEs
  passSMPSelection = true;
  bool computeMEs = theLooper->hasGenMEs();
  if (computeMEs){
    SimpleParticleCollection_t daughters;
    for (auto const& part:genparticles){
      if (PDGHelpers::isALepton(part->pdgId()) && part->extras.isHardProcess){
        daughters.push_back(SimpleParticle_t(part->pdgId(), ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(part->p4())));
      }
    }

    IvyMELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    IvyMELAHelpers::melaHandle->setInputEvent(&daughters, nullptr, nullptr, true);

    {
      MELACandidate* theCand = IvyMELAHelpers::melaHandle->getCurrentCandidate();
      std::vector<MELAParticle*> const& sortedDaus = theCand->getSortedDaughters();
      passSMPSelection &= (sortedDaus.at(0)->p4 + sortedDaus.at(1)->p4).M()>12.;
      passSMPSelection &= (sortedDaus.at(2)->p4 + sortedDaus.at(3)->p4).M()>4.;
      passSMPSelection &= (sortedDaus.at(0)->p4 + sortedDaus.at(3)->p4).M()>4.;
      passSMPSelection &= (sortedDaus.at(2)->p4 + sortedDaus.at(1)->p4).M()>4.;
      std::vector<double> tmp_ptlist;
      for (auto const& dau:sortedDaus){
        if (std::abs(dau->id)==13) passSMPSelection &= dau->pt()>5. && std::abs(dau->eta())<2.4;
        else if (std::abs(dau->id)==11) passSMPSelection &= dau->pt()>7. && std::abs(dau->eta())<2.5;
        else passSMPSelection &= false;
        HelperFunctions::addByHighest(tmp_ptlist, double(dau->pt()), false);
      }
      passSMPSelection &= tmp_ptlist.at(0)>20.;
      passSMPSelection &= tmp_ptlist.at(1)>10.;
    }

    MEblock.computeMELABranches();
    MEblock.pushMELABranches();

    IvyMELAHelpers::melaHandle->resetInputEvent();
  }
  else{
    MEblock.computeMELABranches();
    MEblock.pushMELABranches();
  }
  // Insert the ME values into commonEntry only when the productTreeList collection is empty.
  // Otherwise, the branches are already made.
  if (!theLooper->hasLinkedOutputTrees()){
    std::unordered_map<std::string, float> ME_values;
    MEblock.getBranchValues(ME_values);
    for (auto const& it:ME_values) commonEntry.setNamedVal(it.first, it.second);
  }

  /*********************/
  /* RECORD THE OUTPUT */
  /*********************/
#define BRANCH_COMMAND(TYPE, NAME) if (!HelperFunctions::checkListVariable<TString>(bnames_exclude, #NAME)) commonEntry.setNamedVal(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  return true;

#undef SCALEFACTOR_HANDLER_DIRECTIVES
#undef OBJECT_HANDLER_DIRECTIVES
}

void LooperFunctionHelpers::addTimeDuration(TString const& strname, std::chrono::microseconds const& dur){
  for (auto& pp:type_accTime_pairs){
    if (pp.first==strname){
      pp.second = std::chrono::duration_cast<std::chrono::microseconds>(pp.second + dur);
      return;
    }
  }
  type_accTime_pairs.emplace_back(strname, dur);
}


using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  bool computeMEs=true
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  constexpr bool useArbitraryNorm = true;
  constexpr bool useSkims = false;
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;
  std::string const systName = SystematicsHelpers::getSystName(theGlobalSyst);
  BaseTree::setRobustInputCheck(false);

  if (nchunks==1){ nchunks = 0; ichunk=0; }
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", (useSkims ? "store_skims" : "store"), prodVersion.Data()));

  constexpr float lumi = 1.;

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;
  bool isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  if (isData || theGlobalSyst!=sNominal) return;

  // Set output directory
  TString coutput_main = "output/GenZ4LEvents/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", systName.data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
  IVYout << "Created output file " << stroutput << "..." << endl;
  curdir->cd();

  // Declare handlers
  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(true);

  curdir->cd();

  // Configure the looper
  BaseTreeLooper theLooper;
  // Set chunk index
  if (nchunks>0){
    theLooper.setEventIndexRangeBySampleChunks(true);
    theLooper.setEventIndexRange(ichunk, nchunks);
  }
  // Set systematic
  theLooper.setSystematic(theGlobalSyst);
  // Set looper function
  theLooper.setLooperFunction(LooperFunctionHelpers::looperRule);
  // Set object handlers
  theLooper.addObjectHandler(&genInfoHandler);
  // Set output tree
  theLooper.addOutputTree(tout);
  // Set the MEs
  if (computeMEs) theLooper.setMatrixElementList(
    {
      "Name:SampleHypothesis Alias:<Name> Process:bkgZZ Production:ZZQQB MatrixElement:MCFM isGen:1 Cluster:Common NoBranch:1",
      "Name:QQB_S_BKG_MCFM Alias:<Name> Process:bkgZZ Production:ZZQQB_S MatrixElement:MCFM isGen:1 Cluster:Common Options:DivideP=SampleHypothesis",
      "Name:QQB_TU_BKG_MCFM Alias:<Name> Process:bkgZZ Production:ZZQQB_TU MatrixElement:MCFM isGen:1 Cluster:Common Options:DivideP=SampleHypothesis"
    },
    true
  );

  curdir->cd();

  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    TString strdsetfname = SampleHelpers::getDatasetFileName(sname);
    IVYout << "=> Accessing the input trees from " << strdsetfname << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strdsetfname, "cms3ntuple/Events", "", ""); sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    IVYout << "\t- Sample identifier (is data ? " << isData << "): " << sample_tree->sampleIdentifier << endl;

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getNEvents();
    double sum_wgts = ((isData || useArbitraryNorm) ? 1.f : 0.f);
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
      sample_tree->releaseBranch("xsec");
      xsec *= 1000.;

      bool has_lheMEweights = false;
      bool has_lheparticles = false;
      bool has_genparticles = false;
      for (auto const& bname:allbranchnames){
        if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights = true;
        else if (bname.Contains(GenInfoHandler::colName_lheparticles)) has_lheparticles = true;
        else if (bname.Contains(GenInfoHandler::colName_genparticles)) has_genparticles = true;
      }

      genInfoHandler.setAcquireGenParticles(false);
      genInfoHandler.bookBranches(sample_tree);

      // Get sum of weights
      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      bool hasCounters = true;
      {
        int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
        int bin_period = 1;
        for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
          if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
        }
        IVYout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
        for (auto const& fname:inputfilenames){
          TFile* ftmp = TFile::Open(fname, "read");
          TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
          if (!hCounters){
            hasCounters = false;
            sum_wgts = 0;
            break;
          }
          IVYout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_raw_withveto += hCounters->GetBinContent(0, 0);
          sum_wgts_raw_noveto += hCounters->GetBinContent(0, 0) / (1. - hCounters->GetBinContent(0, 1));
          ftmp->Close();
        }
        if (hasCounters) IVYout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters && useSkims){
        IVYerr << "Skims should have contained counters histograms!" << endl;
        assert(0);
      }
      if (!hasCounters){
        unsigned int n_zero_genwgts=0;
        double frac_zero_genwgts=0;

        if (useArbitraryNorm){
          sum_wgts = sum_wgts_raw_withveto = 1;
        }
        else{
          IVYout << "No counters histograms are found. Initiation loop over " << nEntries << " events to determine the sample normalization:" << endl;

          genInfoHandler.wrapTree(sample_tree);

          for (int ev=0; ev<nEntries; ev++){
            HelperFunctions::progressbar(ev, nEntries);
            sample_tree->getEvent(ev);

            genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
            auto const& genInfo = genInfoHandler.getGenInfo();

            double genwgt = genInfo->getGenWeight(true);
            double genwgt_defaultMemberZero = genInfo->extras.LHEweight_defaultMemberZero;
            if (genwgt==0.){
              n_zero_genwgts++;
              continue;
            }

            sum_wgts_raw_withveto_defaultMemberZero += genwgt_defaultMemberZero;
            sum_wgts_raw_withveto += genwgt;
          }
        }
        sum_wgts = sum_wgts_raw_withveto;
        if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;

      // Reset gen. and LHE particle settings
      genInfoHandler.setAcquireGenParticles(true);
      genInfoHandler.bookBranches(sample_tree);
    }
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> globalWeights;
    double globalWeight = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts; globalWeights[theGlobalSyst] = globalWeight;
    IVYout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;
    IVYout << "\t- Raw xsec = " << xsec << endl;
    IVYout << "\t- xsec scale * BR scale = " << xsec_scale * BR_scale << endl;
    IVYout << "\t- xsec * BR * lumi = " << xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) << endl;
    IVYout << "\t- Global weight = " << globalWeight << endl;

    sample_tree->silenceUnused();

    // Add the input tree to the looper
    theLooper.addTree(sample_tree, globalWeights);
  }

  // Loop over all events
  theLooper.loop(true);

  for (auto const& pp:LooperFunctionHelpers::type_accTime_pairs) IVYout << pp.first << " duration: " << pp.second.count() << endl;

  // No need for the inputs
  for (auto& ss:sample_trees) delete ss;

  // Write output
  foutput->cd();
  tout->writeToFile(foutput);
  delete tout;
  foutput->Close();

  curdir->cd();

  splitFileAndAddForTransfer(stroutput);
}
