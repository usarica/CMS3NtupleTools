#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"


void testHandlers(
  TString strSampleSet, TString period, TString prodVersion,
  bool useSkims, TString skimTreeName="cms3ntuple/SinglePhoton",
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false, bool useJetOverlapStripping=false,
  // MET options
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true
){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  SampleHelpers::configure(period, Form("%s:%s", (useSkims ? "store_skims" : "store"), prodVersion.Data()));
  //if (useSkims) SampleHelpers::setInputDirectory("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/Skims/");

  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);

  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.size()!=1) return;

  BaseTree sample_tree(SampleHelpers::getDatasetFileName(sampledirs.front()), (useSkims ? skimTreeName : TString("cms3ntuple/Events")), "", "");
  sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sampledirs.front());
  bool const isData = SampleHelpers::checkSampleIsData(sample_tree.sampleIdentifier);

  std::vector<TString> allbranchnames;
  sample_tree.getValidBranchNamesWithoutAlias(allbranchnames, false);

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  PFCandidateHandler pfcandidateHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  SuperclusterHandler superclusterHandler;
  FSRHandler fsrHandler;
  JetMETHandler jetHandler; jetHandler.setVerbosity(MiscUtils::DEBUG);
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  OverlapMapHandler<MuonObject, AK4JetObject> overlapMap_muons_ak4jets;
  OverlapMapHandler<MuonObject, AK8JetObject> overlapMap_muons_ak8jets;
  OverlapMapHandler<ElectronObject, AK4JetObject> overlapMap_electrons_ak4jets;
  OverlapMapHandler<ElectronObject, AK8JetObject> overlapMap_electrons_ak8jets;
  OverlapMapHandler<PhotonObject, AK4JetObject> overlapMap_photons_ak4jets;
  OverlapMapHandler<PhotonObject, AK8JetObject> overlapMap_photons_ak8jets;
  if (useJetOverlapStripping){
    muonHandler.registerOverlapMaps(
      overlapMap_muons_ak4jets,
      overlapMap_muons_ak8jets
    );
    electronHandler.registerOverlapMaps(
      overlapMap_electrons_ak4jets,
      overlapMap_electrons_ak8jets
    );
    photonHandler.registerOverlapMaps(
      overlapMap_photons_ak4jets,
      overlapMap_photons_ak8jets
    );
    jetHandler.registerOverlapMaps(
      overlapMap_muons_ak4jets,
      overlapMap_muons_ak8jets,
      overlapMap_electrons_ak4jets,
      overlapMap_electrons_ak8jets,
      overlapMap_photons_ak4jets,
      overlapMap_photons_ak8jets
    );
  }

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;
  PhotonScaleFactorHandler photonSFHandler;
  BtagScaleFactorHandler btagSFHandler;
  METCorrectionHandler metCorrectionHandler;

  double sum_wgts = (isData ? 1.f : 0.f);
  float xsec = 1;
  float xsec_scale = 1;
  if (!isData){
    // Get cross section
    sample_tree.bookBranch<float>("xsec", 0.f);
    sample_tree.getSelectedEvent(0);
    sample_tree.getVal("xsec", xsec);
    sample_tree.releaseBranch("xsec");
    xsec *= 1000.;

    // Reset gen. and lHE particle settings
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
      genInfoHandler.setAcquireLHEParticles(has_lheparticles);
      genInfoHandler.setAcquireGenParticles(has_genparticles);
    }

    // Book branches
    simEventHandler.bookBranches(&sample_tree); simEventHandler.wrapTree(&sample_tree);
    genInfoHandler.bookBranches(&sample_tree); genInfoHandler.wrapTree(&sample_tree);

    // Get sum of weights
    std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sampledirs.front());
    double sum_wgts_raw_noveto = 0;
    double sum_wgts_raw_withveto = 0;
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
    if (!hasCounters){
      IVYerr << "Please use skim ntuples!" << endl;
      assert(0);
    }
    xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
  }
  double globalWeight = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts;
  IVYout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;
  IVYout << "\t- xsec scale = " << xsec_scale << endl;
  IVYout << "\t- Global weight = " << globalWeight << endl;

  // Configure handlers
  pfcandidateHandler.bookBranches(&sample_tree); pfcandidateHandler.wrapTree(&sample_tree);
  muonHandler.bookBranches(&sample_tree); muonHandler.wrapTree(&sample_tree);
  electronHandler.bookBranches(&sample_tree); electronHandler.wrapTree(&sample_tree);
  photonHandler.bookBranches(&sample_tree); photonHandler.wrapTree(&sample_tree);
  jetHandler.bookBranches(&sample_tree); jetHandler.wrapTree(&sample_tree);
  isotrackHandler.bookBranches(&sample_tree); isotrackHandler.wrapTree(&sample_tree);
  vertexHandler.bookBranches(&sample_tree); vertexHandler.wrapTree(&sample_tree);
  eventFilter.bookBranches(&sample_tree); eventFilter.wrapTree(&sample_tree);

  bool hasSuperclusters = false;
  for (auto const& bname:allbranchnames){
    if (bname.BeginsWith(SuperclusterHandler::colName.data())){
      hasSuperclusters = true;
      break;
    }
  }
  if (hasSuperclusters){ superclusterHandler.bookBranches(&sample_tree); superclusterHandler.wrapTree(&sample_tree); }
  fsrHandler.bookBranches(&sample_tree); fsrHandler.wrapTree(&sample_tree);

  overlapMap_muons_ak4jets.bookBranches(&sample_tree); overlapMap_muons_ak4jets.wrapTree(&sample_tree);
  overlapMap_muons_ak8jets.bookBranches(&sample_tree); overlapMap_muons_ak8jets.wrapTree(&sample_tree);
  overlapMap_electrons_ak4jets.bookBranches(&sample_tree); overlapMap_electrons_ak4jets.wrapTree(&sample_tree);
  overlapMap_electrons_ak8jets.bookBranches(&sample_tree); overlapMap_electrons_ak8jets.wrapTree(&sample_tree);
  overlapMap_photons_ak4jets.bookBranches(&sample_tree); overlapMap_photons_ak4jets.wrapTree(&sample_tree);
  overlapMap_photons_ak8jets.bookBranches(&sample_tree); overlapMap_photons_ak8jets.wrapTree(&sample_tree);

  sample_tree.silenceUnused();

#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE* NAME = nullptr; if (isData) sample_tree.getValRef(#NAME, NAME);
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE

  const int nevents = sample_tree.getNEvents();
  size_t ev_acc=0, max_ev_acc=20;
  for (int ev=0; ev<nevents; ev++){
    sample_tree.getEvent(ev);
    HelperFunctions::progressbar(ev, nevents);
    
    //if (!(*RunNumber==283453 && *LuminosityBlock==359 && *EventNumber==516553703)) continue;

    float event_wgt = 1;
    if (!isData){
      if (ev==0){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
        IVYout << "Sample cross section: " << xsec << " fb" << endl;
      }

      genInfoHandler.constructGenInfo(theGlobalSyst);
      simEventHandler.constructSimEvent();
    }

    eventFilter.constructFilters(&simEventHandler);
    if (ev==0){
      IVYerr << "Available trigger paths are as follows:" << endl;
      for (auto const& hltpath:eventFilter.getHLTPaths()) IVYout << "\t- " << hltpath->name << endl;
    }
    /*
    if (!eventFilter.hasMatchingTriggerPath(triggerCheckList)){
      IVYerr << "No matching trigger paths to " << triggerCheckList << endl;
      break;
    }
    if (eventFilter.getTriggerWeight(triggerCheckList) == 0.) continue;
    */

    IVYout << "================" << endl;
    IVYout << "Event " << ev << ":" << endl;
    if (isData){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) IVYout << "\t- " << #NAME << " = " << *NAME << endl;
      RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    }
    IVYout << "================" << endl;

    if (!isData){
      IVYout << "Sample chosen data period: "
        << simEventHandler.getChosenDataPeriod() << " with global and local random numbers "
        << simEventHandler.getRandomNumber(SimEventHandler::kDataPeriod_global) << ", "
        << simEventHandler.getRandomNumber(SimEventHandler::kDataPeriod_local)
        << "." << endl;
      IVYout << "\t- Gen. MET random number: " << simEventHandler.getRandomNumber(SimEventHandler::kGenMETSmear) << endl;
      IVYout << "\t- PU weight: " << simEventHandler.getPileUpWeight(theGlobalSyst) << endl;
      IVYout << "\t- L1 prefiring weight: " << simEventHandler.getL1PrefiringWeight(theGlobalSyst) << endl;

      auto const& genInfo = genInfoHandler.getGenInfo();
      float wgt = genInfo->getGenWeight(true);
      IVYout << "Sample gen. weight: " << wgt << endl;

      auto const& lheparticles = genInfoHandler.getLHEParticles();
      IVYout << "LHE particles:" << endl;
      for (auto const& part:lheparticles) IVYout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", status = " << part->status() << endl;

      auto const& genparticles = genInfoHandler.getGenParticles();
      IVYout << "Gen. particles:" << endl;
      for (auto const& part:genparticles){
        IVYout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", status = " << part->status() << endl;
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) IVYout << "\t\t- " << #NAME << " = " << part->extras.NAME << endl;
        GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE
      }

      float const& genmet_pTmiss = genInfo->extras.genmet_met;
      float const& genmet_phimiss = genInfo->extras.genmet_metPhi;
      IVYout << "Gen. MET (pT, phi) = (" << genmet_pTmiss << ", " << genmet_phimiss << ")" << endl;
    }

    pfcandidateHandler.constructPFCandidates(theGlobalSyst);
    auto const& pfcandidates = pfcandidateHandler.getProducts();

    muonHandler.constructMuons(theGlobalSyst, &pfcandidates);
    electronHandler.constructElectrons(theGlobalSyst, &pfcandidates);
    photonHandler.constructPhotons(theGlobalSyst, &pfcandidates);
    particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

    ParticleObject::LorentzVector_t sump4_looseObjects;

#define PARTP4PRINTCMD(PART) PART->pt() << ", " << PART->eta() << ", " << PART->phi() << ", " << PART->mass()

    auto const& muons = muonHandler.getProducts();
    IVYout << "Muons:" << endl;
    for (auto const& part:muons){
      if (ParticleSelectionHelpers::isLooseParticle(part)) sump4_looseObjects += part->p4();
      IVYout << "\t- [" << part->pdgId() << "]: " << PARTP4PRINTCMD(part) << ", (L, T) = (" << ParticleSelectionHelpers::isLooseParticle(part) << ", " << ParticleSelectionHelpers::isTightParticle(part) << ")" << endl;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) IVYout << "\t\t- " << #NAME << ": " << part->extras.NAME << endl;
      //MUON_VARIABLES;
#undef MUON_VARIABLE
    }

    auto const& electrons = electronHandler.getProducts();
    IVYout << "Electrons:" << endl;
    for (auto const& part:electrons){
      if (ParticleSelectionHelpers::isLooseParticle(part)) sump4_looseObjects += part->p4();
      IVYout << "\t- [" << part->pdgId() << "]: " << PARTP4PRINTCMD(part) << ", (L, T) = (" << ParticleSelectionHelpers::isLooseParticle(part) << ", " << ParticleSelectionHelpers::isTightParticle(part) << ")" << endl;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) IVYout << "\t\t- " << #NAME << ": " << part->extras.NAME << endl;
      //ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
    }

    auto const& photons = photonHandler.getProducts();
    IVYout << "Photons:" << endl;
    for (auto const& part:photons){
      if (ParticleSelectionHelpers::isLooseParticle(part)) sump4_looseObjects += part->p4();
      IVYout << "\t- [" << part->pdgId() << "]: " << PARTP4PRINTCMD(part) << ", (L, T) = (" << ParticleSelectionHelpers::isLooseParticle(part) << ", " << ParticleSelectionHelpers::isTightParticle(part) << ")" << endl;
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) IVYout << "\t\t- " << #NAME << ": " << part->extras.NAME << endl;
      //PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
    }

    jetHandler.constructJetMET(&simEventHandler, theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);
    auto const& ak4jets = jetHandler.getAK4Jets();
    auto const& ak8jets = jetHandler.getAK8Jets();
    auto const& pfmet = jetHandler.getPFMET();
    auto p4_pfmet = pfmet->p4(use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation);
    auto const& puppimet = jetHandler.getPFPUPPIMET();
    auto p4_puppimet = puppimet->p4(use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation);

    IVYout << "ak4 jets:" << endl;
    for (auto const& part:ak4jets){
      if (ParticleSelectionHelpers::isTightJet(part)) sump4_looseObjects += part->p4();
      IVYout << "\t- p4 = " << PARTP4PRINTCMD(part) << ", (L, T) = (" << ParticleSelectionHelpers::isLooseJet(part) << ", " << ParticleSelectionHelpers::isTightJet(part) << "), btag value = " << part->getBtagValue() << ", uncorrected pt = " << part->uncorrected_p4().Pt() << endl;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) IVYout << "\t\t- " << #NAME << ": " << part->extras.NAME << endl;
      //AK4JET_CORE_VARIABLES;
#undef AK4JET_VARIABLE
    }
    IVYout << "ak8 jets:" << endl;
    for (auto const& part:ak8jets){
      IVYout << "\t- p4 = " << PARTP4PRINTCMD(part) << ", (L, T) = (" << ParticleSelectionHelpers::isLooseJet(part) << ", " << ParticleSelectionHelpers::isTightJet(part) << ")" << endl;
    }

    sump4_looseObjects = -sump4_looseObjects;

    IVYout << "PF MET pT, phi = " << p4_pfmet.pt() << ", " << p4_pfmet.phi() << endl;
    IVYout << "\t- PF MET extra; met_Nominal, metPhi_Nominal = " << pfmet->extras.met_Nominal << ", " << pfmet->extras.metPhi_Nominal << endl;
    IVYout << "PUPPI MET pT, phi = " << p4_puppimet.pt() << ", " << p4_puppimet.phi() << endl;
    IVYout << "MET pT, phi calculated from loose objects: " << sump4_looseObjects.Pt() << ", " << sump4_looseObjects.phi() << endl;

#undef PARTP4PRINTCMD

    dileptonHandler.constructDileptons(&muons, &electrons);
    auto const& dileptons = dileptonHandler.getProducts();

    DileptonObject* theChosenDilepton = nullptr;
    for (auto const& dilepton:dileptons){
      if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
        theChosenDilepton = dilepton;
        break;
      }
    }
    IVYout << "Dileptons:" << endl;
    if (!theChosenDilepton) IVYout << "No valid OS TT dilepton pair is found." << endl;
    else IVYout << "Found " << dileptons.size() << " dileptons. The chosen dilepton pair has p4 = " << theChosenDilepton->p4() << endl;

    //IVYout << "Triggers:" << endl;
    //for (auto const& strTrigger:triggerCheckList) IVYout << "\t- Trigger weight(" << strTrigger << ") = " << eventFilter.getTriggerWeight({ strTrigger }) << endl;

    IVYout << "Event " << (eventFilter.passMETFilters(EventFilterHandler::kMETFilters_Standard) ? "passes" : "fails") << " MET filters. Available MET filters are as follows:" << endl;
    for (auto it:eventFilter.getMETFilters()) IVYout << "\t- " << it.first << ": " << it.second << endl;

    ev_acc++;
    if (ev_acc==max_ev_acc) break;
  }
}
