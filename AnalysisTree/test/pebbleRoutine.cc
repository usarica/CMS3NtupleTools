#include "common_includes.h"


void pebbleRoutine(){
  // Run only on UCSD
  if (HostHelpers::GetHostLocation()!=HostHelpers::kUCSDT2) return;

  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  SampleHelpers::setDataPeriod("2018");

  TString strPebbleDir = "/store/user/usarica/Offshell_2L2Nu/Pebbles/PebbleTree";
  strPebbleDir = HostHelpers::GetStandardHostPathToStore(strPebbleDir, HostHelpers::kUCSDT2);
  SampleHelpers::setInputDirectory(strPebbleDir);


  BaseTree sample_tree(strPebbleDir+"/*.root", "cms3ntuple/Dilepton", "", "");
  sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier("/Pebble/PebbleTree/MINIAODSIM");
  bool const isData = SampleHelpers::checkSampleIsData(sample_tree.sampleIdentifier);

  std::vector<TString> allbranchnames;
  sample_tree.getValidBranchNamesWithoutAlias(allbranchnames, false);

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  PFCandidateHandler pfcandidateHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  SuperclusterHandler superclusterHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;

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
  }

  // Configure handlers
  pfcandidateHandler.bookBranches(&sample_tree); pfcandidateHandler.wrapTree(&sample_tree);
  muonHandler.bookBranches(&sample_tree); muonHandler.wrapTree(&sample_tree);
  electronHandler.bookBranches(&sample_tree); electronHandler.wrapTree(&sample_tree);
  photonHandler.bookBranches(&sample_tree); photonHandler.wrapTree(&sample_tree);
  jetHandler.bookBranches(&sample_tree); jetHandler.wrapTree(&sample_tree);
  isotrackHandler.bookBranches(&sample_tree); isotrackHandler.wrapTree(&sample_tree);
  vertexHandler.bookBranches(&sample_tree); vertexHandler.wrapTree(&sample_tree);

  bool hasSuperclusters = false;
  for (auto const& bname:allbranchnames){
    if (bname.BeginsWith(SuperclusterHandler::colName.data())){
      hasSuperclusters = true;
      break;
    }
  }
  if (hasSuperclusters){ superclusterHandler.bookBranches(&sample_tree); superclusterHandler.wrapTree(&sample_tree); }

  sample_tree.silenceUnused();

  const int nevents = sample_tree.getNEvents();
  size_t ev_acc=0, max_ev_acc=300000;
  MELAout << "Number of events: " << nevents << endl;
  for (int ev=0; ev<nevents; ev++){
    sample_tree.getEvent(ev);
    HelperFunctions::progressbar(ev, nevents);
    
    if (!isData){
      genInfoHandler.constructGenInfo(theGlobalSyst);
      simEventHandler.constructSimEvent();
    }

    pfcandidateHandler.constructPFCandidates(theGlobalSyst);
    auto const& pfcandidates = pfcandidateHandler.getProducts();

    muonHandler.constructMuons(theGlobalSyst, &pfcandidates);
    electronHandler.constructElectrons(theGlobalSyst, &pfcandidates);
    photonHandler.constructPhotons(theGlobalSyst, &pfcandidates);
    particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

    auto const& muons = muonHandler.getProducts();
    auto const& electrons = electronHandler.getProducts();
    auto const& photons = photonHandler.getProducts();

    jetHandler.constructJetMET(&simEventHandler, theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);

    ev_acc++;
    if (ev_acc==max_ev_acc) break;
  }
}
