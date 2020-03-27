#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


void produceGammaJetsSkims(TString strSampleSet, TString period, TString prodVersion, TString strdate="", int ichunk=0, int nchunks=0){
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop:"+prodVersion);

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  std::vector<TriggerHelpers::TriggerType> requiredTriggers{ TriggerHelpers::kSinglePho };
  std::vector<std::string> triggerCheckList = TriggerHelpers::getHLTMenus(requiredTriggers);
  //auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);

  std::vector<TString> validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t nValidDataPeriods = validDataPeriods.size();

  // Get handlers
  // Some of these are needed only to enable the branches
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  ParticleDisambiguator particleDisambiguator;

  eventFilter.setTrackDataEvents(false);
  eventFilter.setCheckUniqueDataEvent(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  TString const stroutputcore = Form("output/GammaJetsSkims/%s/%s", SampleHelpers::theDataPeriod.Data(), strdate.Data());

  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    //if (!isData && nchunks>0) return;

    TString const cinputcore = SampleHelpers::getDatasetDirectoryName(strSample);
    TString const cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    const int nEntries = sample_tree.getSelectedNEvents();

    std::vector<float> sum_wgts(nValidDataPeriods+1, 0);
    std::vector<float> sum_wgts_PUDn(nValidDataPeriods+1, 0);
    std::vector<float> sum_wgts_PUUp(nValidDataPeriods+1, 0);
    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(false);
      genInfoHandler.setAcquireGenParticles(false);
      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
      for (int ev=0; ev<nEntries; ev++){
        HelperFunctions::progressbar(ev, nEntries);
        sample_tree.getSelectedEvent(ev);

        genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
        auto const& genInfo = genInfoHandler.getGenInfo();
        float genwgt = genInfo->getGenWeight(true);

        for (size_t idp=0; idp<nValidDataPeriods+1; idp++){
          if (idp==0) SampleHelpers::setDataPeriod(period);
          else SampleHelpers::setDataPeriod(validDataPeriods.at(idp-1));

          float puwgt;
          simEventHandler.constructSimEvent(SystematicsHelpers::sNominal);
          puwgt = simEventHandler.getPileUpWeight();
          sum_wgts.at(idp) += genwgt * puwgt;
          simEventHandler.constructSimEvent(SystematicsHelpers::ePUDn);
          puwgt = simEventHandler.getPileUpWeight();
          sum_wgts_PUDn.at(idp) += genwgt * puwgt;
          simEventHandler.constructSimEvent(SystematicsHelpers::ePUUp);
          puwgt = simEventHandler.getPileUpWeight();
          sum_wgts_PUUp.at(idp) += genwgt * puwgt;
        }
        SampleHelpers::setDataPeriod(period);
      }
      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(false);
      genInfoHandler.setAcquireGenParticles(true);
      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);
    }

    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    isotrackHandler.bookBranches(&sample_tree);
    isotrackHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the rest of the handles..." << endl;
    sample_tree.silenceUnused();
    sample_tree.getSelectedTree()->SetBranchStatus("triggerObjects*", 1); // Needed for trigger matching

    // Create output
    //TString stroutput = stroutputcore; if (!strSample.BeginsWith('/')) stroutput += '/'; stroutput += strSample;
    TString stroutput = cinputcore; HelperFunctions::replaceString(stroutput, (SampleHelpers::theInputDirectory + '/' + SampleHelpers::theSamplesTag).Data(), stroutputcore.Data());
    gSystem->Exec(Form("mkdir -p %s", stroutput.Data()));
    if (nchunks>0) stroutput = stroutput + Form("/allevents_%i_of_%i", ichunk, nchunks);
    else stroutput += "/allevents";
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TDirectory* outdir = foutput->mkdir("cms3ntuple");
    outdir->cd();
    if (!isData){
      TH2F hCounters("Counters", "", 3, 0, 3, nValidDataPeriods+1, 0, nValidDataPeriods+1);
      for (size_t idp=0; idp<nValidDataPeriods+1; idp++){
        hCounters.SetBinContent(1, idp+1, sum_wgts.at(idp));
        hCounters.SetBinContent(2, idp+1, sum_wgts_PUDn.at(idp));
        hCounters.SetBinContent(3, idp+1, sum_wgts_PUUp.at(idp));
      }
      outdir->WriteTObject(&hCounters);
    }
    //sample_tree.unmuteAllBranches();
    TTree* outtree = sample_tree.getSelectedTree()->CloneTree(0);

    // Loop over the tree
    size_t n_acc=0;
    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);
      if (ev%10000==0) MELAout << sample_tree.sampleIdentifier << " events: " << n_acc << " / " << ev << " / " << nEntries << endl;

      eventFilter.constructFilters();
      //if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float trigwgt = eventFilter.getTriggerWeight(triggerCheckList);
      if (trigwgt==0.) continue;

      muonHandler.constructMuons(SystematicsHelpers::sNominal);
      electronHandler.constructElectrons(SystematicsHelpers::sNominal);
      photonHandler.constructPhotons(SystematicsHelpers::sNominal);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      auto const& muons = muonHandler.getProducts();
      size_t n_muons_veto = 0;
      for (auto const& part:muons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++; }
      if (n_muons_veto>0) continue;

      auto const& electrons = electronHandler.getProducts();
      size_t n_electrons_veto = 0;
      for (auto const& part:electrons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++; }
      if (n_electrons_veto>0) continue;

      auto const& photons = photonHandler.getProducts();
      size_t n_photons_tight = 0;
      PhotonObject* theChosenPhoton = nullptr;
      for (auto const& part:photons){
        if (ParticleSelectionHelpers::isTightParticle(part)){
          if (!theChosenPhoton) theChosenPhoton = part;
          n_photons_tight++;
        }
      }
      if (n_photons_tight!=1) continue;

      // Disable advanced trigger rejection
      /*
      jetHandler.constructJetMET(SystematicsHelpers::sNominal, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();

      float trigwgt = eventFilter.getTriggerWeight(triggerPropsCheckList, &muons, &electrons, &photons, &ak4jets, &ak8jets, pfmet);
      if (trigwgt==0.) continue;
      */

      // No b-jet veto at the moment
      /*
      bool hasBtaggedJet=false;
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          if (jet->getBtagValue()>=btagvalue_thr){
            hasBtaggedJet=true;
            break;
          }
        }
      }
      if (hasBtaggedJet) continue;
      */

      // Disable HEM veto for further studies
      //if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons, &photons, &ak4jets, &ak8jets)) continue;

      outtree->Fill();
      n_acc++;
    }

    outdir->WriteTObject(outtree);
    delete outtree;
    foutput->cd();
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  }
}
