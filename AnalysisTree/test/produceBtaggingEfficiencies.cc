#include <cassert>
#include <limits>
#include "common_includes.h"
#include "TStyle.h"


using namespace std;


using namespace SystematicsHelpers;
void produceBtaggingEfficiencies(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<BtagHelpers::BtagWPType> btagwptypes{
    BtagHelpers::kDeepCSV_Loose,
    BtagHelpers::kDeepCSV_Medium,
    BtagHelpers::kDeepCSV_Tight,
    BtagHelpers::kDeepFlav_Loose,
    BtagHelpers::kDeepFlav_Medium,
    BtagHelpers::kDeepFlav_Tight
  };
  std::unordered_map<BtagHelpers::BtagWPType, float> btagwp_type_val_map;
  for (auto const& btagwptype:btagwptypes){
    BtagHelpers::setBtagWPType(btagwptype);
    btagwp_type_val_map[btagwptype] = BtagHelpers::getBtagWP(false);
  }

  std::vector<std::string> triggerCheckList_Dilepton = TriggerHelpers::getHLTMenus(
    {
      TriggerHelpers::kDoubleMu, TriggerHelpers::kDoubleEle, TriggerHelpers::kMuEle,
      TriggerHelpers::kSingleMu, TriggerHelpers::kSingleEle
    }
  );
  std::vector<std::string> triggerCheckList_SinglePhoton = TriggerHelpers::getHLTMenus(TriggerHelpers::kSinglePho);

  // Get handlers
  GenInfoHandler genInfoHandler;
  SimEventHandler simEventHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  PhotonScaleFactorHandler photonSFHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampleList);

  TString const coutput_main = "output/BtaggingEffs/" + strdate + "/" + period;
  gSystem->mkdir(coutput_main, true);

  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (isData) continue;

    TString cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    TString coutput = SampleHelpers::getSampleIdentifier(strSample);
    HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(coutput, "_MINIAOD", "");

    TString selectedTreeName = "cms3ntuple/Events";
    TString failedTreeName;
    double sum_wgts = (isData ? 1.f : 0.f);

    std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(strSample);
    bool hasCounters = true;
    {
      int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
      int bin_period = 1;
      for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
        if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
      }
      MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
      bool firstFile = true;
      for (auto const& fname:inputfilenames){
        TFile* ftmp = TFile::Open(fname, "read");
        TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
        if (!hCounters){
          hasCounters = false;
          sum_wgts = 0;
          break;
        }
        MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
        sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
        if (firstFile){
          if (ftmp->Get("cms3ntuple/Dilepton")) selectedTreeName = "cms3ntuple/Dilepton";
          if (ftmp->Get("cms3ntuple/SinglePhoton")) failedTreeName = "cms3ntuple/SinglePhoton";
          firstFile = false;
        }
        ftmp->Close();
      }
      if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
    }

    BaseTree sample_tree(cinput, selectedTreeName, failedTreeName, "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    const int nEntries = sample_tree.getNEvents();

    simEventHandler.bookBranches(&sample_tree);
    simEventHandler.wrapTree(&sample_tree);

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    sample_tree.silenceUnused();

    {
      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      if (!hasCounters){
        MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
        for (int ev=0; ev<nEntries; ev++){
          HelperFunctions::progressbar(ev, nEntries);
          sample_tree.getSelectedEvent(ev);

          genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
          auto const& genInfo = genInfoHandler.getGenInfo();
          double genwgt = genInfo->getGenWeight(true);

          simEventHandler.constructSimEvent(theGlobalSyst);
          double puwgt = simEventHandler.getPileUpWeight();
          sum_wgts += genwgt * puwgt;
        }
      }
    }

    MELAout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;
    double wgtNorm = static_cast<const float>(nEntries) / sum_wgts;

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

    vertexHandler.bookBranches(&sample_tree);
    vertexHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    sample_tree.silenceUnused();

    MELAout << "Completed getting the rest of the handles..." << endl;

    // Create output tree
    TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    foutput->cd();
    ExtendedBinning ptbins({ 0, 20, 25, 30, 50, 75, 100, 150, 200 });
    ExtendedBinning etabins(24, -2.4, 2.4);
    etabins.addBinBoundary(-2.5); etabins.addBinBoundary(2.5);
    etabins.addBinBoundary(-5.); etabins.addBinBoundary(5.);
    std::vector<TH2F> h_All{
      TH2F("AllJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepCSV_Loose{
      TH2F("DeepCSV_LooseJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepCSV_LooseJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepCSV_LooseJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepCSV_Medium{
      TH2F("DeepCSV_MediumJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepCSV_MediumJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepCSV_MediumJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepCSV_Tight{
      TH2F("DeepCSV_TightJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepCSV_TightJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepCSV_TightJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepFlavor_Loose{
      TH2F("DeepFlavor_LooseJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_LooseJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_LooseJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepFlavor_Medium{
      TH2F("DeepFlavor_MediumJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_MediumJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_MediumJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepFlavor_Tight{
      TH2F("DeepFlavor_TightJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_TightJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_TightJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };

    // Loop over the tree
    unsigned int n_acc=0;
    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events from " << sample_tree.sampleIdentifier << ", starting from " << ev_start << " and ending at " << ev_end << "..." << endl;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);

      bool doDileptons = (ev<sample_tree.getSelectedNEvents());

      genInfoHandler.constructGenInfo(theGlobalSyst);
      auto const& genInfo = genInfoHandler.getGenInfo();
      float genwgt = genInfo->getGenWeight(true);
      if (genwgt==0.f) continue;

      simEventHandler.constructSimEvent(theGlobalSyst);
      float puwgt = simEventHandler.getPileUpWeight();
      if (puwgt==0.f) continue;

      double wgt = genwgt * puwgt * wgtNorm;

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters() || !eventFilter.hasGoodVertex()) continue;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;

      wgt *= eventFilter.getTriggerWeight((doDileptons ? triggerCheckList_Dilepton : triggerCheckList_SinglePhoton));
      if (wgt==0.f) continue;

      muonHandler.constructMuons(theGlobalSyst);
      electronHandler.constructElectrons(theGlobalSyst);
      photonHandler.constructPhotons(theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      auto const& muons = muonHandler.getProducts();
      unsigned int n_muons_veto = 0;
      float SF_muons = 1;
      for (auto const& part:muons){
        float theSF = 1;
        //if (!isData) muonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_muons *= theSF;

        if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++;
      }

      auto const& electrons = electronHandler.getProducts();
      unsigned int n_electrons_veto = 0;
      float SF_electrons = 1;
      for (auto const& part:electrons){
        float theSF = 1;
        //if (!isData) electronSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_electrons *= theSF;

        if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++;
      }

      auto const& photons = photonHandler.getProducts();
      unsigned int n_photons_tight = 0;
      float SF_photons = 1;
      PhotonObject const* theChosenPhoton = nullptr;
      for (auto const& part:photons){
        float theSF = 1;
        if (!isData) photonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_photons *= theSF;

        if (ParticleSelectionHelpers::isTightParticle(part)){
          if (!theChosenPhoton) theChosenPhoton = part;
          n_photons_tight++;
        }
      }

      wgt *= SF_muons*SF_electrons*SF_photons;
      if (wgt==0.f) continue;

      if (!doDileptons && (n_muons_veto+n_electrons_veto)>0) continue;

      dileptonHandler.constructDileptons(&muons, &electrons);
      DileptonObject* theChosenDilepton = nullptr;
      size_t nTightDilep = 0;
      for (auto const& dilepton:dileptonHandler.getProducts()){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          nTightDilep++;
        }
      }
      if (doDileptons && nTightDilep==0) continue;

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& pfmet = jetHandler.getPFMET();

      if (!eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) continue;

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)) ak4jets_tight.push_back(jet);
      }
      size_t n_ak4jets_tight = ak4jets_tight.size();
      if (n_ak4jets_tight==0) continue;

      wgt = std::abs(wgt);
      for (auto const& jet:ak4jets_tight){
        int const& jetFlavor = jet->extras.hadronFlavour;
        unsigned int iflav = 2;
        if (abs(jetFlavor)==5) iflav = 0;
        else if (abs(jetFlavor)==4) iflav = 1;

        h_All.at(iflav).Fill(jet->pt(), jet->eta(), wgt);
        for (auto const& btagwptype:btagwptypes){
          BtagHelpers::setBtagWPType(btagwptype);

          std::vector<TH2F>* hlist = nullptr;
          switch (btagwptype){
          case BtagHelpers::kDeepCSV_Loose:
            hlist = &h_DeepCSV_Loose;
            break;
          case BtagHelpers::kDeepCSV_Medium:
            hlist = &h_DeepCSV_Medium;
            break;
          case BtagHelpers::kDeepCSV_Tight:
            hlist = &h_DeepCSV_Tight;
            break;
          case BtagHelpers::kDeepFlav_Loose:
            hlist = &h_DeepFlavor_Loose;
            break;
          case BtagHelpers::kDeepFlav_Medium:
            hlist = &h_DeepFlavor_Medium;
            break;
          case BtagHelpers::kDeepFlav_Tight:
            hlist = &h_DeepFlavor_Tight;
            break;
          default:
            continue;
          }

          if (jet->getBtagValue()>=btagwp_type_val_map[btagwptype]) hlist->at(iflav).Fill(jet->pt(), jet->eta(), wgt);
        }
      }

    }

    for (unsigned int iflav=0; iflav<3; iflav++){
      foutput->WriteTObject(&h_All.at(iflav));
      //h_DeepCSV_Loose.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepCSV_Loose.at(iflav));
      //h_DeepCSV_Medium.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepCSV_Medium.at(iflav));
      //h_DeepCSV_Tight.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepCSV_Tight.at(iflav));
      //h_DeepFlavor_Loose.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepFlavor_Loose.at(iflav));
      //h_DeepFlavor_Medium.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepFlavor_Medium.at(iflav));
      //h_DeepFlavor_Tight.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepFlavor_Tight.at(iflav));
    }
    foutput->Close();
  }
}
