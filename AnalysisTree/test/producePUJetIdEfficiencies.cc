#include <cassert>
#include <limits>
#include "common_includes.h"
#include "TStyle.h"


using namespace std;


using namespace SystematicsHelpers;
void producePUJetIdEfficiencies(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool applyTightLeptonVetoIdToAK4Jets=false
){
  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  AK4JetSelectionHelpers::setPUIdWP(AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<AK4JetSelectionHelpers::SelectionBits> pujetidwpbits{
    AK4JetSelectionHelpers::kTightPUJetId,
    AK4JetSelectionHelpers::kMediumPUJetId,
    AK4JetSelectionHelpers::kLoosePUJetId
  };

  std::vector<std::string> triggerCheckList_Dilepton = TriggerHelpers::getHLTMenus(
    {
      TriggerHelpers::kDoubleMu,
      TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt,
      TriggerHelpers::kMuEle,
      TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
      TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
    }
  );
  std::vector<std::string> triggerCheckList_SinglePhoton = TriggerHelpers::getHLTMenus(TriggerHelpers::kSinglePho);

  // Binning for efficiencies
  ExtendedBinning ptbins({ AK4JetSelectionHelpers::ptThr_gen, AK4JetSelectionHelpers::ptThr_btag, 25, 30, 40, AK4JetSelectionHelpers::ptThr_PUId });
  ExtendedBinning etabins({ -AK4JetSelectionHelpers::etaThr_PUId, -3, -2.75, -2.5, -2, -1.479, 0, 1.479, 2, 2.5, 2.75, 3, AK4JetSelectionHelpers::etaThr_PUId });

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

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;
  PhotonScaleFactorHandler photonSFHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampleList);

  TString const coutput_main = "output/PUJetIdEffs/" + strdate + "/" + period
    + "/AK4Jets"
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId");
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

    TString selectedTreeName;
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

    if (selectedTreeName == "" && failedTreeName == ""){
      MELAerr << "Both trees from " << cinput << " have an empty name!" << endl;
      continue;
    }

    MELAout << "Acquiring trees \"" << selectedTreeName << "\", \"" << failedTreeName << "\" from " << cinput << endl;
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

          simEventHandler.constructSimEvent();
          double puwgt = simEventHandler.getPileUpWeight(theGlobalSyst);
          sum_wgts += genwgt * puwgt;
        }
      }
    }

    MELAout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;

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
    std::vector<TH2F> h_All{
      TH2F("AllJets_mistagged", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_matched", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_PUJetId_Loose{
      TH2F("PUJetId_LooseJets_mistagged", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("PUJetId_LooseJets_matched", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_PUJetId_Medium{
      TH2F("PUJetId_MediumJets_mistagged", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("PUJetId_MediumJets_matched", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_PUJetId_Tight{
      TH2F("PUJetId_TightJets_mistagged", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("PUJetId_TightJets_matched", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    for (unsigned int iflav=0; iflav<h_All.size(); iflav++){
      h_All.at(iflav).Sumw2();
      h_PUJetId_Loose.at(iflav).Sumw2();
      h_PUJetId_Medium.at(iflav).Sumw2();
      h_PUJetId_Tight.at(iflav).Sumw2();
    }

    // Loop over the tree
    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events from " << sample_tree.sampleIdentifier << ", starting from " << ev_start << " and ending at " << ev_end << "..." << endl;
    size_t n_acc=0;
    double sum_acc_wgts=0;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);

      bool doDileptons = (ev<sample_tree.getSelectedNEvents());

      genInfoHandler.constructGenInfo(theGlobalSyst);
      auto const& genInfo = genInfoHandler.getGenInfo();
      float genwgt = genInfo->getGenWeight(true);
      if (genwgt==0.f) continue;

      simEventHandler.constructSimEvent();
      float pul1wgt = simEventHandler.getPileUpWeight(theGlobalSyst)*simEventHandler.getL1PrefiringWeight(theGlobalSyst);
      if (pul1wgt==0.f) continue;

      double wgt = genwgt * pul1wgt;

      eventFilter.constructFilters(&simEventHandler);
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters(EventFilterHandler::kMETFilters_Standard)) continue;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;

      wgt *= eventFilter.getTriggerWeight((doDileptons ? triggerCheckList_Dilepton : triggerCheckList_SinglePhoton));
      if (wgt==0.f) continue;

      muonHandler.constructMuons(theGlobalSyst, nullptr);
      electronHandler.constructElectrons(theGlobalSyst, nullptr);
      photonHandler.constructPhotons(theGlobalSyst, nullptr);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      auto const& muons = muonHandler.getProducts();
      unsigned int n_muons_veto = 0;
      float SF_muons = 1;
      for (auto const& part:muons){
        float theSF = 1;
        if (!isData) muonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_muons *= theSF;

        if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++;
      }

      auto const& electrons = electronHandler.getProducts();
      unsigned int n_electrons_veto = 0;
      float SF_electrons = 1;
      for (auto const& part:electrons){
        float theSF = 1;
        if (!isData) electronSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_electrons *= theSF;

        if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++;
      }

      auto const& photons = photonHandler.getProducts();
      unsigned int n_photons_tight = 0;
      unsigned int n_photons_veto = 0;
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
        else if (ParticleSelectionHelpers::isVetoParticle(part)) n_photons_veto++;
      }

      wgt *= SF_muons*SF_electrons*SF_photons;
      if (wgt==0.f) continue;

      isotrackHandler.constructIsotracks(&muons, &electrons);
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotrackHandler.getProducts()){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;

      if (!doDileptons && !((n_muons_veto+n_electrons_veto)==0 && n_photons_veto==0 && n_photons_tight==1)) continue;

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

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons, nullptr);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();

      if (!eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) continue;

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      for (auto const& jet:ak4jets){
        if (
          jet->testSelectionBit(AK4JetSelectionHelpers::bit_preselectionTight_id)
          &&
          (!applyTightLeptonVetoIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
          ) ak4jets_tight.push_back(jet);
      }
      size_t n_ak4jets_tight = ak4jets_tight.size();
      if (n_ak4jets_tight==0) continue;

      wgt = std::abs(wgt);
      for (auto const& jet:ak4jets_tight){
        unsigned short iflav = static_cast<unsigned short>(jet->extras.is_genMatched_fullCone);

        float const jetpt = jet->pt();
        float const jeteta = jet->eta();

        h_All.at(iflav).Fill(jetpt, jeteta, wgt);
        for (auto const& wpbit:pujetidwpbits){
          std::vector<TH2F>* hlist = nullptr;
          switch (wpbit){
          case AK4JetSelectionHelpers::kLoosePUJetId:
            hlist = &h_PUJetId_Loose;
            break;
          case AK4JetSelectionHelpers::kMediumPUJetId:
            hlist = &h_PUJetId_Medium;
            break;
          case AK4JetSelectionHelpers::kTightPUJetId:
            hlist = &h_PUJetId_Tight;
            break;
          default:
            continue;
          }
          if (jet->testSelectionBit(wpbit)) hlist->at(iflav).Fill(jetpt, jeteta, wgt);
        }
      }

      sum_acc_wgts += wgt;
      n_acc++;
    }

    for (unsigned int iflav=0; iflav<h_All.size(); iflav++){
      h_All.at(iflav).Scale(((double) n_acc)/sum_acc_wgts);
      h_PUJetId_Loose.at(iflav).Scale(((double) n_acc)/sum_acc_wgts);
      h_PUJetId_Medium.at(iflav).Scale(((double) n_acc)/sum_acc_wgts);
      h_PUJetId_Tight.at(iflav).Scale(((double) n_acc)/sum_acc_wgts);

      foutput->WriteTObject(&h_All.at(iflav));
      foutput->WriteTObject(&h_PUJetId_Loose.at(iflav));
      foutput->WriteTObject(&h_PUJetId_Medium.at(iflav));
      foutput->WriteTObject(&h_PUJetId_Tight.at(iflav));
    }
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  }
}

bool checkGoodHistogram(TH2F const* hist){
  if (!HelperFunctions::checkHistogramIntegrity(hist)) return false;

  for (int ix=0; ix<=hist->GetNbinsX()+1; ix++){
    for (int iy=0; iy<=hist->GetNbinsY()+1; iy++){
      double val=hist->GetBinContent(ix, iy);
      double err=hist->GetBinError(ix, iy);
      double Neff=(err<=0. ? 0. : std::pow(val/err, 2));
      if (err>0. && Neff<1.){
        MELAerr << "checkGoodHistogram: " << hist->GetName() << " integrity check was good, but bin (" << ix << ", " << iy << ") has too little Neff = (" << val << "/" << err << ")^2 = " << Neff << endl;
        return false;
      }
    }
  }

  return true;
}

void getFinalEfficiencies(
  TString period, TString strdate,
  bool applyTightLeptonVetoIdToAK4Jets=false
){
  TString const cinput_main =
    TString("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/")
    + "output/PUJetIdEffs/"
    + strdate + "/" + period
    + "/AK4Jets"
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId");
  TString coutput_main = "output/PUJetIdEffs/" + strdate + "/" + period
    + "/AK4Jets"
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/FinalEffs";
  gSystem->mkdir(coutput_main, true);

  std::vector<TString> hnames{
    "AllJets",
    "PUJetId_LooseJets",
    "PUJetId_MediumJets",
    "PUJetId_TightJets"
  };
  std::vector<TString> strmatches{ "mistagged", "matched" };
  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts{
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp
  };

  TFile* foutput = TFile::Open(coutput_main + "/Final_PUJetId_Efficiencies_AllMC.root", "recreate");

  auto infiles = SampleHelpers::lsdir(cinput_main);
  for (auto const& syst:allowedSysts){
    TString systname = SystematicsHelpers::getSystName(syst).data();
    TString file_suffix = systname + ".root";
    std::vector<std::vector<TH2F*>> hlist(strmatches.size(), std::vector<TH2F*>(hnames.size(), nullptr));

    // Merge all input histograms
    bool firstFile = true;
    for (auto const& fname:infiles){
      if (!fname.Contains(file_suffix)) continue;
      TString cinput = cinput_main + "/" + fname;
      MELAout << "Reading " << cinput << "..." << endl;
      TFile* finput = TFile::Open(cinput, "read");

      std::vector<std::vector<TH2F*>> htmplist(strmatches.size(), std::vector<TH2F*>(hnames.size(), nullptr));
      for (unsigned short im=0; im<strmatches.size(); im++){
        for (unsigned short iwp=0; iwp<hnames.size(); iwp++){
          finput->cd();
          htmplist.at(im).at(iwp) = (TH2F*) finput->Get(Form("%s_%s", hnames.at(iwp).Data(), strmatches.at(im).Data()));
          foutput->cd();
        }
      }

      // Check histogram integrity
      bool omitFile = false;
      for (unsigned short im=0; im<strmatches.size(); im++){
        for (unsigned short iwp=0; iwp<hnames.size(); iwp++){
          TH2F const* htmp = htmplist.at(im).at(iwp);
          if (!checkGoodHistogram(htmp)){
            MELAerr << "\t- Histogram " << htmp->GetName() << " from " << cinput << " failed integrity checks..." << endl;
            omitFile = true;
          }
        }
      }
      
      // If integrity checks fail, skip the entire file.
      if (omitFile){
        finput->Close();
        foutput->cd();
        continue;
      }

      if (firstFile) MELAout << "\t- First file, so copying histograms..." << endl;
      else MELAout << "\t- Adding to existing histograms..." << endl;
      for (unsigned short im=0; im<strmatches.size(); im++){
        for (unsigned short iwp=0; iwp<hnames.size(); iwp++){
          TH2F* const& htmp = htmplist.at(im).at(iwp);
          foutput->cd();
          if (firstFile) hlist.at(im).at(iwp) = (TH2F*) htmp->Clone(Form("%s_%s_%s", hnames.at(iwp).Data(), strmatches.at(im).Data(), systname.Data()));
          else hlist.at(im).at(iwp)->Add(htmp);
          finput->cd();
        }
      }

      finput->Close();
      firstFile = false;

      foutput->cd();
    }

    // Extract the recursive efficiencies and record them
    for (unsigned short im=0; im<strmatches.size(); im++){
      for (unsigned short iwp=hnames.size()-1; iwp>=1; iwp--){
        hlist.at(im).at(iwp)->Divide(hlist.at(im).at(iwp-1));
        hlist.at(im).at(iwp)->GetZaxis()->SetRangeUser(0, 1);
        foutput->WriteTObject(hlist.at(im).at(iwp));
      }
    }

    // Since the histograms are also owned, delete them in a separate loop
    for (unsigned short im=0; im<strmatches.size(); im++){
      for (unsigned short iwp=1; iwp<hnames.size(); iwp++){
        delete hlist.at(im).at(iwp);
      }
    }
  }

  foutput->Close();
}
