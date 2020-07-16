#include "common_includes.h"
#include "offshell_cutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false
){
  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const coutput_main =
    "output/GJetsMETResolution/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Medium);
  const float btag_medium_thr = BtagHelpers::getBtagWP(false);
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btag_loose_thr = BtagHelpers::getBtagWP(false);

  std::vector<TriggerHelpers::TriggerType> requiredTriggers{ TriggerHelpers::kSinglePho };
  auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;
  bool isData = false;
  bool isQCD = false;
  float pTG_true_range[2]={ -1, -1 }; bool haspTGRange = false;
  for (auto const& sname:sampledirs){
    isData = SampleHelpers::checkSampleIsData(sname);
    if (isData && theGlobalSyst!=sNominal) return;
    if (isData && nchunks>0) return;

    isQCD = sname.Contains("QCD") && sname.Contains("HT");
    if ((sname.Contains("ZGTo2NuG") || sname.Contains("ZGTo2LG")) && sname.Contains("amcatnloFXFX") && !sname.Contains("PtG-130")) pTG_true_range[1]=130;

    break;
  }
  haspTGRange = pTG_true_range[0]!=pTG_true_range[1];
  bool needGenParticleChecks = isQCD || haspTGRange;

  gSystem->mkdir(coutput_main, true);

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;

  PhotonScaleFactorHandler photonSFHandler;
  BtagScaleFactorHandler btagSFHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(needGenParticleChecks);

  bool isFirstInputFile=true;
  for (auto const& sname:sampledirs){
    TString coutput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(coutput, "_MINIAOD", "");

    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sname), "cms3ntuple/SinglePhoton", "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    const int nEntries = sample_tree.getSelectedNEvents();

    double sum_wgts = (isData ? 1.f : 0.f);
    float xsec = 1;
    float xsec_scale = 1;
    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      double sum_wgts_raw_noveto = 0;
      double sum_wgts_raw_withveto = 0;
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
            sum_wgts = 0;
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_raw_withveto += hCounters->GetBinContent(0, 0);
          sum_wgts_raw_noveto += hCounters->GetBinContent(0, 0) / (1. - hCounters->GetBinContent(0, 1));
          ftmp->Close();
        }
        if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters){
        MELAerr << "Please use skim ntuples!" << endl;
        assert(0);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
    }
    MELAout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;
    MELAout << "\t- xsec scale = " << xsec_scale << endl;

    // Set data tracking options
    eventFilter.setTrackDataEvents(isData);
    eventFilter.setCheckUniqueDataEvent(isData && !isFirstInputFile);

    // Configure handlers
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

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree* tout = new TTree("SkimTree", "");
    MELAout << "Created output file " << stroutput << "..." << endl;

#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
    // Event variables
    BRANCH_COMMAND(float, event_wgt);
    BRANCH_COMMAND(float, event_wgt_triggers);
    BRANCH_COMMAND(float, event_wgt_SFs);
    BRANCH_COMMAND(float, pfmet_pTmiss_XY_JER_PartMomShifts); BRANCH_COMMAND(float, pfmet_phimiss_XY_JER_PartMomShifts);
    BRANCH_COMMAND(float, pfmet_pTmiss_XY_JER); BRANCH_COMMAND(float, pfmet_phimiss_XY_JER);
    BRANCH_COMMAND(float, pfmet_pTmiss_XY_PartMomShifts); BRANCH_COMMAND(float, pfmet_phimiss_XY_PartMomShifts);
    BRANCH_COMMAND(float, pfmet_pTmiss_JER_PartMomShifts); BRANCH_COMMAND(float, pfmet_phimiss_JER_PartMomShifts);
    BRANCH_COMMAND(float, pfmet_pTmiss_XY); BRANCH_COMMAND(float, pfmet_phimiss_XY);
    BRANCH_COMMAND(float, pfmet_pTmiss_JER); BRANCH_COMMAND(float, pfmet_phimiss_JER);
    BRANCH_COMMAND(float, pfmet_pTmiss_PartMomShifts); BRANCH_COMMAND(float, pfmet_phimiss_PartMomShifts);
    BRANCH_COMMAND(float, pfmet_pTmiss); BRANCH_COMMAND(float, pfmet_phimiss);
    BRANCH_COMMAND(float, puppimet_pTmiss_PartMomShifts); BRANCH_COMMAND(float, puppimet_phimiss_PartMomShifts);
    BRANCH_COMMAND(float, puppimet_pTmiss); BRANCH_COMMAND(float, puppimet_phimiss);
    BRANCH_COMMAND(unsigned int, event_nvtxs_good);
    BRANCH_COMMAND(unsigned int, event_Njets);
    // Photon
    BRANCH_COMMAND(float, pt_gamma);
    BRANCH_COMMAND(float, eta_gamma);
    BRANCH_COMMAND(float, phi_gamma);
    BRANCH_COMMAND(float, mass_gamma);
    BRANCH_COMMAND(bool, isConversionSafe);
    BRANCH_COMMAND(bool, passHGGSelection);
    BRANCH_COMMAND(bool, isGap);
    BRANCH_COMMAND(bool, isEB);
    BRANCH_COMMAND(bool, isEE);
    BRANCH_COMMAND(bool, isEBEEGap);
    // Jets
    BRANCH_COMMAND(float, pt_jets);
    BRANCH_COMMAND(float, eta_jets);
    BRANCH_COMMAND(float, phi_jets);
    BRANCH_COMMAND(float, mass_jets);
    BRANCH_COMMAND(float, HT_jets);
#undef BRANCH_COMMAND
    float genmet = 0; float genmet_phimiss = 0;
    if (!isData){
      tout->Branch("genmet", &genmet);
      tout->Branch("genmet_phimiss", &genmet_phimiss);
    }

    foutput->cd();

    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events from " << sample_tree.sampleIdentifier << ", starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    size_t n_evts_acc=0;
    size_t n_pass_genWeights=0;
    size_t n_pass_leptonVeto=0;
    size_t n_pass_photons=0;
    size_t n_pass_isotrackVeto=0;
    size_t n_pass_uniqueEvent=0;
    size_t n_pass_commonFilters=0;
    size_t n_pass_goodPVFilter=0;
    size_t n_pass_triggers=0;
    size_t n_pass_HEMfilter=0;
    size_t n_pass_btagVeto=0;
    bool firstEvent=true;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      if (!isData && firstEvent){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
        xsec *= 1000.;
      }
      if (firstEvent) firstEvent=false;

      event_wgt = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts;

      if (!isData){
        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();
        event_wgt *= genInfo->getGenWeight(true);
        genmet = genInfo->extras.genmet_met;
        genmet_phimiss = genInfo->extras.genmet_metPhi;
        auto const& genparticles = genInfoHandler.getGenParticles();

        if (needGenParticleChecks){
          if (isQCD){
            auto const& genparticles = genInfoHandler.getGenParticles();
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState && part->pt()>=25.f){
                event_wgt = 0;
                break;
              }
            }
          }
          if (haspTGRange){
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
                if ((pTG_true_range[0]>=0.f && part->pt()<pTG_true_range[0]) || (pTG_true_range[1]>=0.f && part->pt()>=pTG_true_range[1])) event_wgt = 0;
                break;
              }
            }
          }
        }

        simEventHandler.constructSimEvent(theGlobalSyst);
        event_wgt *= simEventHandler.getPileUpWeight()*simEventHandler.getL1PrefiringWeight();

        if (event_wgt==0.f) continue;
      }
      n_pass_genWeights++;

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
      if (n_muons_veto!=0) continue;

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
      if (n_electrons_veto!=0) continue;
      n_pass_leptonVeto++;

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
      if (n_photons_tight!=1) continue;
      pt_gamma = theChosenPhoton->pt();
      eta_gamma = theChosenPhoton->eta();
      phi_gamma = theChosenPhoton->phi();
      mass_gamma = theChosenPhoton->m();
      isConversionSafe = theChosenPhoton->testSelection(PhotonSelectionHelpers::kConversionSafe);
      passHGGSelection = theChosenPhoton->extras.id_cutBased_HGG_Bits;
      isGap = theChosenPhoton->isGap();
      isEBEEGap = theChosenPhoton->isEBEEGap();
      isEB = theChosenPhoton->isEB();
      isEE = theChosenPhoton->isEE();
      if ((theGlobalSyst!=SystematicsHelpers::sNominal && pt_gamma<150.f) || pt_gamma<50.f) continue; // Skip photons that are more likely to be pi^0s
      n_pass_photons++;

      event_wgt_SFs = SF_muons*SF_electrons*SF_photons;

      isotrackHandler.constructIsotracks(&muons, &electrons);
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotrackHandler.getProducts()){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;
      n_pass_isotrackVeto++;

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      auto const& puppimet = jetHandler.getPFPUPPIMET();

#define PTMISS_FILL_COMMAND(SUFFIX, OPT_XY, OPT_JER, OPT_PARTMOMSHIFTS) \
      auto pfmet_p4_##SUFFIX = pfmet->p4(OPT_XY, OPT_JER, OPT_PARTMOMSHIFTS); \
      pfmet_pTmiss_##SUFFIX = pfmet_p4_##SUFFIX.Pt(); pfmet_phimiss_##SUFFIX = pfmet_p4_##SUFFIX.Phi();
      PTMISS_FILL_COMMAND(XY_JER_PartMomShifts, true, true, true);
      PTMISS_FILL_COMMAND(XY_JER, true, true, false);
      PTMISS_FILL_COMMAND(XY_PartMomShifts, true, false, true);
      PTMISS_FILL_COMMAND(JER_PartMomShifts, false, true, true);
      PTMISS_FILL_COMMAND(XY, true, false, false);
      PTMISS_FILL_COMMAND(JER, false, true, false);
      PTMISS_FILL_COMMAND(PartMomShifts, false, false, true);
#undef PTMISS_FILL_COMMAND
      auto pfmet_p4 = pfmet->p4(false, false, false);
      pfmet_pTmiss = pfmet_p4.Pt();
      pfmet_phimiss = pfmet_p4.Phi();

      auto puppimet_p4_PartMomShifts = puppimet->p4(false, false, true);
      puppimet_pTmiss_PartMomShifts = puppimet_p4_PartMomShifts.Pt();
      puppimet_phimiss_PartMomShifts = puppimet_p4_PartMomShifts.Phi();
      auto puppimet_p4 = puppimet->p4(false, false, false);
      puppimet_pTmiss = puppimet_p4.Pt();
      puppimet_phimiss = puppimet_p4.Phi();

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      n_pass_uniqueEvent++;

      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters() || !eventFilter.hasGoodVertex()) continue;
      n_pass_commonFilters++;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;
      n_pass_goodPVFilter++;

      event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList, nullptr, nullptr, &photons, nullptr, nullptr, nullptr);
      if (event_wgt_triggers == 0.f) continue;
      n_pass_triggers++;

      // Test HEM filter
      if (!eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) continue;
      n_pass_HEMfilter++;

      HT_jets = 0;
      ParticleObject::LorentzVector_t ak4jets_sump4(0, 0, 0, 0);
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      unsigned int n_ak4jets_tight_btagged = 0;
      float SF_btagging = 1;
      for (auto* jet:ak4jets){
        float theSF = 1;
        if (!isData) btagSFHandler.getSFAndEff(theGlobalSyst, jet, theSF, nullptr);
        if (theSF != 0.f) SF_btagging *= theSF;

        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btag_loose_thr) n_ak4jets_tight_btagged++;

          HT_jets += jet->pt();
          ak4jets_sump4 += jet->p4();
        }
      }
      if (n_ak4jets_tight_btagged>0) continue;
      event_wgt_SFs *= SF_btagging;
      n_pass_btagVeto++;

      event_Njets = ak4jets_tight.size();
      pt_jets = ak4jets_sump4.Pt();
      eta_jets = ak4jets_sump4.Eta();
      phi_jets = ak4jets_sump4.Phi();
      mass_jets = ak4jets_sump4.M();

      sample_tree.getVal("vtxs_nvtxs_good", event_nvtxs_good);

      tout->Fill();
      n_evts_acc++;
    } // End loop over events

    MELAout << "Number of events accepted from " << sample_tree.sampleIdentifier << ": " << n_evts_acc << " / " << (ev_end - ev_start) << endl;
    MELAout << "\t- Number of events passing each cut:\n"
      << "\t\t- Gen. weights!=0: " << n_pass_genWeights << '\n'
      << "\t\t- Lepton veto: " <<  n_pass_leptonVeto << '\n'
      << "\t\t- Photon selection: " <<  n_pass_photons << '\n'
      << "\t\t- Isotrack veto: " <<  n_pass_isotrackVeto << '\n'
      << "\t\t- Unique event: " <<  n_pass_uniqueEvent << '\n'
      << "\t\t- Common filters: " <<  n_pass_commonFilters << '\n'
      << "\t\t- Good PV filter: " << n_pass_goodPVFilter << '\n'
      << "\t\t- Trigger: " <<  n_pass_triggers << '\n'
      << "\t\t- HEM15/16 veto: " << n_pass_HEMfilter << '\n'
      << "\t\t- b-tag veto: " <<  n_pass_btagVeto
      << endl;

    // Set this flag for data so that subsequent files ensure checking for unique events
    isFirstInputFile=false;

    // Write output
    foutput->WriteTObject(tout);
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples list
}
