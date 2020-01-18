#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


void produceBtaggingEfficiencies(TString strSampleSet, TString period, TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  constexpr SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure(period, "191212");

  std::vector<BtagHelpers::BtagWPType> btagwptypes{
    //BtagHelpers::kDeepCSV_Loose,
    //BtagHelpers::kDeepCSV_Medium,
    //BtagHelpers::kDeepCSV_Tight,
    BtagHelpers::kDeepFlav_Loose,
    BtagHelpers::kDeepFlav_Medium,
    BtagHelpers::kDeepFlav_Tight
  };
  std::unordered_map<BtagHelpers::BtagWPType, float> btagwp_type_val_map;
  for (auto const& btagwptype:btagwptypes){
    BtagHelpers::setBtagWPType(btagwptype);
    btagwp_type_val_map[btagwptype] = BtagHelpers::getBtagWP(false);
  }

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(
    {
      OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kMuEle,
      OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle
    }
  );
  std::vector<std::string> triggerCheckList_SinglePhoton = OffshellTriggerHelpers::getHLTMenus(OffshellTriggerHelpers::kSinglePho);

  // Get handlers
  GenInfoHandler genInfoHandler;
  SimEventHandler simEventHandler;
  EventFilterHandler eventFilter;
  VertexHandler vertexHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  DileptonHandler dileptonHandler;

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

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    float sum_wgts=0;
    const int nEntries = sample_tree.getSelectedNEvents();

    simEventHandler.bookBranches(&sample_tree);
    simEventHandler.wrapTree(&sample_tree);

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    sample_tree.silenceUnused();

    MELAout << "Initial MC loop to determine sample normalization:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
      auto const& genInfo = genInfoHandler.getGenInfo();
      float genwgt = genInfo->getGenWeight(true);

      simEventHandler.constructSimEvent(theGlobalSyst);
      float puwgt = simEventHandler.getPileUpWeight();
      float wgt = genwgt * puwgt;
      sum_wgts += wgt;
    }
    float wgtNorm = static_cast<const float>(nEntries) / sum_wgts;
    MELAout << "Sum of weights is " << sum_wgts << " over the " << nEntries << " entries. The weight normalization becomes " << wgtNorm << "." << endl;

    vertexHandler.bookBranches(&sample_tree);
    vertexHandler.wrapTree(&sample_tree);

    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    sample_tree.silenceUnused();

    MELAout << "Completed getting the rest of the handles..." << endl;

    // Create output tree
    TString cSample = SampleHelpers::getDatasetDirectoryCoreName(strSample.Data()).data();
    HelperFunctions::replaceString(cSample, "/", "_");
    TString stroutput = coutput_main + "/" + cSample + "_btageff_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    ExtendedBinning ptbins({ 0, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 3000 });
    ExtendedBinning etabins(24, -2.4, 2.4);
    etabins.addBinBoundary(-2.5); etabins.addBinBoundary(2.5);
    etabins.addBinBoundary(-5.); etabins.addBinBoundary(5.);
    std::vector<TH2F> h_All{
      TH2F("AllJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
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
    MELAout << "Starting to loop over " << nEntries << " events" << endl;
    unsigned int n_acc=0;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      genInfoHandler.constructGenInfo(theGlobalSyst);
      auto const& genInfo = genInfoHandler.getGenInfo();
      float genwgt = genInfo->getGenWeight(true);

      simEventHandler.constructSimEvent(theGlobalSyst);
      float puwgt = simEventHandler.getPileUpWeight();

      float wgt = genwgt * puwgt * wgtNorm;

      eventFilter.constructFilters();
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float trigwgt = eventFilter.getTriggerWeight(triggerCheckList);
      float trigwgt_singlephoton = eventFilter.getTriggerWeight(triggerCheckList_SinglePhoton);

      bool doSinglePhotonHypothesis = (trigwgt_singlephoton>0.f);
      bool doDileptonHypothesis = (trigwgt>0.f);
      if (!doDileptonHypothesis && !doSinglePhotonHypothesis) continue;

      muonHandler.constructMuons(theGlobalSyst);
      auto const& muons = muonHandler.getProducts();

      electronHandler.constructElectrons(theGlobalSyst);
      auto const& electrons = electronHandler.getProducts();

      photonHandler.constructPhotons(theGlobalSyst);
      auto const& photons = photonHandler.getProducts();

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& pfmet = jetHandler.getPFMET();

      if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons, &photons, &ak4jets, &ak8jets)) continue;

      float abs_dPhi_min_pTj_puppimet_pTmiss = TMath::Pi();
      ParticleObject::LorentzVector_t p4_alljets;
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);

          p4_alljets = p4_alljets + jet->p4();

          float dphi_tmp;
          HelperFunctions::deltaPhi(float(jet->phi()), float(puppimet->phi()), dphi_tmp);
          abs_dPhi_min_pTj_puppimet_pTmiss = std::min(abs_dPhi_min_pTj_puppimet_pTmiss, std::abs(dphi_tmp));
        }
      }
      size_t n_ak4jets_tight = ak4jets_tight.size();
      if (n_ak4jets_tight==0) continue;

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();

      DileptonObject* theChosenDilepton = nullptr;
      size_t nTightDilep = 0;
      for (auto const& dilepton:dileptons){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          nTightDilep++;
        }
      }

      size_t n_muons_veto = 0; for (auto const& part:muons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++; }
      size_t n_electrons_veto = 0; for (auto const& part:electrons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++; }
      size_t n_photons_veto = 0;
      size_t n_photons_tight = 0;
      PhotonObject* theChosenPhoton = nullptr;
      for (auto const& part:photons){
        if (ParticleSelectionHelpers::isVetoParticle(part)) n_photons_veto++;
        if (ParticleSelectionHelpers::isTightParticle(part)){
          if (!theChosenPhoton) theChosenPhoton = part;
          n_photons_tight++;
        }
      }

      if (doDileptonHypothesis && !doSinglePhotonHypothesis){
        if (!(nTightDilep==1 && n_photons_veto==0)) continue;
      }
      else if (!doDileptonHypothesis && doSinglePhotonHypothesis){
        if (!((n_muons_veto+n_electrons_veto)==0 && n_photons_tight==1)) continue;
      }
      else if (doDileptonHypothesis && doSinglePhotonHypothesis){
        if (
          !((nTightDilep==1 && n_photons_veto==0) ^ ((n_muons_veto+n_electrons_veto)==0 && n_photons_tight==1))
          ) continue;
        else if (nTightDilep==1 && n_photons_veto==0) doSinglePhotonHypothesis=false;
        else if ((n_muons_veto+n_electrons_veto)==0 && n_photons_tight==1) doDileptonHypothesis=false;
      }
      float overallTriggerWeight = (doDileptonHypothesis ? trigwgt : trigwgt_singlephoton);
      wgt *= overallTriggerWeight;

      constexpr float met_thr = 85.;
      bool pass_puppimet_thr = puppimet->pt()>=met_thr;
      bool pass_abs_dPhi_min_pTj_puppimet_pTmiss = abs_dPhi_min_pTj_puppimet_pTmiss>=0.6;
      if (!(pass_puppimet_thr && pass_abs_dPhi_min_pTj_puppimet_pTmiss)) continue;
      if (theChosenDilepton){
        /*
        bool is_ee=false, is_mumu=false, is_emu=false;
        if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -121) is_ee=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -143) is_emu=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -169) is_mumu=true;
        */
        ParticleObject* leadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(0) : theChosenDilepton->daughter(1));
        ParticleObject* subleadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(1) : theChosenDilepton->daughter(0));
        float pTl1 = leadingLepton->pt();
        float pTl2 = subleadingLepton->pt();
        float pTll = theChosenDilepton->pt();
        float mll = theChosenDilepton->m();

        float dPhi_pTlljets_puppimet_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(puppimet->phi()), dPhi_pTlljets_puppimet_pTmiss);
        float pTlljets = (theChosenDilepton->p4() + p4_alljets).Pt();

        bool pass_pTl1 = pTl1>=25.;
        bool pass_pTl2 = pTl2>=20.;
        bool pass_mll = mll>=81.;
        bool pass_pTll = pTll>=35.;
        bool pass_puppimet_over_pTlljets_thr = (puppimet->pt() / pTlljets)>=std::pow(met_thr / pTlljets, 1.5);
        bool pass_puppimet_dPhilljets_thr = std::abs(dPhi_pTlljets_puppimet_pTmiss)>=2.6;
        if (
          !(
            pass_pTll
            &&
            pass_puppimet_over_pTlljets_thr && pass_puppimet_dPhilljets_thr
            &&
            pass_pTl1 && pass_pTl2 && pass_mll
            )
          ) continue;
      }
      if (theChosenPhoton){
        float pTG = theChosenPhoton->pt();

        float dPhi_pTGjets_puppimet_pTmiss; HelperFunctions::deltaPhi(float((theChosenPhoton->p4() + p4_alljets).Phi()), float(puppimet->phi()), dPhi_pTGjets_puppimet_pTmiss);
        float pTGjets = (theChosenPhoton->p4() + p4_alljets).Pt();

        bool pass_pTG = pTG>=35.;
        bool pass_puppimet_over_pTGjets_thr = (puppimet->pt() / pTGjets)>=std::pow(met_thr / pTGjets, 1.5);
        bool pass_puppimet_dPhiGjets_thr = std::abs(dPhi_pTGjets_puppimet_pTmiss)>=2.6;
        if (
          !(
            pass_pTG
            &&
            pass_puppimet_over_pTGjets_thr && pass_puppimet_dPhiGjets_thr
            )
          ) continue;
      }

      for (auto const& jet:ak4jets_tight){
        int const& jetFlavor = jet->extras.hadronFlavour;
        unsigned int iflav = 2;
        if (abs(jetFlavor)==5) iflav = 0;
        else if (abs(jetFlavor)==4) iflav = 1;

        float jetpt_onlyjec = jet->pt() * (jet->extras.JECNominal / jet->currentSystScale);

        h_All.at(iflav).Fill(jetpt_onlyjec, jet->eta(), wgt);
        for (auto const& btagwptype:btagwptypes){
          BtagHelpers::setBtagWPType(btagwptype);

          std::vector<TH2F>* hlist = nullptr;
          switch (btagwptype){
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

          if (jet->getBtagValue()>=btagwp_type_val_map[btagwptype]) hlist->at(iflav).Fill(jetpt_onlyjec, jet->eta(), wgt);
        }
      }

    }

    for (unsigned int iflav=0; iflav<3; iflav++){
      foutput->WriteTObject(&h_All.at(iflav));
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
