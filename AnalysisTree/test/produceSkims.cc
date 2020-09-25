#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


void produceSkims(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate="",
  int ichunk=0, int nchunks=0,
  bool doDilepton=true, bool doDilepton_Control=true, bool doSingleLepton=true, bool doGJets=true
){
  enum FinalStateType{
    kDilepton,
    kDilepton_Control,
    kSingleLepton,
    kGJets,
    nFinalStateTypes
  };

  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop:"+prodVersion);

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton{
    TriggerHelpers::kTripleLep,
    TriggerHelpers::kDoubleMu, TriggerHelpers::kDoubleMu_Prescaled,
    TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt, TriggerHelpers::kSingleEle_L1EG,
    TriggerHelpers::kMuEle, TriggerHelpers::kMuEle_Extra
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleLepton{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt, TriggerHelpers::kSingleMu_Prescaled, TriggerHelpers::kSingleMu_Control,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt, TriggerHelpers::kSingleEle_Prescaled, TriggerHelpers::kSingleEle_Control
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_GJets{
    TriggerHelpers::kSinglePho
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_OrtogonalControl{
    TriggerHelpers::kAK8PFJet_Control,
    TriggerHelpers::kPFHT_Control,
    TriggerHelpers::kPFMET_MHT_Control
  };
  std::vector<std::string> triggerCheckList_Dilepton = TriggerHelpers::getHLTMenus(requiredTriggers_Dilepton);
  std::vector<std::string> triggerCheckList_SingleLepton = TriggerHelpers::getHLTMenus(requiredTriggers_SingleLepton);
  std::vector<std::string> triggerCheckList_GJets = TriggerHelpers::getHLTMenus(requiredTriggers_GJets);
  std::vector<std::string> triggerCheckList_OrtogonalControl = TriggerHelpers::getHLTMenus(requiredTriggers_OrtogonalControl);

  std::vector<TString> validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t nValidDataPeriods = validDataPeriods.size();

  // Get handlers
  // Some of these are needed only to enable the branches
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  PFCandidateHandler pfcandidateHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  SuperclusterHandler superclusterHandler;
  FSRHandler fsrHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;

  OverlapMapHandler<MuonObject, AK4JetObject> overlapMap_muons_ak4jets;
  OverlapMapHandler<MuonObject, AK8JetObject> overlapMap_muons_ak8jets;
  OverlapMapHandler<ElectronObject, AK4JetObject> overlapMap_electrons_ak4jets;
  OverlapMapHandler<ElectronObject, AK8JetObject> overlapMap_electrons_ak8jets;
  OverlapMapHandler<PhotonObject, AK4JetObject> overlapMap_photons_ak4jets;
  OverlapMapHandler<PhotonObject, AK8JetObject> overlapMap_photons_ak8jets;
  // No need to register jet overlaps because we do not reconstruct jets in the skimmer.
  /*
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
  */

  genInfoHandler.setAllowLargeGenWeightRemoval(true);

  eventFilter.setTrackDataEvents(false);
  eventFilter.setCheckUniqueDataEvent(false);
  eventFilter.setCheckHLTPathRunRanges(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  TString const stroutputcore = Form("output/Skims/%s", strdate.Data());

  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    //if (!isData && nchunks>0) return;

    std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts(1, SystematicsHelpers::sNominal);
    if (!isData){
      using namespace SystematicsHelpers;
      allowedSysts = std::vector<SystematicsHelpers::SystematicVariationTypes>{
        sNominal,

        eEleScaleDn, eEleScaleUp,
        eEleResDn, eEleResUp,

        eMuScaleDn, eMuScaleUp,
        eMuResDn, eMuResUp,

        ePhoScaleDn, ePhoScaleUp,
        ePhoResDn, ePhoResUp,
      };
    }

    TString const cinputcore = SampleHelpers::getDatasetDirectoryName(strSample);
    TString const cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    std::vector<TString> branchnames;
    sample_tree.getValidBranchNamesWithoutAlias(branchnames, false);

    const int nEntries = sample_tree.getNEvents();
    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    unsigned int n_zero_genwgts=0;
    double frac_zero_genwgts=0;
    double sum_wgts_noPU=0;
    double sum_abs_wgts_noPU=0;
    std::vector<double> sum_wgts(nValidDataPeriods+1, 0);
    std::vector<double> sum_wgts_PUDn(nValidDataPeriods+1, 0);
    std::vector<double> sum_wgts_PUUp(nValidDataPeriods+1, 0);
    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(false);
      genInfoHandler.setAcquireGenParticles(false);
      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      MELAout << "Initial MC loop over " << ev_end-ev_start << " / " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
      for (int ev=ev_start; ev<ev_end; ev++){
        HelperFunctions::progressbar(ev, nEntries);
        sample_tree.getEvent(ev);

        genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
        auto const& genInfo = genInfoHandler.getGenInfo();
        double genwgt = genInfo->getGenWeight(true);
        if (genwgt==0.){
          n_zero_genwgts++;
          continue;
        }

        sum_wgts_noPU += genwgt;
        sum_abs_wgts_noPU += std::abs(genwgt);
        for (size_t idp=0; idp<nValidDataPeriods+1; idp++){
          if (idp==0) SampleHelpers::setDataPeriod(period);
          else SampleHelpers::setDataPeriod(validDataPeriods.at(idp-1));

          double puwgt;
          simEventHandler.constructSimEvent(SystematicsHelpers::sNominal);
          puwgt = simEventHandler.getPileUpWeight();
          sum_wgts.at(idp) += genwgt * puwgt;
          simEventHandler.resetCache();
          simEventHandler.constructSimEvent(SystematicsHelpers::ePUDn);
          puwgt = simEventHandler.getPileUpWeight();
          sum_wgts_PUDn.at(idp) += genwgt * puwgt;
          simEventHandler.resetCache();
          simEventHandler.constructSimEvent(SystematicsHelpers::ePUUp);
          puwgt = simEventHandler.getPileUpWeight();
          sum_wgts_PUUp.at(idp) += genwgt * puwgt;
          simEventHandler.resetCache();
        }
        SampleHelpers::setDataPeriod(period);
      }
      if (ev_end-ev_start>0) frac_zero_genwgts = double(n_zero_genwgts)/double(ev_end-ev_start);
      {
        bool has_lheMEweights = false;
        bool has_lheparticles = false;
        bool has_genparticles = false;
        std::vector<TString> allbranchnames; sample_tree.getValidBranchNamesWithoutAlias(allbranchnames, false);
        for (auto const& bname:allbranchnames){
          if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights=true;
          else if (bname.Contains(GenInfoHandler::colName_lheparticles)) has_lheparticles = true;
          else if (bname.Contains(GenInfoHandler::colName_genparticles)) has_genparticles = true;
        }
        genInfoHandler.setAcquireLHEMEWeights(has_lheMEweights);
        genInfoHandler.setAcquireLHEParticles(has_lheparticles);
        genInfoHandler.setAcquireGenParticles(has_genparticles);
      }
      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);
    }

    pfcandidateHandler.bookBranches(&sample_tree);
    pfcandidateHandler.wrapTree(&sample_tree);

    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    bool hasSuperclusters = false;
    for (auto const& bname:branchnames){
      if (bname.BeginsWith(SuperclusterHandler::colName.data())){
        hasSuperclusters = true;
        break;
      }
    }
    if (hasSuperclusters){
      superclusterHandler.bookBranches(&sample_tree);
      superclusterHandler.wrapTree(&sample_tree);
    }

    fsrHandler.bookBranches(&sample_tree);
    fsrHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    isotrackHandler.bookBranches(&sample_tree);
    isotrackHandler.wrapTree(&sample_tree);

    vertexHandler.bookBranches(&sample_tree);
    vertexHandler.wrapTree(&sample_tree);

    overlapMap_muons_ak4jets.bookBranches(&sample_tree); overlapMap_muons_ak4jets.wrapTree(&sample_tree);
    overlapMap_muons_ak8jets.bookBranches(&sample_tree); overlapMap_muons_ak8jets.wrapTree(&sample_tree);

    overlapMap_electrons_ak4jets.bookBranches(&sample_tree); overlapMap_electrons_ak4jets.wrapTree(&sample_tree);
    overlapMap_electrons_ak8jets.bookBranches(&sample_tree); overlapMap_electrons_ak8jets.wrapTree(&sample_tree);

    overlapMap_photons_ak4jets.bookBranches(&sample_tree); overlapMap_photons_ak4jets.wrapTree(&sample_tree);
    overlapMap_photons_ak8jets.bookBranches(&sample_tree); overlapMap_photons_ak8jets.wrapTree(&sample_tree);

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
    MELAout << "Creating output file " << stroutput << "..." << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TDirectory* outdir = foutput->mkdir("cms3ntuple");
    outdir->cd();
    if (!isData){
      MELAout << "Writing the sum of gen. weights:" << endl;
      MELAout << "\t- Fraction of discarded events with null weight: " << frac_zero_genwgts << endl;
      MELAout << "\t- Before PU reweighting: " << setprecision(15) << sum_wgts_noPU << endl;
      MELAout << "\t- Fraction of negative weights before PU reweighting: " << setprecision(15) << (sum_abs_wgts_noPU-sum_wgts_noPU)*0.5/sum_abs_wgts_noPU << endl;
      MELAout << "\t- PU nominal: " << setprecision(15) << sum_wgts << endl;
      MELAout << "\t- PU down: " << setprecision(15) << sum_wgts_PUDn << endl;
      MELAout << "\t- PU up: " << setprecision(15) << sum_wgts_PUUp << endl;

      TH2D hCounters("Counters", "", 3, 0, 3, nValidDataPeriods+1, 0, nValidDataPeriods+1);
      hCounters.SetBinContent(0, 0, sum_wgts_noPU); // Sum with no PU reweighting
      hCounters.SetBinContent(0, 1, frac_zero_genwgts); // Fraction of discarded events
      hCounters.SetBinContent(0, 2, (ev_end-ev_start)); // Total number of events iterated
      for (size_t idp=0; idp<nValidDataPeriods+1; idp++){
        hCounters.SetBinContent(1, idp+1, sum_wgts.at(idp));
        hCounters.SetBinContent(2, idp+1, sum_wgts_PUDn.at(idp));
        hCounters.SetBinContent(3, idp+1, sum_wgts_PUUp.at(idp));
      }
      outdir->WriteTObject(&hCounters);
    }
    //sample_tree.unmuteAllBranches();

    std::vector<TTree*> outtree(nFinalStateTypes, nullptr);
    if (doDilepton){
      outtree[kDilepton] = sample_tree.getSelectedTree()->CloneTree(0);
      outtree[kDilepton]->SetName("Dilepton");
      outtree[kDilepton]->SetAutoSave(0);
    }
    if (doDilepton_Control){
      outtree[kDilepton_Control] = sample_tree.getSelectedTree()->CloneTree(0);
      outtree[kDilepton_Control]->SetName("Dilepton_Control");
      outtree[kDilepton_Control]->SetAutoSave(0);
    }
    if (doSingleLepton){
      outtree[kSingleLepton] = sample_tree.getSelectedTree()->CloneTree(0);
      outtree[kSingleLepton]->SetName("SingleLepton");
      outtree[kSingleLepton]->SetAutoSave(0);
    }
    if (doGJets){
      outtree[kGJets] = sample_tree.getSelectedTree()->CloneTree(0);
      outtree[kGJets]->SetName("SinglePhoton");
      outtree[kGJets]->SetAutoSave(0);
    }

    // Loop over the tree
    std::vector<size_t> n_acc(nFinalStateTypes, 0);
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);
      if (ev%10000==0) MELAout << sample_tree.sampleIdentifier << " events: " << n_acc << " / " << ev << " / " << nEntries << endl;

      if (!isData){
        genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
        auto const& genInfo = genInfoHandler.getGenInfo();
        double genwgt = genInfo->getGenWeight(true);
        if (genwgt==0.) continue;
      }

      eventFilter.constructFilters(nullptr);
      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters(EventFilterHandler::kMETFilters_Standard)) continue;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;

      float trigwgt_dilepton = eventFilter.getTriggerWeight(triggerCheckList_Dilepton);
      float trigwgt_singlelepton = eventFilter.getTriggerWeight(triggerCheckList_SingleLepton);
      float trigwgt_gjets = eventFilter.getTriggerWeight(triggerCheckList_GJets);
      float trigwgt_orthogonalcontrol = eventFilter.getTriggerWeight(triggerCheckList_OrtogonalControl);
      trigwgt_singlelepton += trigwgt_orthogonalcontrol;
      trigwgt_gjets += trigwgt_orthogonalcontrol;
      trigwgt_dilepton += trigwgt_singlelepton;
      if ((trigwgt_dilepton + trigwgt_singlelepton + trigwgt_gjets)==0.f) continue;

      // Categorize events based on how many leptons or photons they have
      // Do not apply any b-jet or HEM vetoes because their application is analysis-specific

      bool doRecordTree[nFinalStateTypes]={ 0 };
      for (auto const& syst:allowedSysts){
        pfcandidateHandler.constructPFCandidates(syst);
        auto const& pfcandidates = pfcandidateHandler.getProducts();

        muonHandler.constructMuons(syst, &pfcandidates);
        electronHandler.constructElectrons(syst, &pfcandidates);
        photonHandler.constructPhotons(syst, &pfcandidates);

        // Particle counts before disambiguation
        size_t n_muons_tight_predisambiguation = 0;
        size_t n_muons_tnp_predisambiguation = 0;
        size_t n_muons_fakeable_predisambiguation = 0;
        size_t n_electrons_tight_predisambiguation = 0;
        size_t n_electrons_tnp_predisambiguation = 0;
        size_t n_electrons_fakeable_predisambiguation = 0;
        size_t n_photons_tight_predisambiguation = 0;
        {
          auto const& muons = muonHandler.getProducts();
          for (auto const& part:muons){
            if (ParticleSelectionHelpers::isTightParticle(part)) n_muons_tight_predisambiguation++;
            else{
              if (part->testSelectionBit(MuonSelectionHelpers::kProbeId)) n_muons_tnp_predisambiguation++;
              if (part->testSelectionBit(MuonSelectionHelpers::kFakeableBase)) n_muons_fakeable_predisambiguation++;
            }
          }
          auto const& electrons = electronHandler.getProducts();
          for (auto const& part:electrons){
            if (ParticleSelectionHelpers::isTightParticle(part)) n_electrons_tight_predisambiguation++;
            else{
              if (part->testSelectionBit(ElectronSelectionHelpers::kProbeId)) n_electrons_tnp_predisambiguation++;
              if (part->testSelectionBit(ElectronSelectionHelpers::kFakeableBase)) n_electrons_fakeable_predisambiguation++;
            }
          }
          auto const& photons = photonHandler.getProducts();
          for (auto const& part:photons){
            if (ParticleSelectionHelpers::isTightParticle(part)) n_photons_tight_predisambiguation++;
          }
        }
        size_t n_leptons_tight_predisambiguation = n_muons_tight_predisambiguation + n_electrons_tight_predisambiguation;
        size_t n_leptons_fakeable_predisambiguation = n_muons_fakeable_predisambiguation + n_electrons_fakeable_predisambiguation;

        // Particle counts after disambiguation
        size_t n_muons_tight = 0;
        size_t n_muons_tnp = 0;
        size_t n_muons_fakeable = 0;
        size_t n_electrons_tight = 0;
        size_t n_electrons_tnp = 0;
        size_t n_electrons_fakeable = 0;
        size_t n_photons_tight = 0;
        {
          particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);
          auto const& muons = muonHandler.getProducts();
          for (auto const& part:muons){
            if (ParticleSelectionHelpers::isTightParticle(part)) n_muons_tight++;
            else{
              if (part->testSelectionBit(MuonSelectionHelpers::kProbeId)) n_muons_tnp++;
              if (part->testSelectionBit(MuonSelectionHelpers::kFakeableBase)) n_muons_fakeable++;
            }
          }
          auto const& electrons = electronHandler.getProducts();
          for (auto const& part:electrons){
            if (ParticleSelectionHelpers::isTightParticle(part)) n_electrons_tight++;
            else{
              if (part->testSelectionBit(ElectronSelectionHelpers::kProbeId)) n_electrons_tnp++;
              if (part->testSelectionBit(ElectronSelectionHelpers::kFakeableBase)) n_electrons_fakeable++;
            }
          }
          auto const& photons = photonHandler.getProducts();
          for (auto const& part:photons){
            if (ParticleSelectionHelpers::isTightParticle(part)) n_photons_tight++;
          }
        }
        size_t n_leptons_tight = n_muons_tight + n_electrons_tight;
        size_t n_leptons_fakeable = n_muons_fakeable + n_electrons_fakeable;

        doRecordTree[kDilepton] |= (n_leptons_tight>=2 || n_leptons_tight_predisambiguation>=2);
        doRecordTree[kDilepton_Control] |= (
          (n_muons_tight==1 && n_muons_tnp==1) || (n_muons_tight_predisambiguation==1 && n_muons_tnp_predisambiguation==1)
          ||
          (n_electrons_tight==1 && n_electrons_tnp==1) || (n_electrons_tight_predisambiguation==1 && n_electrons_tnp_predisambiguation==1)
          ||
          (n_muons_tight==1 && n_muons_fakeable==1) || (n_muons_tight_predisambiguation==1 && n_muons_fakeable_predisambiguation==1)
          ||
          (n_electrons_tight==1 && n_electrons_fakeable==1) || (n_electrons_tight_predisambiguation==1 && n_electrons_fakeable_predisambiguation==1)
          );
        doRecordTree[kSingleLepton] |= (
          (n_muons_tight==1 && n_electrons_tight<2) || (n_muons_tight_predisambiguation==1 && n_electrons_tight_predisambiguation<2)
          ||
          (n_muons_tight<2 && n_electrons_tight==1) || (n_muons_tight_predisambiguation<2 && n_electrons_tight_predisambiguation==1)
          ||
          (n_leptons_tight==0 && (n_muons_fakeable==1 || n_electrons_fakeable==1))
          ||
          (n_leptons_tight_predisambiguation==0 && (n_muons_fakeable_predisambiguation==1 || n_electrons_fakeable_predisambiguation==1))
          );
        doRecordTree[kGJets] |= (
          (n_photons_tight==1 && n_leptons_tight==0)
          ||
          (n_photons_tight_predisambiguation==1 && n_leptons_tight_predisambiguation==0)
          );

        // Reset the caches because the next systematic should re-reconstruct the particles from scratch
        pfcandidateHandler.resetCache();
        muonHandler.resetCache();
        electronHandler.resetCache();
        isotrackHandler.resetCache();
        photonHandler.resetCache();
        fsrHandler.resetCache();
        jetHandler.resetCache();
        overlapMap_muons_ak4jets.resetCache();
        overlapMap_muons_ak8jets.resetCache();
        overlapMap_electrons_ak4jets.resetCache();
        overlapMap_electrons_ak8jets.resetCache();
        overlapMap_photons_ak4jets.resetCache();
        overlapMap_photons_ak8jets.resetCache();
      }

      // Combine triggers after object combinations
      doRecordTree[kDilepton] &= (outtree[kDilepton] && trigwgt_dilepton!=0.f);
      doRecordTree[kDilepton_Control] &= (outtree[kDilepton_Control] && trigwgt_dilepton!=0.f);
      doRecordTree[kSingleLepton] &= (outtree[kSingleLepton] && trigwgt_singlelepton!=0.f);
      doRecordTree[kGJets] &= (outtree[kGJets] && trigwgt_gjets!=0.f);
      if (doRecordTree[kDilepton]){
        outtree[kDilepton]->Fill();
        n_acc[kDilepton]++;
      }
      else if (doRecordTree[kDilepton_Control]){
        outtree[kDilepton_Control]->Fill();
        n_acc[kDilepton_Control]++;
      }
      else if (doRecordTree[kSingleLepton]){
        outtree[kSingleLepton]->Fill();
        n_acc[kSingleLepton]++;
      }

      if (doRecordTree[kGJets]){
        outtree[kGJets]->Fill();
        n_acc[kGJets]++;
      }
    }
    MELAout << sample_tree.sampleIdentifier << " total number of accumulated events: "
      << "(Dilepton, dilepton control, single lepton, single photon) = ( " << n_acc << " )"
      << " / " << (ev_end - ev_start) << " / " << nEntries << endl;

    for (unsigned int i=0; i<(unsigned int) nFinalStateTypes; i++){ outdir->WriteTObject(outtree.at(i)); delete outtree.at(i); }
    foutput->cd();
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  }
}
