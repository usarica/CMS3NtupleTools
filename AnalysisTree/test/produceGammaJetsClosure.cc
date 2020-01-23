#include <cassert>
#include <limits>
#include "common_includes.h"
#include "RooMsgService.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooVoigtian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TNumericUtil.hh"
//#include <HiggsAnalysis/CombinedLimit/interface/AsymPow.h>
#include "TMatrixDSym.h"
#include "TChain.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLegend.h"


using namespace std;
using namespace RooFit;


void produceGammaJetsClosure(TString strSampleSet, TString period, unsigned int istep, TString strdate=""){
  if (istep>2) return;

  constexpr bool doZZselections = true;

  gStyle->SetOptStat(0);
  TDirectory* curdir = gDirectory;

  bool isGammaJetsSample = false;
  {
    TString strSampleSetLower = strSampleSet; strSampleSetLower.ToLower();
    isGammaJetsSample = strSampleSetLower.Contains("gjets");
  }
  if (isGammaJetsSample) MELAout << "Sample " << strSampleSet << " is a gamma+jets - like sample" << endl;
  else MELAout << "Sample " << strSampleSet << " is a DY-like sample" << endl;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure(period, "hadoop:200101");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);
  float lumi = SampleHelpers::getIntegratedLuminosity(period);

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(
    {
      OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kMuEle,
      OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle,
      OffshellTriggerHelpers::kSinglePho
    }
  );

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  VertexHandler vertexHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampleList);

  std::vector<TFile*> finput_DYJets;
  std::vector<TH1F*> h1D_rewgt_0j_DYJets;
  std::vector<TH1F*> h1D_rewgt_1j_DYJets;
  std::vector<TH1F*> h1D_rewgt_2j_DYJets;
  std::vector<TH2F*> h2D_rewgt_0j_DYJets;
  std::vector<TH2F*> h2D_rewgt_1j_DYJets;
  std::vector<TH2F*> h2D_rewgt_2j_DYJets;
  std::vector<TFile*> finput_GJets;
  std::vector<TH1F*> h1D_rewgt_0j_GJets;
  std::vector<TH1F*> h1D_rewgt_1j_GJets;
  std::vector<TH1F*> h1D_rewgt_2j_GJets;
  std::vector<TH2F*> h2D_rewgt_0j_GJets;
  std::vector<TH2F*> h2D_rewgt_1j_GJets;
  std::vector<TH2F*> h2D_rewgt_2j_GJets;
  if (istep>=1 && isGammaJetsSample){
    TString const strinputcore = Form("output/GammaJetsClosure/%s/%s/Step0", SampleHelpers::theDataPeriod.Data(), strdate.Data());
    TString cinput_rewgt_DYJets = strinputcore + "/DYJets_all.root";
    TString cinput_rewgt_GJets = strinputcore + "/GJets_all.root";

    finput_DYJets.push_back(TFile::Open(cinput_rewgt_DYJets, "read"));
    h1D_rewgt_0j_DYJets.push_back((TH1F*) finput_DYJets.back()->Get("Nvtxs_0j_rewgt")); HelperFunctions::wipeOverUnderFlows(h1D_rewgt_0j_DYJets.back(), false, true);
    h1D_rewgt_1j_DYJets.push_back((TH1F*) finput_DYJets.back()->Get("Nvtxs_1j_rewgt")); HelperFunctions::wipeOverUnderFlows(h1D_rewgt_1j_DYJets.back(), false, true);
    h1D_rewgt_2j_DYJets.push_back((TH1F*) finput_DYJets.back()->Get("Nvtxs_2j_rewgt")); HelperFunctions::wipeOverUnderFlows(h1D_rewgt_2j_DYJets.back(), false, true);

    finput_GJets.push_back(TFile::Open(cinput_rewgt_GJets, "read"));
    h1D_rewgt_0j_GJets.push_back((TH1F*) finput_GJets.back()->Get("Nvtxs_0j_rewgt")); HelperFunctions::wipeOverUnderFlows(h1D_rewgt_0j_GJets.back(), false, true);
    h1D_rewgt_1j_GJets.push_back((TH1F*) finput_GJets.back()->Get("Nvtxs_1j_rewgt")); HelperFunctions::wipeOverUnderFlows(h1D_rewgt_1j_GJets.back(), false, true);
    h1D_rewgt_2j_GJets.push_back((TH1F*) finput_GJets.back()->Get("Nvtxs_2j_rewgt")); HelperFunctions::wipeOverUnderFlows(h1D_rewgt_2j_GJets.back(), false, true);

    MELAout << "Extracted the weighting histograms from step 0!" << endl;
    curdir->cd();
  }
  if (istep>=2 && isGammaJetsSample){
    TString const strinputcore = Form("output/GammaJetsClosure/%s/%s/Step1", SampleHelpers::theDataPeriod.Data(), strdate.Data());
    TString cinput_rewgt_DYJets = strinputcore + "/DYJets_all.root";
    TString cinput_rewgt_GJets = strinputcore + "/GJets_all.root";

    finput_DYJets.push_back(TFile::Open(cinput_rewgt_DYJets, "read"));
    h2D_rewgt_0j_DYJets.push_back((TH2F*) finput_DYJets.back()->Get("pTll_etall_0j_rewgt")); HelperFunctions::wipeOverUnderFlows(h2D_rewgt_0j_DYJets.back(), false, true);
    h2D_rewgt_1j_DYJets.push_back((TH2F*) finput_DYJets.back()->Get("pTll_etall_1j_rewgt")); HelperFunctions::wipeOverUnderFlows(h2D_rewgt_1j_DYJets.back(), false, true);
    h2D_rewgt_2j_DYJets.push_back((TH2F*) finput_DYJets.back()->Get("pTll_etall_2j_rewgt")); HelperFunctions::wipeOverUnderFlows(h2D_rewgt_2j_DYJets.back(), false, true);

    finput_GJets.push_back(TFile::Open(cinput_rewgt_GJets, "read"));
    h2D_rewgt_0j_GJets.push_back((TH2F*) finput_GJets.back()->Get("pTll_etall_0j_rewgt")); HelperFunctions::wipeOverUnderFlows(h2D_rewgt_0j_GJets.back(), false, true);
    h2D_rewgt_1j_GJets.push_back((TH2F*) finput_GJets.back()->Get("pTll_etall_1j_rewgt")); HelperFunctions::wipeOverUnderFlows(h2D_rewgt_1j_GJets.back(), false, true);
    h2D_rewgt_2j_GJets.push_back((TH2F*) finput_GJets.back()->Get("pTll_etall_2j_rewgt")); HelperFunctions::wipeOverUnderFlows(h2D_rewgt_2j_GJets.back(), false, true);

    MELAout << "Extracted the weighting histograms from step 1!" << endl;
    curdir->cd();
  }

  TString const stroutputcore = Form("output/GammaJetsClosure/%s/%s/Step%i", SampleHelpers::theDataPeriod.Data(), strdate.Data(), istep);
  gSystem->Exec(Form("mkdir -p %s", stroutputcore.Data()));

  bool isFirstSample = true;
  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (isData && theGlobalSyst != SystematicsHelpers::sNominal) continue;

    TString cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    float xsec = 1;
    float sum_wgts = (isData ? 1 : 0);
    const int nEntries = sample_tree.getSelectedNEvents();

    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);
      
      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      sample_tree.silenceUnused();

      MELAout << "Initial MC loop over " << nEntries << " events to determine sample normalization:" << endl;
      for (int ev=0; ev<nEntries; ev++){
        HelperFunctions::progressbar(ev, nEntries);
        sample_tree.getSelectedEvent(ev);

        if (ev==0){
          sample_tree.getVal("xsec", xsec);
          sample_tree.releaseBranch("xsec");
        }

        genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
        auto const& genInfo = genInfoHandler.getGenInfo();
        float genwgt = genInfo->getGenWeight(true);

        simEventHandler.constructSimEvent(theGlobalSyst);
        float puwgt = simEventHandler.getPileUpWeight();
        float wgt = genwgt * puwgt;
        sum_wgts += wgt;
      }
      MELAout << "Sum of weights: " << sum_wgts << endl;
    }

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

    MELAout << "Completed getting the rest of the handles..." << endl;
    sample_tree.silenceUnused();

    // Create output tree
    TString cSample = SampleHelpers::getDatasetDirectoryCoreName(strSample.Data()).data();
    HelperFunctions::replaceString(cSample, "/", "_");
    TString stroutput = stroutputcore + "/" + cSample + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    ExtendedBinning bins_Nvtxs(50, 0., 100.);
    ExtendedBinning bins_pTll({ 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 250, 500, 1000, 3000 });
    ExtendedBinning bins_pTll_coarse({ 0, 5, 35, 50, 75, 100, 150, 200, 300, 500, 1000, 3000 });
    ExtendedBinning bins_etall({ -2.5, -2.0, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5 });
    ExtendedBinning bins_mZZ({ 0, 100, 180, 200, 225, 250, 300, 400, 500, 1000, 3000 });
    ExtendedBinning bins_Njets({ 0., 1., 2., 3. }, "N_{j}");
    ExtendedBinning bins_pt_jet_leadingpt({ 30., 40., 50., 60., 80., 100., 130., 165., 200., 250., 300., 400., 500., 1000., 3000. });
    ExtendedBinning bins_pt_jet_subleadingpt({ 30., 40., 50., 75., 100., 150., 200., 300., 500., 1000., 3000. });
    ExtendedBinning bins_min_abs_dPhi_j_puppimet(50, 0, TMath::Pi());
    ExtendedBinning bins_abs_dPhi_lljets_puppimet(50, 0, TMath::Pi());
    ExtendedBinning bins_puppimet_pTmiss(50, 0., 500.);
    ExtendedBinning bins_puppimet_pTmiss_over_pTlljets(40, 0., 8.);

    TH1F h1D_Nvtxs_0j_rewgt("Nvtxs_0j_rewgt", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()); h1D_Nvtxs_0j_rewgt.Sumw2(); h1D_Nvtxs_0j_rewgt.GetXaxis()->SetTitle("N_{vtx}");
    TH1F h1D_Nvtxs_1j_rewgt("Nvtxs_1j_rewgt", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()); h1D_Nvtxs_1j_rewgt.Sumw2(); h1D_Nvtxs_1j_rewgt.GetXaxis()->SetTitle("N_{vtx}");
    TH1F h1D_Nvtxs_2j_rewgt("Nvtxs_2j_rewgt", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()); h1D_Nvtxs_2j_rewgt.Sumw2(); h1D_Nvtxs_2j_rewgt.GetXaxis()->SetTitle("N_{vtx}");

    TH2F h2D_pTll_etall_0j_rewgt("pTll_etall_0j_rewgt", "", bins_pTll.getNbins(), bins_pTll.getBinning(), bins_etall.getNbins(), bins_etall.getBinning()); h2D_pTll_etall_0j_rewgt.Sumw2(); h2D_pTll_etall_0j_rewgt.GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h2D_pTll_etall_0j_rewgt.GetYaxis()->SetTitle("#eta_{ll}");
    TH2F h2D_pTll_etall_1j_rewgt("pTll_etall_1j_rewgt", "", bins_pTll.getNbins(), bins_pTll.getBinning(), bins_etall.getNbins(), bins_etall.getBinning()); h2D_pTll_etall_1j_rewgt.Sumw2(); h2D_pTll_etall_1j_rewgt.GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h2D_pTll_etall_1j_rewgt.GetYaxis()->SetTitle("#eta_{ll}");
    TH2F h2D_pTll_etall_2j_rewgt("pTll_etall_2j_rewgt", "", bins_pTll.getNbins(), bins_pTll.getBinning(), bins_etall.getNbins(), bins_etall.getBinning()); h2D_pTll_etall_2j_rewgt.Sumw2(); h2D_pTll_etall_2j_rewgt.GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h2D_pTll_etall_2j_rewgt.GetYaxis()->SetTitle("#eta_{ll}");

    std::unordered_map<TString, std::vector<TH1F>> var_histcoll_map;
    std::vector<std::pair<float, float>> met_thrs; met_thrs.reserve(4);
    met_thrs.emplace_back(-1.f, 20.f); // No MET
    met_thrs.emplace_back(20.f, 50.f); // Loose MET
    met_thrs.emplace_back(50.f, 85.f); // Medium MET
    met_thrs.emplace_back(85.f, -1.f); // Tight MET

    var_histcoll_map["pTll_0j"] = std::vector<TH1F>{
      TH1F("pTll_0j_noMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_0j_looseMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_0j_mediumMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_0j_tightMET", "", bins_pTll_coarse.getNbins(), bins_pTll_coarse.getBinning())
    };
    for (auto& hist:var_histcoll_map["pTll_0j"]) hist.GetXaxis()->SetTitle("p_{T}^{ll} (GeV)");
    var_histcoll_map["pTll_1j"] = std::vector<TH1F>{
      TH1F("pTll_1j_noMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_1j_looseMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_1j_mediumMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_1j_tightMET", "", bins_pTll_coarse.getNbins(), bins_pTll_coarse.getBinning())
    };
    for (auto& hist:var_histcoll_map["pTll_1j"]) hist.GetXaxis()->SetTitle("p_{T}^{ll} (GeV)");
    var_histcoll_map["pTll_2j"] = std::vector<TH1F>{
      TH1F("pTll_2j_noMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_2j_looseMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_2j_mediumMET", "", bins_pTll.getNbins(), bins_pTll.getBinning()),
      TH1F("pTll_2j_tightMET", "", bins_pTll_coarse.getNbins(), bins_pTll_coarse.getBinning())
    };
    for (auto& hist:var_histcoll_map["pTll_2j"]) hist.GetXaxis()->SetTitle("p_{T}^{ll} (GeV)");

    var_histcoll_map["mZZ_0j"] = std::vector<TH1F>{
      TH1F("mZZ_0j_noMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_0j_looseMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_0j_mediumMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_0j_tightMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning())
    };
    for (auto& hist:var_histcoll_map["mZZ_0j"]) hist.GetXaxis()->SetTitle("m_{ZZ} (GeV)");
    var_histcoll_map["mZZ_1j"] = std::vector<TH1F>{
      TH1F("mZZ_1j_noMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_1j_looseMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_1j_mediumMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_1j_tightMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning())
    };
    for (auto& hist:var_histcoll_map["mZZ_1j"]) hist.GetXaxis()->SetTitle("m_{ZZ} (GeV)");
    var_histcoll_map["mZZ_2j"] = std::vector<TH1F>{
      TH1F("mZZ_2j_noMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_2j_looseMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_2j_mediumMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning()),
      TH1F("mZZ_2j_tightMET", "", bins_mZZ.getNbins(), bins_mZZ.getBinning())
    };
    for (auto& hist:var_histcoll_map["mZZ_2j"]) hist.GetXaxis()->SetTitle("m_{ZZ} (GeV)");

    var_histcoll_map["Njets"] = std::vector<TH1F>{
      TH1F("Njets_noMET", "", bins_Njets.getNbins(), bins_Njets.getBinning()),
      TH1F("Njets_looseMET", "", bins_Njets.getNbins(), bins_Njets.getBinning()),
      TH1F("Njets_mediumMET", "", bins_Njets.getNbins(), bins_Njets.getBinning()),
      TH1F("Njets_tightMET", "", bins_Njets.getNbins(), bins_Njets.getBinning())
    };
    for (auto& hist:var_histcoll_map["Njets"]) hist.GetXaxis()->SetTitle("N_{j}");

    var_histcoll_map["Nvtxs_0j"] = std::vector<TH1F>{
      TH1F("Nvtxs_0j_noMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_0j_looseMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_0j_mediumMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_0j_tightMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning())
    };
    for (auto& hist:var_histcoll_map["Nvtxs_0j"]) hist.GetXaxis()->SetTitle("N_{vtx}");
    var_histcoll_map["Nvtxs_1j"] = std::vector<TH1F>{
      TH1F("Nvtxs_1j_noMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_1j_looseMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_1j_mediumMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_1j_tightMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning())
    };
    for (auto& hist:var_histcoll_map["Nvtxs_1j"]) hist.GetXaxis()->SetTitle("N_{vtx}");
    var_histcoll_map["Nvtxs_2j"] = std::vector<TH1F>{
      TH1F("Nvtxs_2j_noMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_2j_looseMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_2j_mediumMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning()),
      TH1F("Nvtxs_2j_tightMET", "", bins_Nvtxs.getNbins(), bins_Nvtxs.getBinning())
    };
    for (auto& hist:var_histcoll_map["Nvtxs_2j"]) hist.GetXaxis()->SetTitle("N_{vtx}");

    var_histcoll_map["pt_jet_leadingpt_1j"] = std::vector<TH1F>{
      TH1F("pt_jet_leadingpt_1j_noMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning()),
      TH1F("pt_jet_leadingpt_1j_looseMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning()),
      TH1F("pt_jet_leadingpt_1j_mediumMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning()),
      TH1F("pt_jet_leadingpt_1j_tightMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning())
    };
    for (auto& hist:var_histcoll_map["pt_jet_leadingpt_1j"]) hist.GetXaxis()->SetTitle("p_{T}^{j1} (GeV)");
    var_histcoll_map["pt_jet_leadingpt_2j"] = std::vector<TH1F>{
      TH1F("pt_jet_leadingpt_2j_noMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning()),
      TH1F("pt_jet_leadingpt_2j_looseMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning()),
      TH1F("pt_jet_leadingpt_2j_mediumMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning()),
      TH1F("pt_jet_leadingpt_2j_tightMET", "", bins_pt_jet_leadingpt.getNbins(), bins_pt_jet_leadingpt.getBinning())
    };
    for (auto& hist:var_histcoll_map["pt_jet_leadingpt_2j"]) hist.GetXaxis()->SetTitle("p_{T}^{j1} (GeV)");

    var_histcoll_map["pt_jet_subleadingpt_2j"] = std::vector<TH1F>{
      TH1F("pt_jet_subleadingpt_2j_noMET", "", bins_pt_jet_subleadingpt.getNbins(), bins_pt_jet_subleadingpt.getBinning()),
      TH1F("pt_jet_subleadingpt_2j_looseMET", "", bins_pt_jet_subleadingpt.getNbins(), bins_pt_jet_subleadingpt.getBinning()),
      TH1F("pt_jet_subleadingpt_2j_mediumMET", "", bins_pt_jet_subleadingpt.getNbins(), bins_pt_jet_subleadingpt.getBinning()),
      TH1F("pt_jet_subleadingpt_2j_tightMET", "", bins_pt_jet_subleadingpt.getNbins(), bins_pt_jet_subleadingpt.getBinning())
    };
    for (auto& hist:var_histcoll_map["pt_jet_subleadingpt_2j"]) hist.GetXaxis()->SetTitle("p_{T}^{j2} (GeV)");

    var_histcoll_map["min_abs_dPhi_j_puppimet_1j"] = std::vector<TH1F>{
      TH1F("min_abs_dPhi_j_puppimet_1j_noMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning()),
      TH1F("min_abs_dPhi_j_puppimet_1j_looseMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning()),
      TH1F("min_abs_dPhi_j_puppimet_1j_mediumMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning()),
      TH1F("min_abs_dPhi_j_puppimet_1j_tightMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning())
    };
    for (auto& hist:var_histcoll_map["min_abs_dPhi_j_puppimet_1j"]) hist.GetXaxis()->SetTitle("Min. |#phi_{j} - #phi_{miss}^{PUPPI}|");
    var_histcoll_map["min_abs_dPhi_j_puppimet_2j"] = std::vector<TH1F>{
      TH1F("min_abs_dPhi_j_puppimet_2j_noMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning()),
      TH1F("min_abs_dPhi_j_puppimet_2j_looseMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning()),
      TH1F("min_abs_dPhi_j_puppimet_2j_mediumMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning()),
      TH1F("min_abs_dPhi_j_puppimet_2j_tightMET", "", bins_min_abs_dPhi_j_puppimet.getNbins(), bins_min_abs_dPhi_j_puppimet.getBinning())
    };
    for (auto& hist:var_histcoll_map["min_abs_dPhi_j_puppimet_2j"]) hist.GetXaxis()->SetTitle("Min. |#phi_{j} - #phi_{miss}^{PUPPI}|");

    var_histcoll_map["abs_dPhi_lljets_puppimet_0j"] = std::vector<TH1F>{
      TH1F("abs_dPhi_lljets_puppimet_0j_noMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_0j_looseMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_0j_mediumMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_0j_tightMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning())
    };
    for (auto& hist:var_histcoll_map["abs_dPhi_lljets_puppimet_0j"]) hist.GetXaxis()->SetTitle("|#phi_{ll} - #phi_{miss}^{PUPPI}|");
    var_histcoll_map["abs_dPhi_lljets_puppimet_1j"] = std::vector<TH1F>{
      TH1F("abs_dPhi_lljets_puppimet_1j_noMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_1j_looseMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_1j_mediumMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_1j_tightMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning())
    };
    for (auto& hist:var_histcoll_map["abs_dPhi_lljets_puppimet_1j"]) hist.GetXaxis()->SetTitle("|#phi_{ll+j} - #phi_{miss}^{PUPPI}|");
    var_histcoll_map["abs_dPhi_lljets_puppimet_2j"] = std::vector<TH1F>{
      TH1F("abs_dPhi_lljets_puppimet_2j_noMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_2j_looseMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_2j_mediumMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning()),
      TH1F("abs_dPhi_lljets_puppimet_2j_tightMET", "", bins_abs_dPhi_lljets_puppimet.getNbins(), bins_abs_dPhi_lljets_puppimet.getBinning())
    };
    for (auto& hist:var_histcoll_map["abs_dPhi_lljets_puppimet_2j"]) hist.GetXaxis()->SetTitle("|#phi_{ll+jets} - #phi_{miss}^{PUPPI}|");

    TH1F puppimet_pTmiss_0j = TH1F("puppimet_pTmiss_0j", "", bins_puppimet_pTmiss.getNbins(), bins_puppimet_pTmiss.getBinning()); puppimet_pTmiss_0j.GetXaxis()->SetTitle("p_{T}^{miss,PUPPI} (GeV)");
    TH1F puppimet_pTmiss_1j = TH1F("puppimet_pTmiss_1j", "", bins_puppimet_pTmiss.getNbins(), bins_puppimet_pTmiss.getBinning()); puppimet_pTmiss_1j.GetXaxis()->SetTitle("p_{T}^{miss,PUPPI} (GeV)");
    TH1F puppimet_pTmiss_2j = TH1F("puppimet_pTmiss_2j", "", bins_puppimet_pTmiss.getNbins(), bins_puppimet_pTmiss.getBinning()); puppimet_pTmiss_2j.GetXaxis()->SetTitle("p_{T}^{miss,PUPPI} (GeV)");

    TH1F puppimet_pTmiss_over_pTlljets_0j = TH1F("puppimet_pTmiss_over_pTlljets_0j_noMET", "", bins_puppimet_pTmiss_over_pTlljets.getNbins(), bins_puppimet_pTmiss_over_pTlljets.getBinning()); puppimet_pTmiss_over_pTlljets_0j.GetXaxis()->SetTitle("p_{T}^{miss,PUPPI} / p_{T}^{ll}");
    TH1F puppimet_pTmiss_over_pTlljets_1j = TH1F("puppimet_pTmiss_over_pTlljets_1j_noMET", "", bins_puppimet_pTmiss_over_pTlljets.getNbins(), bins_puppimet_pTmiss_over_pTlljets.getBinning()); puppimet_pTmiss_over_pTlljets_1j.GetXaxis()->SetTitle("p_{T}^{miss,PUPPI} / p_{T}^{ll}");
    TH1F puppimet_pTmiss_over_pTlljets_2j = TH1F("puppimet_pTmiss_over_pTlljets_2j_noMET", "", bins_puppimet_pTmiss_over_pTlljets.getNbins(), bins_puppimet_pTmiss_over_pTlljets.getBinning()); puppimet_pTmiss_over_pTlljets_2j.GetXaxis()->SetTitle("p_{T}^{miss,PUPPI} / p_{T}^{ll}");

    // Setup stats tracking in histograms
    for (auto& it:var_histcoll_map){
      for (auto& hist:it.second){
        hist.Sumw2();
        hist.GetYaxis()->SetTitle("Number of events / bin");
      }
    }

    // Loop over the tree
    MELAout << "Starting to loop over " << nEntries << " events" << endl;
    unsigned int n_acc=0;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float trigwgt = eventFilter.getTriggerWeight(triggerCheckList);
      float genwgt = 1;
      float puwgt = 1;
      float wgt = 1;
      if (!isData){
        simEventHandler.constructSimEvent(theGlobalSyst);

        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();

        genwgt = genInfo->getGenWeight(true);
        puwgt = simEventHandler.getPileUpWeight();
      }

      wgt = genwgt * puwgt * trigwgt * xsec * 1000.f * (isData ? 1.f : lumi) / sum_wgts;
      if (wgt==0.f) continue;

      vertexHandler.constructVertices();
      unsigned int n_vertices_good = vertexHandler.getNGoodVertices();
      float n_vertices_good_f = n_vertices_good;

      muonHandler.constructMuons(theGlobalSyst);
      electronHandler.constructElectrons(theGlobalSyst);
      photonHandler.constructPhotons(theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      auto const& muons = muonHandler.getProducts();
      size_t n_muons_veto = 0;
      for (auto const& part:muons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++; }

      auto const& electrons = electronHandler.getProducts();
      size_t n_electrons_veto = 0;
      for (auto const& part:electrons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++; }

      auto const& photons = photonHandler.getProducts();
      size_t n_photons_tight = 0;
      PhotonObject* theChosenPhoton = nullptr;
      for (auto const& part:photons){
        if (ParticleSelectionHelpers::isTightParticle(part)){
          if (!theChosenPhoton) theChosenPhoton = part;
          n_photons_tight++;
        }
      }

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();

      if (!eventFilter.test2018HEMFilter((!isData ? &simEventHandler : nullptr), &electrons, &photons, &ak4jets, &ak8jets)) continue;

      auto const& puppimet = jetHandler.getPFPUPPIMET();
      float pt_puppimet = puppimet->pt();
      float phi_puppimet = puppimet->phi();

      ParticleObject::LorentzVector_t p4_ak4jets_tight;
      AK4JetObject* jet_leadingpt = nullptr;
      AK4JetObject* jet_subleadingpt = nullptr;
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);

          p4_ak4jets_tight = p4_ak4jets_tight + jet->p4();
          if (!jet_leadingpt) jet_leadingpt = jet;
          if (jet_leadingpt && !jet_subleadingpt && jet!=jet_leadingpt) jet_subleadingpt = jet;
        }
      }
      if (!ak4jets_tight_btagged.empty()) continue;
      unsigned int n_ak4jets_tight = ak4jets_tight.size();
      float n_ak4jets_tight_f = n_ak4jets_tight;
      float pt_jet_leadingpt = -1; if (jet_leadingpt) pt_jet_leadingpt = jet_leadingpt->pt();
      float pt_jet_subleadingpt = -1; if (jet_subleadingpt) pt_jet_subleadingpt = jet_subleadingpt->pt();

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();
      size_t n_dileptons_tight = 0;
      DileptonObject* theChosenDilepton = nullptr;
      for (auto const& dilepton:dileptons){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          n_dileptons_tight++;
        }
      }

      ParticleObject* theChosenDileptonProxy = nullptr;
      if (isGammaJetsSample && n_photons_tight==1) theChosenDileptonProxy = theChosenPhoton;
      if (!isGammaJetsSample && n_dileptons_tight==1) theChosenDileptonProxy = theChosenDilepton;
      //MELAout << "HERE: " << __LINE__ << endl;
      if (!theChosenDileptonProxy) continue;

      bool is_ee=false, is_mumu=false, is_emu=false;
      ParticleObject* leadingLepton = nullptr;
      ParticleObject* subleadingLepton = nullptr;
      float pTl1 = -1;
      float pTl2 = -1;
      if (!isGammaJetsSample){
        if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -121) is_ee=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -143) is_emu=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -169) is_mumu=true;
        leadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(0) : theChosenDilepton->daughter(1));
        subleadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(1) : theChosenDilepton->daughter(0));
        pTl1 = leadingLepton->pt();
        pTl2 = subleadingLepton->pt();
      }
      //MELAout << "HERE: " << __LINE__ << endl;
      if (!isGammaJetsSample && (pTl1<25.f || pTl2<20.f)) continue;

      float mll = theChosenDileptonProxy->m();
      float pTll = theChosenDileptonProxy->pt();
      float etall = theChosenDileptonProxy->eta();

      //MELAout << "HERE: " << __LINE__ << endl;
      if (pTll<35.f) continue;
      //MELAout << "HERE: " << __LINE__ << endl;
      if (
        !isGammaJetsSample && (
          (doZZselections && (is_emu || mll<81.f || mll>=101.f))
          ||
          (!doZZselections && mll<101.f)
          )
        ) continue;

      float pZmiss_approx = theChosenDileptonProxy->p4().Z();
      float etamiss_approx = theChosenDileptonProxy->eta();

      float min_abs_dPhi_j_puppimet = TMath::Pi();
      for (AK4JetObject* jet:ak4jets_tight){
        float dphi_tmp; HelperFunctions::deltaPhi(float(jet->phi()), phi_puppimet, dphi_tmp);
        min_abs_dPhi_j_puppimet = std::min(min_abs_dPhi_j_puppimet, std::abs(dphi_tmp));
      }
      bool pass_min_abs_dPhi_j_puppimet = (min_abs_dPhi_j_puppimet>=0.6);

      ParticleObject::LorentzVector_t p4_dileptonProxy;
      if (!isGammaJetsSample) p4_dileptonProxy = theChosenDileptonProxy->p4();
      else p4_dileptonProxy = ParticleObject::PolarLorentzVector_t(theChosenDileptonProxy->pt(), theChosenDileptonProxy->eta(), theChosenDileptonProxy->phi(), PDGHelpers::Zmass);

      float dPhi_lljets_puppimet; HelperFunctions::deltaPhi(float((p4_dileptonProxy + p4_ak4jets_tight).Phi()), phi_puppimet, dPhi_lljets_puppimet);
      float abs_dPhi_lljets_puppimet = std::abs(dPhi_lljets_puppimet);
      bool pass_abs_dPhi_lljets_puppimet = (abs_dPhi_lljets_puppimet>=2.6);

      float pTlljets = (p4_dileptonProxy + p4_ak4jets_tight).Pt();
      ParticleObject::LorentzVector_t puppimet_p4_ZZapprox; puppimet_p4_ZZapprox = ParticleObject::PolarLorentzVector_t(pt_puppimet, etamiss_approx, phi_puppimet, PDGHelpers::Zmass);
      float mZZ_puppimet = (puppimet_p4_ZZapprox + p4_dileptonProxy).M();

      if (istep>0){
        std::vector<TH1F*>* h1D_rewgt_DYJets=nullptr;
        std::vector<TH1F*>* h1D_rewgt_GJets=nullptr;
        std::vector<TH2F*>* h2D_rewgt_DYJets=nullptr;
        std::vector<TH2F*>* h2D_rewgt_GJets=nullptr;
        if (n_ak4jets_tight==0){
          h1D_rewgt_DYJets = &h1D_rewgt_0j_DYJets;
          h1D_rewgt_GJets = &h1D_rewgt_0j_GJets;
          h2D_rewgt_DYJets = &h2D_rewgt_0j_DYJets;
          h2D_rewgt_GJets = &h2D_rewgt_0j_GJets;
        }
        else if (n_ak4jets_tight==1){
          h1D_rewgt_DYJets = &h1D_rewgt_1j_DYJets;
          h1D_rewgt_GJets = &h1D_rewgt_1j_GJets;
          h2D_rewgt_DYJets = &h2D_rewgt_1j_DYJets;
          h2D_rewgt_GJets = &h2D_rewgt_1j_GJets;
        }
        else if (n_ak4jets_tight==2){
          h1D_rewgt_DYJets = &h1D_rewgt_2j_DYJets;
          h1D_rewgt_GJets = &h1D_rewgt_2j_GJets;
          h2D_rewgt_DYJets = &h2D_rewgt_2j_DYJets;
          h2D_rewgt_GJets = &h2D_rewgt_2j_GJets;
        }
        else continue;
        for (size_t iwgt=0; iwgt<h1D_rewgt_DYJets->size(); iwgt++){
          float extra_wgt=0;
          int nx = h1D_rewgt_DYJets->at(iwgt)->GetNbinsX(); int ix = std::min(nx, std::max(1, h1D_rewgt_DYJets->at(iwgt)->GetXaxis()->FindBin(pTll)));
          float extra_wgt_DYJets = h1D_rewgt_DYJets->at(iwgt)->GetBinContent(ix);
          float extra_wgt_GJets = h1D_rewgt_GJets->at(iwgt)->GetBinContent(ix);
          if (extra_wgt_GJets>0.f && extra_wgt_DYJets>0.f) extra_wgt = extra_wgt_DYJets / extra_wgt_GJets;
          wgt *= extra_wgt;
        }
        for (size_t iwgt=0; iwgt<h2D_rewgt_DYJets->size(); iwgt++){
          float extra_wgt=0;
          int nx = h2D_rewgt_DYJets->at(iwgt)->GetNbinsX(); int ix = std::min(nx, std::max(1, h2D_rewgt_DYJets->at(iwgt)->GetXaxis()->FindBin(pTll)));
          int ny = h2D_rewgt_DYJets->at(iwgt)->GetNbinsY(); int iy = std::min(ny, std::max(1, h2D_rewgt_DYJets->at(iwgt)->GetYaxis()->FindBin(etall)));
          float extra_wgt_DYJets = h2D_rewgt_DYJets->at(iwgt)->GetBinContent(ix, iy);
          float extra_wgt_GJets = h2D_rewgt_GJets->at(iwgt)->GetBinContent(ix, iy);
          if (extra_wgt_GJets>0.f && extra_wgt_DYJets>0.f) extra_wgt = extra_wgt_DYJets / extra_wgt_GJets;
          wgt *= extra_wgt;
        }
      }

      float puppimet_over_pTlljets = pt_puppimet / pTlljets;
      unsigned int i_metbin = 0;
      for (auto const& met_thr_pair:met_thrs){
        bool pass_puppimet_low = met_thr_pair.first<0.f || (puppimet_over_pTlljets>=std::pow(met_thr_pair.first / pTlljets, 1.5) && pt_puppimet>=met_thr_pair.first);
        bool pass_puppimet_high = met_thr_pair.second<0.f || (puppimet_over_pTlljets<std::pow(met_thr_pair.second / pTlljets, 1.5) || pt_puppimet<met_thr_pair.second);
        if (pass_puppimet_low && pass_puppimet_high){
          // Fill 1D histograms
          for (auto& var_histcoll_pair:var_histcoll_map){
            /*
            if (n_acc%100000==0) MELAout << "Trying to fill " << var_histcoll_pair.first
              << " with Nj=" << n_ak4jets_tight
              << ", pt_puppimet=" << pt_puppimet
              << ", puppimet_over_pTlljets=" << puppimet_over_pTlljets 
              << " in MET bin " << i_metbin
              << endl;
            */

            if (var_histcoll_pair.first.Contains("0j") && n_ak4jets_tight!=0) continue;
            else if (var_histcoll_pair.first.Contains("1j") && n_ak4jets_tight!=1) continue;
            else if (var_histcoll_pair.first.Contains("2j") && n_ak4jets_tight!=2) continue;
            else if (n_ak4jets_tight>2) continue;

            float* var = nullptr;
            if (var_histcoll_pair.first.Contains("Njets")) var = &n_ak4jets_tight_f;
            else if (var_histcoll_pair.first.Contains("Nvtxs")) var = &n_vertices_good_f;
            else if (var_histcoll_pair.first.Contains("mZZ")) var = &mZZ_puppimet;
            else if (var_histcoll_pair.first.Contains("pTll")) var = &pTll;
            else if (var_histcoll_pair.first.Contains("min_abs_dPhi_j_puppimet")) var = &min_abs_dPhi_j_puppimet;
            else if (var_histcoll_pair.first.Contains("abs_dPhi_lljets_puppimet")) var = &abs_dPhi_lljets_puppimet;
            else if (var_histcoll_pair.first.Contains("pt_jet_leadingpt")){ if (jet_leadingpt) var = &pt_jet_leadingpt; else continue; }
            else if (var_histcoll_pair.first.Contains("pt_jet_subleadingpt")){ if (jet_subleadingpt) var = &pt_jet_subleadingpt; else continue; }

            // Do not apply cuts on the variables plotted
            if (
              !(pass_min_abs_dPhi_j_puppimet || var==&min_abs_dPhi_j_puppimet)
              ||
              !(pass_abs_dPhi_lljets_puppimet || var==&abs_dPhi_lljets_puppimet)
              ) continue;

            //if (n_acc%100000==0) MELAout << "Filling " << var_histcoll_pair.first << " with Nj=" << n_ak4jets_tight << " (var = " << *var << ", weight = " << wgt << ")" << endl;

            if (var){
              var_histcoll_pair.second.at(i_metbin).Fill(*var, wgt);
              //if (n_acc%100000==0) MELAout << "Fill successful!" << endl;
            }
            else{
              MELAerr << "Histogram group " << var_histcoll_pair.first << " does not have a variable associated at Nj=" << n_ak4jets_tight << "!" << endl;
              assert(0);
            }
          }

          if (!(istep<1 || (pass_min_abs_dPhi_j_puppimet && pass_abs_dPhi_lljets_puppimet))) continue;

          // Fill MET-inclusive histograms
          if (n_ak4jets_tight == 0){
            puppimet_pTmiss_0j.Fill(pt_puppimet, wgt);
            puppimet_pTmiss_over_pTlljets_0j.Fill(puppimet_over_pTlljets, wgt);
          }
          else if (n_ak4jets_tight == 1){
            puppimet_pTmiss_1j.Fill(pt_puppimet, wgt);
            puppimet_pTmiss_over_pTlljets_1j.Fill(puppimet_over_pTlljets, wgt);
          }
          else if (n_ak4jets_tight == 2){
            puppimet_pTmiss_2j.Fill(pt_puppimet, wgt);
            puppimet_pTmiss_over_pTlljets_2j.Fill(puppimet_over_pTlljets, wgt);
          }

          if (!(pass_min_abs_dPhi_j_puppimet && pass_abs_dPhi_lljets_puppimet)) continue;

          // Fill 2D reweighting histogram only for bins except the last
          if (i_metbin != met_thrs.size()-1){
            if (n_ak4jets_tight == 0){
              h1D_Nvtxs_0j_rewgt.Fill(n_vertices_good_f, wgt);
              h2D_pTll_etall_0j_rewgt.Fill(pTll, etall, wgt);
            }
            else if (n_ak4jets_tight == 1){
              h1D_Nvtxs_1j_rewgt.Fill(n_vertices_good_f, wgt);
              h2D_pTll_etall_1j_rewgt.Fill(pTll, etall, wgt);
            }
            else if (n_ak4jets_tight == 2){
              h1D_Nvtxs_2j_rewgt.Fill(n_vertices_good_f, wgt);
              h2D_pTll_etall_2j_rewgt.Fill(pTll, etall, wgt);
            }
          }
        }

        i_metbin++;
      }

      n_acc++;
    }

    foutput->WriteTObject(&h1D_Nvtxs_0j_rewgt);
    foutput->WriteTObject(&h1D_Nvtxs_1j_rewgt);
    foutput->WriteTObject(&h1D_Nvtxs_2j_rewgt);
    foutput->WriteTObject(&h2D_pTll_etall_0j_rewgt);
    foutput->WriteTObject(&h2D_pTll_etall_1j_rewgt);
    foutput->WriteTObject(&h2D_pTll_etall_2j_rewgt);
    for (auto& it:var_histcoll_map){
      for (auto& hist:it.second){
        MELAout << "Histogram " << hist.GetName() << " has integral " << hist.Integral() << endl;
        foutput->WriteTObject(&hist);
      }
    }
    foutput->WriteTObject(&puppimet_pTmiss_0j);
    foutput->WriteTObject(&puppimet_pTmiss_1j);
    foutput->WriteTObject(&puppimet_pTmiss_2j);
    foutput->WriteTObject(&puppimet_pTmiss_over_pTlljets_0j);
    foutput->WriteTObject(&puppimet_pTmiss_over_pTlljets_1j);
    foutput->WriteTObject(&puppimet_pTmiss_over_pTlljets_2j);

    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  }

  for (auto* finput:finput_GJets) finput->Close();
  for (auto* finput:finput_DYJets) finput->Close();
}
