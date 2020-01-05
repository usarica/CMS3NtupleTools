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


void get2DParallelAndPerpendicularComponents(TVector3 axis, TVector3 ref, float& parallel, float& perp){
  TVector3 unitAxis = TVector3(axis.X(), axis.Y(), 0).Unit();
  TVector3 refPerp = TVector3(ref.X(), ref.Y(), 0);
  parallel = unitAxis.Dot(refPerp);
  perp = unitAxis.Cross(refPerp).Z();
}

void addDataTreeList(std::vector<TString>& list){
  std::vector<TString> tmplist;
  tmplist.push_back("EGamma");
  tmplist.push_back("MuonEG");
  tmplist.push_back("SingleMuon");
  tmplist.push_back("DoubleMuon");
  for (auto const& s:tmplist) SampleHelpers::constructSamplesList(s, SystematicsHelpers::sNominal, list);
}
void addDataTreeList_Combined(std::vector<TString>& list){
  switch (SampleHelpers::theDataYear){
  case 2016:
    list.push_back("Run2016B-17Jul2018_ver2");
    list.push_back("Run2016C-17Jul2018");
    list.push_back("Run2016D-17Jul2018");
    list.push_back("Run2016E-17Jul2018");
    list.push_back("Run2016F-17Jul2018");
    list.push_back("Run2016G-17Jul2018");
    list.push_back("Run2016H-17Jul2018");
    break;
  case 2017:
    list.push_back("Run2017B-31Mar2018-v1");
    list.push_back("Run2017C-31Mar2018-v1");
    list.push_back("Run2017D-31Mar2018-v1");
    list.push_back("Run2017E-31Mar2018-v1");
    list.push_back("Run2017F-31Mar2018-v1");
    list.push_back("Run2017F-09May2018-v1");
    break;
  case 2018:
    list.push_back("Run2018A-17Sep2018");
    list.push_back("Run2018B-17Sep2018");
    list.push_back("Run2018C-17Sep2018");
    list.push_back("Run2018D-PromptReco");
    break;
  }
}
void addMCTreeList(std::vector<TString>& list){
  std::vector<TString> tmplist;
  tmplist.push_back("GJets_topology");
  for (auto const& s:tmplist) SampleHelpers::constructSamplesList(s, SystematicsHelpers::sNominal, list);
}


void createGammaTrees(TString strSampleSet, TString period, SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  SampleHelpers::configure(period, "191212");

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
  VertexHandler vertexHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  EventFilterHandler eventFilter;

  std::vector<TString> sampleList;
  addDataTreeList(sampleList);
  addMCTreeList(sampleList);
  {
    std::vector<TString> tmplist;
    for (TString const& strSample:sampleList){
      bool const isData = SampleHelpers::checkSampleIsData(strSample);
      if (
        (strSampleSet == "Data" && isData)
        ||
        (strSampleSet == "MC" && !isData)
        ||
        (strSample.Contains(strSampleSet))
        ) tmplist.push_back(strSample);
    }
    std::swap(sampleList, tmplist);
  }

  TString const stroutputcore = Form("output/GammaTrees/%s", SampleHelpers::theDataPeriod.Data());
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

    float xsec=1;
    float sum_wgts=(isData ? 1 : 0);
    const int nEntries = sample_tree.getSelectedNEvents();

    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);
      
      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      sample_tree.silenceUnused();

      MELAout << "Initial MC loop to determine sample normalization:" << endl;
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
    TString stroutput = stroutputcore + "/" + cSample + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    BaseTree* outtree = new BaseTree("AnalysisTree");
    outtree->setAutoSave(0);
    bool firstEvent = true;
    SimpleEntry commonEntry;

    // Loop over the tree
    MELAout << "Starting to loop over " << nEntries << " events" << endl;
    unsigned int n_acc=0;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);
      if (ev%10000==0) MELAout << "Event " << n_acc << " / " << ev << "..." << endl;
      //if (ev==10000) break;

      float genmet_met=0, genmet_metPhi=0;

      float genwgt = 1;
      float puwgt = 1;
      float trigwgt = 1;
      float wgt = 1;
      if (!isData){
        simEventHandler.constructSimEvent(theGlobalSyst);

        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();

        genwgt = genInfo->getGenWeight(true);
        genmet_met = genInfo->extras.genmet_met;
        genmet_metPhi = genInfo->extras.genmet_metPhi;
        puwgt = simEventHandler.getPileUpWeight();
      }

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      trigwgt = eventFilter.getTriggerWeight(triggerCheckList);

      wgt = genwgt * puwgt * trigwgt * xsec * (isData ? 1.f : lumi) / sum_wgts;

      if (wgt==0.f) continue;

      vertexHandler.constructVertices();
      unsigned int n_vertices_good = vertexHandler.getNGoodVertices();

      muonHandler.constructMuons(theGlobalSyst);
      auto const& muons = muonHandler.getProducts();
      size_t n_muons_veto = 0;
      for (auto const& part:muons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++; }
      if (n_muons_veto>0) continue;

      electronHandler.constructElectrons(theGlobalSyst);
      auto const& electrons = electronHandler.getProducts();
      size_t n_electrons_veto = 0;
      for (auto const& part:electrons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++; }
      if (n_electrons_veto>0) continue;

      photonHandler.constructPhotons(theGlobalSyst);
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

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();

      if (!eventFilter.test2018HEMFilter((!isData ? &simEventHandler : nullptr), &electrons, &photons, &ak4jets, &ak8jets)) continue;

      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& pfmet = jetHandler.getPFMET();

      float jetHT = 0;
      ParticleObject::LorentzVector_t ak4jets_sumP;
      AK4JetObject* jet_leadingpt = nullptr;
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);

          ak4jets_sumP = ak4jets_sumP + jet->p4();
          jetHT += jet->pt();
          if (!jet_leadingpt) jet_leadingpt = jet;
        }
      }
      if (!ak4jets_tight_btagged.empty()) continue;

      unsigned int n_ak4jets_tight = ak4jets_tight.size();
      if (n_ak4jets_tight==0) continue;

      float photon_pt = theChosenPhoton->pt();
      float photon_eta = theChosenPhoton->eta();
      float photon_phi = theChosenPhoton->phi();
      float photon_mass = theChosenPhoton->m();

      float ak4jet_leadingpt_pt = jet_leadingpt->pt();
      float ak4jet_leadingpt_eta = jet_leadingpt->eta();
      float ak4jet_leadingpt_phi = jet_leadingpt->phi();
      float ak4jet_leadingpt_mass = jet_leadingpt->m();

      float ak4jets_sumP_pt = ak4jets_sumP.pt();
      float ak4jets_sumP_eta = ak4jets_sumP.eta();
      float ak4jets_sumP_phi = ak4jets_sumP.phi();
      float ak4jets_sumP_mass = ak4jets_sumP.mass();

      float puppimet_met = puppimet->pt();
      float puppimet_metPhi = puppimet->phi();
      float pfmet_met = pfmet->pt();
      float pfmet_metPhi = pfmet->phi();

      // TVectors
      TLorentzVector vec_puppimet; vec_puppimet.SetPtEtaPhiM(puppimet_met, 0, puppimet_metPhi, 0);
      TLorentzVector vec_pfmet; vec_pfmet.SetPtEtaPhiM(pfmet_met, 0, pfmet_metPhi, 0);
      TLorentzVector vec_ak4jets_sumP; vec_ak4jets_sumP.SetPtEtaPhiM(ak4jets_sumP_pt, ak4jets_sumP_eta, ak4jets_sumP_phi, ak4jets_sumP_mass);
      TVector3 photonAxis; photonAxis.SetPtEtaPhi(theChosenPhoton->pt(), theChosenPhoton->eta(), theChosenPhoton->phi());

      // Parallel and perpendicular components
      float uParallel, uPerp;
      float puppimet_parallel, puppimet_perp;
      float pfmet_parallel, pfmet_perp;
      get2DParallelAndPerpendicularComponents(photonAxis, vec_ak4jets_sumP.Vect(), uParallel, uPerp);
      get2DParallelAndPerpendicularComponents(photonAxis, vec_puppimet.Vect(), puppimet_parallel, puppimet_perp);
      get2DParallelAndPerpendicularComponents(photonAxis, vec_pfmet.Vect(), pfmet_parallel, pfmet_perp);

      // Record
      commonEntry.setNamedVal("n_vertices_good", n_vertices_good);

      commonEntry.setNamedVal("genwgt", genwgt / sum_wgts);
      commonEntry.setNamedVal("puwgt", puwgt);
      commonEntry.setNamedVal("trigwgt", trigwgt);
      commonEntry.setNamedVal("xsec", xsec);
      commonEntry.setNamedVal("wgt", wgt);

      commonEntry.setNamedVal("n_ak4jets_tight", n_ak4jets_tight);

      commonEntry.setNamedVal("jetHT", jetHT);

      commonEntry.setNamedVal("photon_pt", photon_pt);
      commonEntry.setNamedVal("photon_eta", photon_eta);
      commonEntry.setNamedVal("photon_phi", photon_phi);
      commonEntry.setNamedVal("photon_mass", photon_mass);

      commonEntry.setNamedVal("ak4jet_leadingpt_pt", ak4jet_leadingpt_pt);
      commonEntry.setNamedVal("ak4jet_leadingpt_eta", ak4jet_leadingpt_eta);
      commonEntry.setNamedVal("ak4jet_leadingpt_phi", ak4jet_leadingpt_phi);
      commonEntry.setNamedVal("ak4jet_leadingpt_mass", ak4jet_leadingpt_mass);

      commonEntry.setNamedVal("ak4jets_sumP_pt", ak4jets_sumP_pt);
      commonEntry.setNamedVal("ak4jets_sumP_eta", ak4jets_sumP_eta);
      commonEntry.setNamedVal("ak4jets_sumP_phi", ak4jets_sumP_phi);
      commonEntry.setNamedVal("ak4jets_sumP_mass", ak4jets_sumP_mass);

      commonEntry.setNamedVal("genmet_met", genmet_met);
      commonEntry.setNamedVal("genmet_metPhi", genmet_metPhi);

      commonEntry.setNamedVal("puppimet_met", puppimet_met);
      commonEntry.setNamedVal("puppimet_metPhi", puppimet_metPhi);
      commonEntry.setNamedVal("puppimet_parallel", puppimet_parallel);
      commonEntry.setNamedVal("puppimet_perp", puppimet_perp);

      commonEntry.setNamedVal("pfmet_met", pfmet_met);
      commonEntry.setNamedVal("pfmet_metPhi", pfmet_metPhi);
      commonEntry.setNamedVal("pfmet_parallel", pfmet_parallel);
      commonEntry.setNamedVal("pfmet_perp", pfmet_perp);

      commonEntry.setNamedVal("uParallel", uParallel);
      commonEntry.setNamedVal("uPerp", uPerp);

      if (firstEvent){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) outtree->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) outtree->putBranch(itb->first, &(itb->second));
        SIMPLE_DATA_OUTPUT_DIRECTIVES;
        VECTOR_DATA_OUTPUT_DIRECTIVES;
        DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES;
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
      }

      // Record whatever is in commonEntry into the tree.
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) outtree->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) outtree->setVal(itb->first, &(itb->second));
      SIMPLE_DATA_OUTPUT_DIRECTIVES;
      VECTOR_DATA_OUTPUT_DIRECTIVES;
      DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES;
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

      outtree->fill();
      n_acc++;
      firstEvent = false;
    }

    outtree->writeToFile(foutput);
    delete outtree;
    foutput->Close();
  }
}
