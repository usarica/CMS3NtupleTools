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


void getHistograms(TString strSampleSet="EGamma", TString period="2018"){
  gStyle->SetOptStat(0);

  SampleHelpers::configure(period, "200101");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);
  float lumi = SampleHelpers::getIntegratedLuminosity(period);

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(OffshellTriggerHelpers::kSinglePho);

  // Get handlers
  PhotonHandler photonHandler;
  EventFilterHandler eventFilter;

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  TString const stroutputcore = Form("output/PhotonTriggerScan/%s", SampleHelpers::theDataPeriod.Data());
  gSystem->Exec(Form("mkdir -p %s", stroutputcore.Data()));

  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (!isData) continue;

    TString cinput = SampleHelpers::getDatasetFileName(strSample);
    IVYout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    float xsec=1;
    const int nEntries = sample_tree.getSelectedNEvents();

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    IVYout << "Completed getting the rest of the handles..." << endl;
    sample_tree.silenceUnused();

    // Create output
    TString cSample = SampleHelpers::getDatasetDirectoryCoreName(strSample.Data()).data();
    HelperFunctions::replaceString(cSample, "/", "_");
    TString stroutput = stroutputcore + "/" + cSample + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    ExtendedBinning ptbins(200, 0, 500);
    ExtendedBinning etabins(25, -2.5, 2.5);
    std::unordered_map<std::string, TH1F> trigger_photonPtDistribution_map;
    std::unordered_map<std::string, TH2F> trigger_photonDistribution_map;
    std::unordered_map<std::string, unsigned int> trigger_nFilledPhotons_map;
    for (auto const& strtrig:triggerCheckList){
      TString hnamecore = strtrig.data();
      HelperFunctions::replaceString(hnamecore, "*", "");
      trigger_photonPtDistribution_map[strtrig] = TH1F(Form("h1D_Pt_%s", hnamecore.Data()), "", ptbins.getNbins(), ptbins.getBinning());
      trigger_photonPtDistribution_map[strtrig].Sumw2();
      trigger_photonDistribution_map[strtrig] = TH2F(Form("h2D_%s", hnamecore.Data()), "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning());
      trigger_photonDistribution_map[strtrig].Sumw2();
      trigger_nFilledPhotons_map[strtrig] = 0;
    }
    // Loop over the tree
    IVYout << "Starting to loop over " << nEntries << " events" << endl;
    unsigned int n_acc=0;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);
      if (ev%10000==0) IVYout << "Event " << n_acc << " / " << ev << "..." << endl;
      //if (ev==10000) break;

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float trigwgt = eventFilter.getTriggerWeight(triggerCheckList);
      if (trigwgt==0.) continue;

      photonHandler.constructPhotons(SystematicsHelpers::sNominal);
      auto const& photons = photonHandler.getProducts();

      if (!eventFilter.test2018HEMFilter(nullptr, nullptr, &photons, nullptr, nullptr)) continue;

      for (auto const& hltpath:eventFilter.getHLTPaths()){
        if (std::find(triggerCheckList.cbegin(), triggerCheckList.cend(), hltpath->name) == triggerCheckList.cend()) continue;
        if (hltpath->L1prescale>1) IVYout << hltpath->name << " has L1 and HLT prescales of " << hltpath->L1prescale << ", " << hltpath->HLTprescale << endl;
      }

      for (auto const& part:photons){
        if (!ParticleSelectionHelpers::isTightParticle(part)) continue;
        for (auto const& strtrig:triggerCheckList){
          trigwgt = eventFilter.getTriggerWeight({ strtrig });
          if (trigwgt==0.) continue;
          trigger_photonPtDistribution_map[strtrig].Fill(part->pt(), trigwgt);
          trigger_photonDistribution_map[strtrig].Fill(part->pt(), part->eta(), trigwgt);
          trigger_nFilledPhotons_map[strtrig]++;
        }
      }

      bool needMorePhotons=false;
      for (auto const& strtrig:triggerCheckList) needMorePhotons |= (trigger_nFilledPhotons_map[strtrig]<10000);
      if (!needMorePhotons) break;

      n_acc++;
      if (n_acc==100000) break;
    }

    for (auto it:trigger_photonPtDistribution_map){
      IVYout << it.first << " 1D integral: " << it.second.Integral() << endl;
      foutput->WriteTObject(&(it.second));
    }
    for (auto it:trigger_photonDistribution_map) foutput->WriteTObject(&(it.second));
    foutput->Close();
  }
}
