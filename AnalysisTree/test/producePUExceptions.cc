#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


void producePUExceptions(TString strSampleSet, TString period, TString prodVersion, TString strdate=""){
  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop:"+prodVersion);
  MELAout << "Input directory: " << SampleHelpers::theInputDirectory << endl;
  MELAout << "Input tag: " << SampleHelpers::theSamplesTag << endl;

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  TString const stroutputcore = Form("output/PUExceptions/%s", strdate.Data());

  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (isData) return;

    TString const cinputcore = SampleHelpers::getDatasetDirectoryName(strSample);
    TString const cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);
    sample_tree.bookBranch<float>("n_true_int", 0.f);
    float* n_true_int=nullptr;
    sample_tree.getValRef("n_true_int", n_true_int);
    sample_tree.silenceUnused();

    // Create output
    TString stroutput = stroutputcore;
    gSystem->Exec(Form("mkdir -p %s", stroutput.Data()));
    stroutput += "/" + sample_tree.sampleIdentifier + ".root";
    MELAout << "Creating output file " << stroutput << "..." << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TH1F hpu("pileup", sample_tree.sampleIdentifier, 100, 0, 100);

    const int nEntries = sample_tree.getSelectedNEvents();
    int ev_start = 0;
    int ev_end = nEntries;
    MELAout << "Looping over " << nEntries << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      if (HelperFunctions::checkVarNanInf(*n_true_int) && (*n_true_int)>=0) hpu.Fill(*n_true_int);
    }

    foutput->WriteTObject(&hpu);
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  }
}
