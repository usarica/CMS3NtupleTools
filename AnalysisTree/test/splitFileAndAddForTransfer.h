// No includes, should be added at the end of common_includes.h
void splitFileAndAddForTransfer(TString const& stroutput){
  using namespace std;
  using namespace MELAStreamHelpers;

  // Trivial case: If not running on condor, there is no need to transfer. Just exit.
  if (!SampleHelpers::checkRunOnCondor()){
    SampleHelpers::addToCondorTransferList(stroutput);
    return;
  }

  TDirectory* curdir = gDirectory;
  size_t const size_limit = std::pow(1024, 3);

  TFile* finput = TFile::Open(stroutput, "read");
  curdir->cd();

  size_t const size_input = finput->GetSize();
  size_t const nchunks = size_input/size_limit+1;
  std::vector<TString> fnames; fnames.reserve(nchunks);
  if (nchunks>1){
    std::vector<TFile*> foutputlist; foutputlist.reserve(nchunks);

    MELAout << "splitFileAndAddForTransfer: Splitting " << stroutput << " into " << nchunks << " chunks:" << endl;
    for (size_t ichunk=0; ichunk<nchunks; ichunk++){
      TString fname = stroutput;
      TString strchunk = Form("_chunk_%zu_of_%zu%s", ichunk, nchunks, ".root");
      HelperFunctions::replaceString<TString, TString const>(fname, ".root", strchunk);
      MELAout << "\t- Making new file " << fname << "..." << endl;
      TFile* foutput = TFile::Open(fname, "recreate");
      foutputlist.push_back(foutput);
      fnames.push_back(fname);
    }

    std::vector<TDirectory*> outputdirs; outputdirs.reserve(nchunks);
    for (auto& ff:foutputlist) outputdirs.push_back(ff);
    HelperFunctions::distributeObjects(finput, outputdirs, TVar::INFO);

    for (auto& ff:foutputlist) ff->Close();
    MELAout << "\t- Splitting is completed." << endl;
  }
  else{
    MELAout << "splitFileAndAddForTransfer: " << stroutput << " will not be split into chunks." << endl;
    fnames.push_back(stroutput);
  }

  finput->Close();
  curdir->cd();

  for (auto const& fname:fnames) SampleHelpers::addToCondorTransferList(fname);
}
