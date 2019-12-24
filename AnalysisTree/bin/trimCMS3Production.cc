#include <iostream>
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "HelperFunctions.h"


using namespace std;


void trimCMS3Production(TString infile, TString outfile, TString strfilter="(@muons_pt.size() + @electrons_pt.size())>1 || @photons_pt>0"){
  if (strfilter == "") return;

  cout << "Cloning and trimming " << infile << " to " << outfile << " with selection \"" << strfilter << "\"." << endl;

  TDirectory* curdir = gDirectory;

  TFile* finput = TFile::Open(infile, "read");
  TTree* intree = (TTree*) finput->Get("cms3ntuple/Events");
  intree->SetAutoSave(0);
  int nEntries = intree->GetEntries();
  cout << "\t- Found " << nEntries << " events in total." << endl;

  TFile* foutput = TFile::Open(outfile, "recreate");
  foutput->cd();
  TDirectory* subdir = foutput->mkdir("cms3ntuple");
  subdir->cd();

  cout << "\t- Beginning to trim and copy..." << endl;

  TTree* outtree = intree->CopyTree(strfilter);
  if (outtree) subdir->WriteTObject(outtree);
  subdir->Close();
  foutput->Close();
  finput->Close();

  curdir->cd();
}

int main(int argc, char** argv){
  if (argc!=4) return 1;
  cout << "Running " << argv[0] << "(\"" << argv[1] << "\", \"" << argv[2] << "\", \"" << argv[3] << "\")" << endl;
  trimCMS3Production(argv[1], argv[2], argv[3]);
  return 0;
}
