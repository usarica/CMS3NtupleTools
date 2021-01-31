#include <iostream>
#include <sstream>
#include "OffshellSampleHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


void printCMS3SampleGroup(TString strSampleSet, TString period, TString prodVersion, TString outfile){
  // Back up the stdout stream buffer and redirect it
  streambuf* stdoutbuf_bkp = cout.rdbuf();
  stringstream ss;
  cout.rdbuf(ss.rdbuf());

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  std::vector<TString> sgroups;
  HelperFunctions::splitOptionRecursive(strSampleSet, sgroups, ',', true);

  std::vector<TString> samples;
  for (auto const& sgroup:sgroups){
    for (unsigned int isyst=0; isyst<(unsigned int) SystematicsHelpers::nSystematicVariations; isyst++){
      SystematicsHelpers::SystematicVariationTypes syst = (SystematicsHelpers::SystematicVariationTypes) isyst;
      std::vector<TString> snames;
      SampleHelpers::constructSamplesList(sgroup, syst, snames);
      for (auto const& sname:snames){
        if (!HelperFunctions::checkListVariable(samples, sname)) samples.push_back(sname);
      }
    }
  }

  // Restore stdout stream buffer either now or after writing the output file
  bool const useStdout = (outfile=="stdout");
  if (!useStdout) MELAout.open(outfile.Data());
  else cout.rdbuf(stdoutbuf_bkp);
  for (auto const& sname:samples) MELAout << sname << endl;
  if (!useStdout){
    MELAout.close();
    cout.rdbuf(stdoutbuf_bkp);
  }
}

int main(int argc, char** argv){
  if (argc!=5){
    cout << "printCMS3SampleGroup usage:\n\tprintCMS3SampleGroup [sample_set] [period] [production_tag] [output_file]" << endl;
    cout << "\t- [sample_set] could be a comma-separated list of different sample groups." << endl;
    cout << "\t- [output_file] could be a file name or 'stdout'." << endl;
    return 1;
  }

  printCMS3SampleGroup(argv[1], argv[2], argv[3], argv[4]);

  return 0;
}
