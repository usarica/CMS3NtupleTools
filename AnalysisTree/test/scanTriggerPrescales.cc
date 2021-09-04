#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


struct TriggerPropInfo{
  unsigned int run;
  unsigned int lumi;
  float prescale;

  TriggerPropInfo();
  TriggerPropInfo(TriggerPropInfo const& other);
  TriggerPropInfo(unsigned int const& run_, unsigned int const& lumi_, float const& prescale_);

  static bool isEarlier(TriggerPropInfo const& tpi_first, TriggerPropInfo const& tpi_second);

};

TriggerPropInfo::TriggerPropInfo() :
  run(0),
  lumi(0),
  prescale(0)
{}
TriggerPropInfo::TriggerPropInfo(TriggerPropInfo const& other) :
  run(other.run),
  lumi(other.lumi),
  prescale(other.prescale)
{}
TriggerPropInfo::TriggerPropInfo(unsigned int const& run_, unsigned int const& lumi_, float const& prescale_) :
  run(run_),
  lumi(lumi_),
  prescale(prescale_)
{}

bool TriggerPropInfo::isEarlier(TriggerPropInfo const& tpi_first, TriggerPropInfo const& tpi_second){
  return !(tpi_first.run>tpi_second.run || (tpi_first.run==tpi_second.run && tpi_first.lumi>=tpi_second.lumi));
}


void scanTriggerPrescales(TString period, TString prodVersion){
  SampleHelpers::configure(period, "store:"+prodVersion);
  if (!SampleHelpers::testDataPeriodIsLikeData()) return;

  std::vector<std::string> triggerCheckList;
  for (unsigned int itype=0; itype<(unsigned int) TriggerHelpers::nTriggerTypes; itype++){
    auto tmpvec = TriggerHelpers::getHLTMenus((TriggerHelpers::TriggerType) itype);
    for (auto const& ss:tmpvec){
      if (!HelperFunctions::checkListVariable(triggerCheckList, ss)) triggerCheckList.push_back(ss);
    }
  }

  auto const& valid_run_lumi_pairs = SampleHelpers::getRunNumberLumiPairsForDataPeriod(SampleHelpers::getDataPeriod());

  EventFilterHandler eventFilter;

  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(false);
  eventFilter.setCheckHLTPathRunRanges(false);

  TString strSampleSet = Form("Run%s", SampleHelpers::getDataPeriod().Data());
  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  std::unordered_map<unsigned int, std::vector<unsigned int>> run_lumilist_map;
  std::vector<std::vector<TriggerPropInfo>> triggerpropinfolists(triggerCheckList.size());

  IVYout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);

    TString const cinputcore = SampleHelpers::getDatasetDirectoryName(strSample, true);
    TString const cinput = SampleHelpers::getDatasetFileName(strSample, true);
    if (cinput=="") continue;
    IVYout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    const int nEntries = sample_tree.getSelectedNEvents();
    int ev_start = 0;
    int ev_end = nEntries;
    IVYout << "Looping over " << nEntries << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    unsigned int* RunNumber;
    unsigned int* LuminosityBlock;
    sample_tree.getValRef("RunNumber", RunNumber);
    sample_tree.getValRef("LuminosityBlock", LuminosityBlock);

    IVYout << "Completed getting the rest of the handles..." << endl;
    sample_tree.silenceUnused();

    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);
      if (ev%10000==0) IVYout << sample_tree.sampleIdentifier << " events: " << ev << " / " << nEntries << endl;

      auto it_run = run_lumilist_map.find(*RunNumber);
      if (it_run == run_lumilist_map.end()){
        run_lumilist_map[*RunNumber] = std::vector<unsigned int>();
        it_run = run_lumilist_map.find(*RunNumber);
      }
      if (!HelperFunctions::checkListVariable(it_run->second, *LuminosityBlock)) it_run->second.push_back(*LuminosityBlock);
      else continue;

      eventFilter.constructFilters(nullptr);
      auto const& hltpaths = eventFilter.getHLTPaths();

      {
        auto it_triggerpropinfolists = triggerpropinfolists.begin();
        for (auto const& strig:triggerCheckList){
          float trigwgt=0;
          for (auto const& hltpath:hltpaths){
            if (hltpath->name.find(strig)!=std::string::npos){
              trigwgt = hltpath->HLTprescale * hltpath->L1prescale; // Combines L1 and HLT prescales
              break;
            }
          }
          it_triggerpropinfolists->emplace_back(*RunNumber, *LuminosityBlock, trigwgt);
          it_triggerpropinfolists++;
        }
      }

    }

  }

  TString coutput_main = Form("output/TriggerPrescales_%s", prodVersion.Data());
  gSystem->Exec(Form("mkdir -p %s", coutput_main.Data()));
  {
    auto it_triggerpropinfolists = triggerpropinfolists.begin();
    for (auto const& strig:triggerCheckList){
      std::sort(it_triggerpropinfolists->begin(), it_triggerpropinfolists->end(), TriggerPropInfo::isEarlier);

      IVYout << "*****" << endl;
      IVYout << strig << endl;
      IVYout << "*****" << endl;

      TString stroutput = Form("%s/%s_%s.txt", coutput_main.Data(), strig.data(), SampleHelpers::getDataPeriod().Data());
      IVYout.open(stroutput.Data());

      IVYout << "Run,lumi,prescale" << endl;
      for (auto const& triggerpropinfo:(*it_triggerpropinfolists)){
        if (triggerpropinfo.prescale == 1.f) continue;
        IVYout << triggerpropinfo.run << ',' << triggerpropinfo.lumi << ',' << triggerpropinfo.prescale << endl;
      }
      for (auto const& run_lumi_pair:valid_run_lumi_pairs){
        auto const& run = run_lumi_pair.first;
        bool found = false;
        for (auto const& triggerpropinfo:(*it_triggerpropinfolists)){
          if (triggerpropinfo.run == run){
            found = true;
            break;
          }
        }
        if (!found) IVYout << run << ',' << -1 << ',' << 0 << endl;
      }

      IVYout.close();
      IVYout << "*****" << endl;

      it_triggerpropinfolists++;

      SampleHelpers::addToCondorTransferList(stroutput);
    }
  }

}
