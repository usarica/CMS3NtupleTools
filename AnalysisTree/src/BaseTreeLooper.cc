#include <cassert>

#include <algorithm>
#include <utility>
#include <iterator>
#include <fstream>

#include "BaseTreeLooper.h"
#include "SampleHelpersCore.h"
#include "SamplesCore.h"

#include "SimEventHandler.h"
#include "GenInfoHandler.h"
#include "EventFilterHandler.h"
#include "RunLumiEventBlock.h"

#include "HelperFunctions.h"
#include "HostHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


BaseTreeLooper::BaseTreeLooper() :
  IvyBase(),

  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),

  maxNEvents(-1),
  eventIndex_begin(-1),
  eventIndex_end(-1),
  useChunkIndices(false),

  isData_currentTree(false),
  isQCD_currentTree(false),
  isGJets_HT_currentTree(false)
{
  set_pTG_exception_range(-1, -1);
  setExternalProductList();
  setCurrentOutputTree();
}
BaseTreeLooper::BaseTreeLooper(BaseTree* inTree, double wgt) :
  IvyBase(),

  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),

  maxNEvents(-1),
  eventIndex_begin(-1),
  eventIndex_end(-1),
  useChunkIndices(false),

  isData_currentTree(false),
  isQCD_currentTree(false),
  isGJets_HT_currentTree(false)
{
  this->addTree(inTree, wgt);
  set_pTG_exception_range(-1, -1);
  setExternalProductList();
  setCurrentOutputTree();
}
BaseTreeLooper::BaseTreeLooper(std::vector<BaseTree*> const& inTreeList) :
  IvyBase(),

  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),

  maxNEvents(-1),
  eventIndex_begin(-1),
  eventIndex_end(-1),
  useChunkIndices(false),

  isData_currentTree(false),
  isQCD_currentTree(false),
  isGJets_HT_currentTree(false),

  treeList(inTreeList)
{
  set_pTG_exception_range(-1, -1);
  setExternalProductList();
  setCurrentOutputTree();
}
BaseTreeLooper::~BaseTreeLooper(){}

void BaseTreeLooper::addTree(BaseTree* tree, double wgt){
  if (tree && !HelperFunctions::checkListVariable(this->treeList, tree)){
    this->treeList.push_back(tree);
    setExternalWeight(tree, wgt);
  }
}
void BaseTreeLooper::addTree(BaseTree* tree, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& syst_wgt_map){
  if (tree && !HelperFunctions::checkListVariable(this->treeList, tree)){
    this->treeList.push_back(tree);
    setExternalWeights(tree, syst_wgt_map);
  }
}

void BaseTreeLooper::addExternalFunction(TString fcnname, BaseTreeLooper::LooperExtFunction_t fcn){
  if (!fcn) return;
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) MELAerr << "BaseTreeLooper::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}
void BaseTreeLooper::addObjectHandler(IvyBase* handler){
  if (handler && !HelperFunctions::checkListVariable(this->registeredHandlers, handler)) this->registeredHandlers.push_back(handler);
}
void BaseTreeLooper::addSFHandler(ScaleFactorHandlerBase* handler){
  if (handler && !HelperFunctions::checkListVariable(this->registeredSFHandlers, handler)) this->registeredSFHandlers.push_back(handler);
}

void BaseTreeLooper::addReweightingBuilder(TString name, BulkReweightingBuilder* rewgtBuilder){
  if (!rewgtBuilder){
    MELAerr << "BaseTreeLooper::addReweightingBuilder: Reweighting builder " << name << " is null." << endl;
    return;
  }
  if (registeredRewgtBuilders.find(name)!=registeredRewgtBuilders.end()) MELAerr << "BaseTreeLooper::addReweightingBuilder: Reweighting builder " << name << " already exists but will override it regardless." << endl;
  registeredRewgtBuilders[name] = rewgtBuilder;
}

void BaseTreeLooper::addHLTMenu(TString name, std::vector< std::string > const& hltmenu){
  if (registeredHLTMenus.find(name)!=registeredHLTMenus.end()) MELAerr << "BaseTreeLooper::addHLTMenu: Simple HLT menu " << name << " already exists but will override it regardless." << endl;
  registeredHLTMenus[name] = hltmenu;
}
void BaseTreeLooper::addHLTMenu(TString name, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > const& hltmenu){
  if (registeredHLTMenuProperties.find(name)!=registeredHLTMenuProperties.end()) MELAerr << "BaseTreeLooper::addHLTMenu: HLT menu properties " << name << " already exists but will override it regardless." << endl;
  registeredHLTMenuProperties[name] = hltmenu;
}

void BaseTreeLooper::setMatrixElementList(std::vector<std::string> const& MElist, bool const& isGen){
  MELAout << "BaseTreeLooper::setMatrixElementList: Setting " << (isGen ? "gen." : "reco.") << " matrix elements:" << endl;
  for (auto const& sme:MElist) MELAout << '\t' << sme << endl;
  if (isGen) lheMElist = MElist;
  else recoMElist = MElist;
}
void BaseTreeLooper::setMatrixElementListFromFile(std::string fname, std::string const& MElistTypes, bool const& isGen){
  if (MElistTypes==""){
    if (this->verbosity>=TVar::ERROR) MELAerr << "BaseTreeLooper::setMatrixElementListFromFile: The ME list types must be specified." << endl;
    assert(0);
  }
  HostHelpers::ExpandEnvironmentVariables(fname);
  if (!HostHelpers::FileReadable(fname.data())){
    if (this->verbosity>=TVar::ERROR) MELAerr << "BaseTreeLooper::setMatrixElementListFromFile: File " << fname << " is not readable." << endl;
    assert(0);
  }

  std::vector<std::string> MEtypes;
  HelperFunctions::splitOptionRecursive(MElistTypes, MEtypes, ',', true);

  // Read the file and collect the MEs
  std::vector<std::string> MElist;
  ifstream fin;
  fin.open(fname.c_str());
  if (fin.good()){
    bool acceptString = false;
    while (!fin.eof()){
      std::string str_in="";
      getline(fin, str_in);
      HelperFunctions::lstrip(str_in);
      HelperFunctions::lstrip(str_in, "\"\'");
      HelperFunctions::rstrip(str_in); HelperFunctions::rstrip(str_in, ",\"\'");

      if (str_in=="" || str_in.find('#')==0) continue;
      else if (str_in.find(']')!=std::string::npos){
        acceptString = false;
        continue;
      }

      bool isMEline = (str_in.find("Name")!=std::string::npos);
      for (auto const& MEtype:MEtypes){
        if (str_in.find(MEtype)!=std::string::npos){
          isMEline = false;
          acceptString = true;
        }
      }

      if (isMEline && acceptString){
        if (isGen && str_in.find("isGen:1")==std::string::npos){
          if (this->verbosity>=TVar::ERROR) MELAerr << "BaseTreeLooper::setMatrixElementListFromFile: ME string " << str_in << " is not a gen. ME while the acquisition is done for gen. MEs!" << endl;
          continue;
        }
        MElist.push_back(str_in);
      }
    }
  }
  fin.close();

  if (!MElist.empty()) setMatrixElementList(MElist, isGen);
  else{
    if (this->verbosity>=TVar::ERROR) MELAerr << "BaseTreeLooper::setMatrixElementListFromFile: File " << fname << " does not contain any of the ME types " << MEtypes << "." << endl;
  }
}

void BaseTreeLooper::setExternalWeight(BaseTree* tree, double const& wgt){
  if (!tree) return;
  if (this->verbosity>=TVar::INFO && !HelperFunctions::checkListVariable(treeList, tree)) MELAout
    << "BaseTreeLooper::setExternalWeight: Warning! Tree " << tree->sampleIdentifier
    << " is not in the list of input trees, but a weight of " << wgt << " is being assigned to it."
    << endl;
  auto it = globalWeights.find(tree);
  if (it==globalWeights.end()) globalWeights[tree] = std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double>();
  globalWeights[tree][SystematicsHelpers::nSystematicVariations] = wgt;
}
void BaseTreeLooper::setExternalWeights(BaseTree* tree, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& wgts){
  if (!tree) return;
  if (this->verbosity>=TVar::INFO && !HelperFunctions::checkListVariable(treeList, tree)) MELAout
    << "BaseTreeLooper::setExternalWeight: Warning! Tree " << tree->sampleIdentifier
    << " is not in the list of input trees, but a set of weights is being assigned to it."
    << endl;
  globalWeights[tree] = wgts;
}

void BaseTreeLooper::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef = extProductListRef;
  else this->productListRef = &(this->productList);
}

void BaseTreeLooper::setCurrentOutputTree(BaseTree* extTree){
  // Set current product tree to this output tree
  this->currentProductTree = extTree;
  // Print warning if this tree is not in the output tree collection, but do not add.
  if (extTree && !HelperFunctions::checkListVariable(this->productTreeList, extTree)){
    if (this->verbosity>=TVar::INFO) MELAout << "BaseTreeLooper::setCurrentOutputTree: Current output tree is not in the output tree collection." << endl;
  }
  // Make sure product list collects some events before flushing
  this->productListRef = &(this->productList);
}

void BaseTreeLooper::addOutputTree(BaseTree* extTree){
  if (extTree){
    if (!HelperFunctions::checkListVariable(this->productTreeList, extTree)) this->productTreeList.push_back(extTree);
    this->setCurrentOutputTree(extTree);
  }
}
void BaseTreeLooper::addOutputTrees(std::vector<BaseTree*> trees){
  for (auto const& tt:trees) addOutputTree(tt);
}

void BaseTreeLooper::setMaximumEvents(int n){ maxNEvents = n; }
void BaseTreeLooper::setEventIndexRange(int istart, int iend){ eventIndex_begin = istart; eventIndex_end = iend; }


void BaseTreeLooper::addProduct(SimpleEntry& product, unsigned int* ev_rec){
  this->productListRef->push_back(product);
  if (ev_rec) (*ev_rec)++;

  // Record products to external tree
  this->recordProductsToTree();
}

void BaseTreeLooper::recordProductsToTree(){
  if (!this->currentProductTree) return;

  auto it_tree = firstTreeOutput.find(this->currentProductTree);
  if (it_tree == firstTreeOutput.cend()){
    firstTreeOutput[this->currentProductTree] = true;
    it_tree = firstTreeOutput.find(this->currentProductTree);
  }

  BaseTree::writeSimpleEntries(this->productListRef->cbegin(), this->productListRef->cend(), this->currentProductTree, it_tree->second);

  it_tree->second = false;
  this->clearProducts();
}

bool BaseTreeLooper::wrapTree(BaseTree* tree){
  if (!tree) return false;
  bool res = true;

  // Sample flags
  TString const& sid = tree->sampleIdentifier;
  this->isData_currentTree = SampleHelpers::checkSampleIsData(sid);
  this->isQCD_currentTree = !this->isData_currentTree && sid.Contains("QCD") && sid.Contains("HT");
  this->isGJets_HT_currentTree = !this->isData_currentTree && sid.Contains("GJets_HT");
  if (!this->isData_currentTree && (sid.Contains("ZGTo2NuG") || sid.Contains("ZGTo2LG")) && sid.Contains("amcatnloFXFX") && !sid.Contains("PtG-130")) set_pTG_exception_range(-1, 130);
  else set_pTG_exception_range(-1, -1);

  if (this->isData_currentTree){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL); this->addConsumed<TYPE>(#NAME); this->defineConsumedSloppy(#NAME);
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
  }

  for (auto const& handler:registeredHandlers){
    bool isHandlerForSim = (dynamic_cast<GenInfoHandler*>(handler) != nullptr || dynamic_cast<SimEventHandler*>(handler) != nullptr);
    EventFilterHandler* eventFilter = dynamic_cast<EventFilterHandler*>(handler);
    if (eventFilter){
      bool isFirstInputFile = (tree == treeList.front());
      eventFilter->setTrackDataEvents(this->isData_currentTree);
      eventFilter->setCheckUniqueDataEvent(this->isData_currentTree && !isFirstInputFile);
    }
    if (!this->isData_currentTree || (this->isData_currentTree && !isHandlerForSim)) res &= handler->wrapTree(tree);
  }
  res &= IvyBase::wrapTree(tree);
  return res;
}

void BaseTreeLooper::incrementSelection(TString const& strsel, unsigned int inc){
  bool isFound = false;
  for (auto& pp:selection_string_count_pairs){
    if (pp.first == strsel){
      pp.second += inc;
      isFound = true;
    }
  }
  if (!isFound) selection_string_count_pairs.emplace_back(strsel, inc);
}

void BaseTreeLooper::loop(bool keepProducts){
  if (!looperFunction){
    MELAerr << "BaseTreeLooper::loop: The looper function is not registered. Please register it using BadeTreeLooper::setLooperFunction." << endl;
    return;
  }

  // Count total number of events
  int nevents_total = 0;
  bool hasDataTrees = false;
  bool hasSimTrees = false;
  for (auto& tree:treeList){
    nevents_total += tree->getNEvents();
    if (SampleHelpers::checkSampleIsData(tree->sampleIdentifier)) hasDataTrees = true;
    else hasSimTrees = true;
  }
  if (hasDataTrees && hasSimTrees){
    MELAerr << "BaseTreeLooper::loop: Looping over both real data and simulation at the same time is forbidden. Please run to separate loopers." << endl;
    assert(0);
  }

  // Adjust event ranges to actual event indices
  int const eventIndex_begin_orig = eventIndex_begin;
  int const eventIndex_end_orig = eventIndex_end;
  if (this->useChunkIndices && eventIndex_end>0){
    const int ichunk = eventIndex_begin;
    const int nchunks = eventIndex_end;
    if (hasSimTrees){
      // Assign the range over total number of events
      int ev_inc = static_cast<int>(float(nevents_total)/float(nchunks));
      int ev_rem = nevents_total - ev_inc*nchunks;
      eventIndex_begin = ev_inc*ichunk + std::min(ev_rem, ichunk);
      eventIndex_end = ev_inc*(ichunk+1) + std::min(ev_rem, ichunk+1);
      MELAout << "BaseTreeLooper::loop: A simulation loop will proceed. The requested event range is [" << eventIndex_begin << ", " << eventIndex_end << ")." << endl;
    }
    else if (hasDataTrees){
      // Assign the range over run numbers
      double const lumi_total = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
      auto const& runnumber_lumi_pairs = SampleHelpers::getRunNumberLumiPairsForDataPeriod(SampleHelpers::getDataPeriod());
      int const nruns_total = runnumber_lumi_pairs.size();
      int const nruns_inc = static_cast<int>(static_cast<double>(nruns_total) / static_cast<double>(nchunks));
      int const nruns_rem = nruns_total - nruns_inc*nchunks;

      int const idx_firstRun = nruns_inc*ichunk + std::min(nruns_rem, ichunk);
      int const idx_firstRun_next = nruns_inc*(ichunk+1) + std::min(nruns_rem, ichunk+1);
      eventIndex_begin = (idx_firstRun>=nruns_total ? 0 : (int) runnumber_lumi_pairs.at(idx_firstRun).first);
      eventIndex_end = (idx_firstRun_next-1>=nruns_total || idx_firstRun_next<=idx_firstRun ? 0 : (int) runnumber_lumi_pairs.at(idx_firstRun_next-1).first);

      double lumi_acc = 0;
      for (auto const& runnumber_lumi_pair:runnumber_lumi_pairs){
        if (
          (eventIndex_begin<0 || eventIndex_begin<=(int) runnumber_lumi_pair.first)
          &&
          (eventIndex_end<0 || eventIndex_end>=(int) runnumber_lumi_pair.first)
          ) lumi_acc += runnumber_lumi_pair.second;
      }

      MELAout << "BaseTreeLooper::loop: A real data loop will proceed. The requested run range is [" << eventIndex_begin << ", " << eventIndex_end << "]. Total luminosity covered will be " << lumi_acc << " / " << lumi_total << "." << endl;
    }
  }

  // Build the MEs if they are specified
  if (!lheMElist.empty() || !recoMElist.empty()){
    // Set up MELA (done only once inside CMS3MELAHelpers)
    CMS3MELAHelpers::setupMela(SampleHelpers::getDataYear(), 125., TVar::ERROR);
    // If there are output trees, set the output trees of the MEblock.
    // Do this before building the branches.
    MELAout << "Setting up ME block output trees..." << endl;
    for (auto const& outtree_:productTreeList){
      MELAout << "\t- Extracting valid output trees" << endl;
      if (!outtree_) MELAerr << "Output tree is NULL!" << endl;
      std::vector<TTree*> const& treelist_ = outtree_->getValidTrees();
      MELAout << "\t- Registering " << treelist_.size() << " trees" << endl;
      for (auto const& tree_:treelist_) MEblock.addRefTree(tree_);
      MELAout << "\t- Done" << endl;
    }
    // Build the MEs
    MELAout << "Building the MEs..." << endl;
    if (!lheMElist.empty()) this->MEblock.buildMELABranches(lheMElist, true);
    if (!recoMElist.empty()) this->MEblock.buildMELABranches(recoMElist, false);
  }

  // Loop over the trees
  unsigned int ev_traversed=0;
  unsigned int ev_acc=0;
  unsigned int ev_rec=0;
  for (auto& tree:treeList){
    // Skip the tree if it cannot be wrapped
    if (!(this->wrapTree(tree))) continue;

#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr;
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    float MHval = -1;
    SampleIdStorageType sampleIdOpt = kNoStorage;
    if (this->isData_currentTree){
      sampleIdOpt = kStoreByRunAndEventNumber;
      bool rlenPresent = true;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) rlenPresent &= this->getConsumed(#NAME, NAME);
      RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
      if (!rlenPresent){
        if (this->verbosity>=TVar::ERROR) MELAerr << "BaseTreeLooper::loop: Run number, lumi block, or event number are not consumed properly..." << endl;
        assert(0);
      }
    }
    else{
      MHval = SampleHelpers::findPoleMass(tree->sampleIdentifier);
      if (MHval>0.f) sampleIdOpt = kStoreByMH;
    }

    auto it_globalWgt = globalWeights.find(tree);
    if (it_globalWgt==globalWeights.cend()){
      if (this->verbosity>=TVar::ERROR) MELAerr << "BaseTreeLooper::loop: " << tree->sampleIdentifier << " does not have any weights assigned..." << endl;
      assert(0);
    }

    const int nevents = tree->getNEvents();
    MELAout << "BaseTreeLooper::loop: Looping over " << nevents << " events in " << tree->sampleIdentifier << "..." << endl;
    for (int ev=0; ev<nevents; ev++){
      if (
        SampleHelpers::doSignalInterrupt==1
        ||
        (maxNEvents>=0 && (int) ev_rec==maxNEvents)
        ) break;

      bool doAccumulate = true;
      if (this->isData_currentTree){
        if (eventIndex_begin>0 || eventIndex_end>0) doAccumulate = (
          tree->updateBranch(ev, "RunNumber", false)
          &&
          (eventIndex_begin<0 || static_cast<int>(*RunNumber)>=eventIndex_begin)
          &&
          (eventIndex_end<0 || static_cast<int>(*RunNumber)<=eventIndex_end)
          );
      }
      else doAccumulate = (
        (eventIndex_begin<0 || (int) ev_traversed>=eventIndex_begin)
        &&
        (eventIndex_end<0 || (int) ev_traversed<eventIndex_end)
        );

      if (doAccumulate){
        if (tree->getEvent(ev)){
          SimpleEntry product;
          if (sampleIdOpt==kStoreByRunAndEventNumber){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) product.setNamedVal<TYPE>(#NAME, *NAME);
            RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
          }
          else if (sampleIdOpt==kStoreByMH) product.setNamedVal("SampleMHVal", MHval);
          if (tree->isValidEvent()){
            if (this->looperFunction(this, it_globalWgt->second, product)){
              if (keepProducts) this->addProduct(product, &ev_rec);
            }
          }
        }
        ev_acc++;
      }

      HelperFunctions::progressbar(ev, nevents);
      ev_traversed++;
    }

    if (!selection_string_count_pairs.empty()){
      MELAout << "BaseTreeLooper::loop: Number of events passing each selection type:" << endl;
      for (auto& pp:selection_string_count_pairs){
        MELAout << (pp.first.BeginsWith("\t") ? "\t" : "\t- ") << pp.first << ": " << pp.second << endl;
      }
    }
    resetSelectionCounts();
  } // End loop over the trees
  MELAout << "BaseTreeLooper::loop: Total number of products: " << ev_rec << " / " << ev_acc << " / " << ev_traversed << endl;

  // Restore original event index values
  eventIndex_begin = eventIndex_begin_orig;
  eventIndex_end = eventIndex_end_orig;
}

std::vector<SimpleEntry> const& BaseTreeLooper::getProducts() const{ return *productListRef; }

void BaseTreeLooper::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "BaseTreeLooper::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "BaseTreeLooper::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void BaseTreeLooper::clearProducts(){ productListRef->clear(); }
