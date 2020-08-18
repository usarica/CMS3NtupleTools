#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>

#include "BaseTreeLooper.h"
#include "SamplesCore.h"

#include "SimEventHandler.h"
#include "GenInfoHandler.h"
#include "EventFilterHandler.h"

#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


BaseTreeLooper::BaseTreeLooper() :
  IvyBase(),
  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),
  maxNEvents(-1),
  eventIndex_begin(-1),
  eventIndex_end(-1)
{
  setExternalProductList();
  setCurrentOutputTree();
}
BaseTreeLooper::BaseTreeLooper(BaseTree* inTree) :
  IvyBase(),
  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),
  maxNEvents(-1),
  eventIndex_begin(-1),
  eventIndex_end(-1)
{
  this->addTree(inTree);
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

  treeList(inTreeList)
{
  setExternalProductList();
  setCurrentOutputTree();
}
BaseTreeLooper::~BaseTreeLooper(){}

void BaseTreeLooper::addTree(BaseTree* tree){
  if (tree && !HelperFunctions::checkListVariable(this->treeList, tree)) this->treeList.push_back(tree);
}

void BaseTreeLooper::addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, BaseTree*, SimpleEntry&)){
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

void BaseTreeLooper::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef = extProductListRef;
  else this->productListRef = &(this->productList);
}

void BaseTreeLooper::setCurrentOutputTree(BaseTree* extTree){
  this->currentProductTree = extTree;
  this->productListRef = &(this->productList); // To make sure product list collects some events before flushing
}

void BaseTreeLooper::addOutputTree(BaseTree* extTree){
  if (extTree && !HelperFunctions::checkListVariable(productTreeList, extTree)){
    productTreeList.push_back(extTree);
    setCurrentOutputTree(extTree);
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
  this->clearProducts();
}

bool BaseTreeLooper::wrapTree(BaseTree* tree){
  if (!tree) return false;
  bool res = true;
  bool const isData = SampleHelpers::checkSampleIsData(tree->sampleIdentifier);
  for (auto const& handler:registeredHandlers){
    bool isHandlerForSim = (dynamic_cast<GenInfoHandler*>(handler) != nullptr || dynamic_cast<SimEventHandler*>(handler) != nullptr);
    EventFilterHandler* eventFilter = dynamic_cast<EventFilterHandler*>(handler);
    if (eventFilter){
      bool isFirstInputFile = (tree == treeList.front());
      eventFilter->setTrackDataEvents(isData);
      eventFilter->setCheckUniqueDataEvent(isData && !isFirstInputFile);
    }
    if (!isData || (isData && !isHandlerForSim)) res &= handler->wrapTree(tree);
  }
  res &= IvyBase::wrapTree(tree);
  return res;
}


void BaseTreeLooper::loop(bool keepProducts){
  if (!looperFunction){
    MELAerr << "BaseTreeLooper::loop: The looper function is not registered. Please register it using BadeTreeLooper::setLooperFunction." << endl;
    return;
  }

  // Loop over the trees
  unsigned int ev_traversed=0;
  unsigned int ev_acc=0;
  unsigned int ev_rec=0;
  for (auto& tree:treeList){
    // Skip the tree if it cannot be wrapped
    if (!(this->wrapTree(tree))) continue;

    float globalTreeWeight = 1;
    auto it_globalWgt = globalWeights.find(tree);
    if (it_globalWgt!=globalWeights.cend()) globalTreeWeight = it_globalWgt->second;

    const int nevents = tree->getNEvents();
    MELAout << "BaseTreeLooper::loop: Looping over " << nevents << " events in " << tree->sampleIdentifier << "..." << endl;
    for (int ev=0; ev<nevents; ev++){
      if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;
      if (
        (eventIndex_begin<0 || (int) ev_traversed>=eventIndex_begin)
        &&
        (eventIndex_end<0 || (int) ev_traversed<eventIndex_end)
        ){
        if (tree->getEvent(ev)){
          SimpleEntry product;
          if (tree->isValidEvent()){
            if (this->looperFunction(this, tree, globalTreeWeight, product)){
              if (keepProducts) this->addProduct(product, &ev_rec);
            }
          }
        }
        ev_acc++;
      }
      HelperFunctions::progressbar(ev, nevents);
      ev_traversed++;
    }
  } // End loop over the trees
  MELAout << "BaseTreeLooper::loop: Total number of products: " << ev_rec << " / " << ev_acc << " / " << ev_traversed << endl;
}

std::vector<SimpleEntry> const& BaseTreeLooper::getProducts() const{ return *productListRef; }

void BaseTreeLooper::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "BaseTreeLooper::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "BaseTreeLooper::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void BaseTreeLooper::clearProducts(){ productListRef->clear(); }
