#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>
#include "BaseTreeLooper.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


BaseTreeLooper::BaseTreeLooper() :
  IvyBase(),
  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),
  maxNEvents(-1),
  firstTreeOutput(true)
{
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::BaseTreeLooper(BaseTree* inTree) :
  IvyBase(),
  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),
  maxNEvents(-1),
  firstTreeOutput(true)
{
  this->addTree(inTree);
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::BaseTreeLooper(std::vector<BaseTree*> const& inTreeList) :
  IvyBase(),

  looperFunction(nullptr),
  registeredSyst(SystematicsHelpers::nSystematicVariations),
  maxNEvents(-1),
  firstTreeOutput(true),

  treeList(inTreeList)
{
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::~BaseTreeLooper(){}

void BaseTreeLooper::addTree(BaseTree* tree){
  if (tree && !HelperFunctions::checkListVariable(this->treeList, tree)) this->treeList.push_back(tree);
}

void BaseTreeLooper::addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper const*, BaseTree*, SimpleEntry&)){
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
  if (extProductListRef) this->productListRef=extProductListRef;
  else this->productListRef=&(this->productList);
}

void BaseTreeLooper::setExternalProductTree(BaseTree* extTree){
  this->productTree = extTree;
  this->productListRef = &(this->productList); // To make sure product list collects some events before flushing
  if (extTree) firstTreeOutput = true;
}

void BaseTreeLooper::setMaximumEvents(int n){ maxNEvents=n; }

void BaseTreeLooper::addProduct(SimpleEntry& product, unsigned int* ev_rec){
  this->productListRef->push_back(product);
  if (ev_rec) (*ev_rec)++;
}

void BaseTreeLooper::recordProductsToTree(){
  if (!this->productTree) return;
  BaseTree::writeSimpleEntries(this->productListRef->begin(), this->productListRef->end(), this->productTree, firstTreeOutput);
  this->clearProducts();
}

void BaseTreeLooper::loop(bool keepProducts){
  if (!looperFunction){
    MELAerr << "BaseTreeLooper::loop: The looper function is not registered. Please register it using BadeTreeLooper::setLooperFunction." << endl;
    return;
  }

  // Loop over the trees
  unsigned int ev_acc=0;
  unsigned int ev_rec=0;
  for (auto& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->linkConsumes(tree))) continue;

    float globalTreeWeight = 1;
    auto it_globalWgt = globalWeights.find(tree);
    if (it_globalWgt!=globalWeights.cend()) globalTreeWeight = it_globalWgt->second;

    MELAout << "BaseTreeLooper::loop: Looping over " << tree->sampleIdentifier << " events" << endl;
    int ev=0;
    const int nevents = tree->getNEvents();
    while (tree->getEvent(ev)){
      if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;
      SimpleEntry product;
      if (tree->isValidEvent()){
        if (this->looperFunction(this, tree, globalTreeWeight, product)){
          if (keepProducts) this->addProduct(product, &ev_rec);
        }
      }
      HelperFunctions::progressbar(ev, nevents);
      ev++; ev_acc++;
    }

    // Record products to external tree
    this->recordProductsToTree();
  } // End loop over the trees
  MELAout << "BaseTreeLooper::loop: Total number of products: " << ev_rec << " / " << ev_acc << endl;
}

std::vector<SimpleEntry> const& BaseTreeLooper::getProducts() const{ return *productListRef; }

void BaseTreeLooper::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "BaseTreeLooper::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "BaseTreeLooper::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void BaseTreeLooper::clearProducts(){ std::vector<SimpleEntry> emptyList; std::swap(emptyList, *productListRef); }
