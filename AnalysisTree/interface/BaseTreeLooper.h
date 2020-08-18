#ifndef BASETREELOOPER_H
#define BASETREELOOPER_H

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include "IvyBase.h"
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "TriggerHelpersCore.h"


class BaseTreeLooper : IvyBase{
protected:
  // Function to determine if event should be included
  bool(*looperFunction)(BaseTreeLooper const*, BaseTree*, float const&, SimpleEntry&);

  // Systematics type
  SystematicsHelpers::SystematicVariationTypes registeredSyst;

  // Max. events to process
  int maxNEvents;

  // Event index ranges
  int eventIndex_begin;
  int eventIndex_end;

  // Flag for output tree
  bool firstTreeOutput;

  // Input trees
  std::vector<BaseTree*> treeList;

  // Registered object handlers
  std::vector<IvyBase*> registeredHandlers;

  // Registered SF handlers
  std::vector<ScaleFactorHandlerBase*> registeredSFHandlers;

  // Registered triggers
  std::unordered_map<TString, std::vector< std::string > > registeredHLTMenus;
  std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > > registeredHLTMenuProperties;

  // External dependencies
  std::unordered_map<BaseTree*, float> globalWeights;
  std::unordered_map<TString, void(*)(BaseTreeLooper const*, BaseTree*, SimpleEntry&)> externalFunctions;

  // List of products
  std::vector<SimpleEntry> productList;
  std::vector<SimpleEntry>* productListRef;
  BaseTree* productTree;
  void addProduct(SimpleEntry& product, unsigned int* ev_rec=nullptr);

  // Flush product list into tree
  void recordProductsToTree();

public:
  // Constructors
  BaseTreeLooper();
  BaseTreeLooper(BaseTree* inTree);
  BaseTreeLooper(std::vector<BaseTree*> const& inTreeList);
  void addTree(BaseTree* tree);

  // Destructors
  virtual ~BaseTreeLooper();

  // Add the necessary objects
  void addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper const*, BaseTree*, SimpleEntry&));
  void addObjectHandler(IvyBase* handler);
  void addSFHandler(ScaleFactorHandlerBase* handler);

  void setLooperFunction(bool(*fcn)(BaseTreeLooper const*, BaseTree*, float const&, SimpleEntry&)){ looperFunction = fcn; }
  void setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){ registeredSyst = syst; }

  void setExternalProductList(std::vector<SimpleEntry>* extProductListRef=nullptr);
  void setExternalProductTree(BaseTree* extTree=nullptr);

  // Max. events
  void setMaximumEvents(int n);

  // Event index offset
  void setEventIndexRange(int istart, int iend);

  // Get-functions
  int const& getMaximumEvents() const{ return maxNEvents; }
  SystematicsHelpers::SystematicVariationTypes const& getSystematic() const{ return registeredSyst; }
  std::vector<IvyBase*> const& getObjectHandlers() const{ return registeredHandlers; }
  std::vector<ScaleFactorHandlerBase*> const& getSFHandlers() const{ return registeredSFHandlers; }
  std::unordered_map<TString, std::vector< std::string > > const& getHLTMenus(){ return registeredHLTMenus; }
  std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > > const& getHLTMenuProperties(){ return registeredHLTMenuProperties; }

  bool hasSimpleHLTMenus() const{ return registeredHLTMenus.size()>0; }
  bool hasHLTMenuProperties() const{ return registeredHLTMenuProperties.size()>0; }

  // Function to loop over the tree list
  virtual void loop(bool keepProducts);

  // Get the products
  std::vector<SimpleEntry> const& getProducts() const;
  // Move the products
  void moveProducts(std::vector<SimpleEntry>& targetColl);
  // Clear the products
  void clearProducts();

};


#endif
