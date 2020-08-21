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
#include "ParticleDisambiguator.h"
#include "DileptonHandler.h"


class BaseTreeLooper : public IvyBase{
public:
  typedef bool(*LooperCoreFunction_t)(BaseTreeLooper*, float const&, SimpleEntry&);
  typedef void(*LooperExtFunction_t)(BaseTreeLooper*, SimpleEntry&);

protected:
  // Function to determine if event should be included
  LooperCoreFunction_t looperFunction;

  // Systematics type
  SystematicsHelpers::SystematicVariationTypes registeredSyst;

  // Max. events to process
  int maxNEvents;

  // Event index ranges
  int eventIndex_begin;
  int eventIndex_end;

  // Variables set per tree
  bool isData_currentTree;
  bool isQCD_currentTree;
  float pTG_true_exception_range[2];

  // Some ready-made stuff
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

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
  std::unordered_map<TString, LooperExtFunction_t> externalFunctions;

  // Output trees
  std::vector<BaseTree*> productTreeList;

  // Flags for output trees
  std::unordered_map<BaseTree*, bool> firstTreeOutput;

  // List of products
  std::vector<SimpleEntry> productList;
  std::vector<SimpleEntry>* productListRef;
  BaseTree* currentProductTree;
  void addProduct(SimpleEntry& product, unsigned int* ev_rec=nullptr);

  // Set pTG exception range
  void set_pTG_exception_range(float const& vlow, float const& vhigh){ pTG_true_exception_range[0] = vlow; pTG_true_exception_range[1] = vhigh; }

  // Flush product list into tree
  void recordProductsToTree();

  bool wrapTree(BaseTree* tree);

public:
  // Constructors
  BaseTreeLooper();
  BaseTreeLooper(BaseTree* inTree);
  BaseTreeLooper(std::vector<BaseTree*> const& inTreeList);
  void addTree(BaseTree* tree);

  // Destructors
  virtual ~BaseTreeLooper();

  // Add the necessary objects
  void addExternalFunction(TString fcnname, BaseTreeLooper::LooperExtFunction_t fcn);
  void addObjectHandler(IvyBase* handler);
  void addSFHandler(ScaleFactorHandlerBase* handler);

  void setLooperFunction(BaseTreeLooper::LooperCoreFunction_t fcn){ looperFunction = fcn; }
  void setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){ registeredSyst = syst; }

  void setExternalProductList(std::vector<SimpleEntry>* extProductListRef=nullptr);
  void setCurrentOutputTree(BaseTree* extTree=nullptr);
  void addOutputTree(BaseTree* extTree);
  void addOutputTrees(std::vector<BaseTree*> trees);

  // Max. events
  void setMaximumEvents(int n);

  // Event index offset
  void setEventIndexRange(int istart, int iend);

  // Get-functions
  int const& getMaximumEvents() const{ return maxNEvents; }
  bool const& getCurrentTreeFlag_IsData() const{ return isData_currentTree; }
  bool const& getCurrentTreeFlag_QCDException() const{ return isQCD_currentTree; }
  bool getPTGExceptionRange(float& vlow, float& vhigh) const{ vlow = pTG_true_exception_range[0]; vhigh = pTG_true_exception_range[1]; return (vlow!=vhigh); }
  SystematicsHelpers::SystematicVariationTypes const& getSystematic() const{ return registeredSyst; }
  std::vector<IvyBase*> const& getObjectHandlers() const{ return registeredHandlers; }
  std::vector<ScaleFactorHandlerBase*> const& getSFHandlers() const{ return registeredSFHandlers; }
  std::unordered_map<TString, std::vector< std::string > > const& getHLTMenus() const{ return registeredHLTMenus; }
  std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > > const& getHLTMenuProperties() const{ return registeredHLTMenuProperties; }
  ParticleDisambiguator& getParticleDisambiguator(){ return particleDisambiguator; }
  ParticleDisambiguator const& getParticleDisambiguator() const{ return particleDisambiguator; }
  DileptonHandler& getDileptonHandler(){ return dileptonHandler; }
  DileptonHandler const& getDileptonHandler() const{ return dileptonHandler; }

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
