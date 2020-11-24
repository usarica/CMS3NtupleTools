#ifndef BASETREELOOPER_H
#define BASETREELOOPER_H

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include "IvyBase.h"
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "BulkReweightingBuilder.h"
#include "TriggerHelpersCore.h"
#include "ParticleDisambiguator.h"
#include "DileptonHandler.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


class BaseTreeLooper : public IvyBase{
public:
  typedef bool(*LooperCoreFunction_t)(BaseTreeLooper*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const&, SimpleEntry&);
  typedef void(*LooperExtFunction_t)(BaseTreeLooper*, SimpleEntry&);

protected:
  enum SampleIdStorageType{
    kNoStorage,
    kStoreByRunAndEventNumber,
    kStoreByMH
  };

  // Function to determine if event should be included
  LooperCoreFunction_t looperFunction;

  // Systematics type
  SystematicsHelpers::SystematicVariationTypes registeredSyst;

  // Max. events to process
  int maxNEvents;

  // Event index ranges
  int eventIndex_begin;
  int eventIndex_end;
  bool useChunkIndices;

  // Variables set per tree
  bool isData_currentTree;
  bool isQCD_currentTree;
  bool isGJets_HT_currentTree;
  float pTG_true_exception_range[2];

  // Some ready-made stuff
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  // ME lists
  std::vector<std::string> lheMElist;
  std::vector<std::string> recoMElist;
  CMS3MELAHelpers::GMECBlock MEblock;

  // Selection counts
  std::vector<std::pair<TString, unsigned int>> selection_string_count_pairs;

  // Input trees
  std::vector<BaseTree*> treeList;

  // Registered object handlers
  std::vector<IvyBase*> registeredHandlers;

  // Registered SF handlers
  std::vector<ScaleFactorHandlerBase*> registeredSFHandlers;

  // Registered reweighting builders
  std::unordered_map<TString, BulkReweightingBuilder*> registeredRewgtBuilders;

  // Registered triggers
  std::unordered_map<TString, std::vector< std::string > > registeredHLTMenus;
  std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > > registeredHLTMenuProperties;

  // External dependencies
  std::unordered_map<BaseTree*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double>> globalWeights;
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

  void sigint_callback_handler(int snum);

  void resetSelectionCounts(){ selection_string_count_pairs.clear(); }

public:
  // Constructors
  BaseTreeLooper();
  BaseTreeLooper(BaseTree* inTree, double wgt=1);
  BaseTreeLooper(std::vector<BaseTree*> const& inTreeList);

  // Destructors
  virtual ~BaseTreeLooper();

  // Add the necessary objects
  void addTree(BaseTree* tree, double wgt=1); // Adds a new input tree
  void addTree(BaseTree* tree, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& syst_wgt_map); // Adds a new input tree
  void addExternalFunction(TString fcnname, BaseTreeLooper::LooperExtFunction_t fcn);
  void addObjectHandler(IvyBase* handler);
  void addSFHandler(ScaleFactorHandlerBase* handler);
  void addReweightingBuilder(TString name, BulkReweightingBuilder* rewgtBuilder);
  void addHLTMenu(TString name, std::vector< std::string > const& hltmenu);
  void addHLTMenu(TString name, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > const& hltmenu);

  void setLooperFunction(BaseTreeLooper::LooperCoreFunction_t fcn){ looperFunction = fcn; }
  void setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){ registeredSyst = syst; }
  void setExternalWeight(BaseTree* tree, double const& wgt);
  void setExternalWeights(BaseTree* tree, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& wgts);
  void setMatrixElementList(std::vector<std::string> const& MElist, bool const& isGen);
  void setMatrixElementListFromFile(std::string fname, std::string const& MElistTypes, bool const& isGen); // MElistTypes is comman-separated

  void setExternalProductList(std::vector<SimpleEntry>* extProductListRef=nullptr);
  void setCurrentOutputTree(BaseTree* extTree=nullptr);
  void addOutputTree(BaseTree* extTree);
  void addOutputTrees(std::vector<BaseTree*> trees);

  // Max. events
  void setMaximumEvents(int n);

  // Event index range
  void setEventIndexRange(int istart, int iend);
  void setEventIndexRangeBySampleChunks(bool flag){ useChunkIndices = flag; } // Set if the function above uses chunk indices instead of actual event ranges

  // Get-functions
  int const& getMaximumEvents() const{ return maxNEvents; }
  bool const& getCurrentTreeFlag_IsData() const{ return isData_currentTree; }
  bool const& getCurrentTreeFlag_QCDException() const{ return isQCD_currentTree; }
  bool const& getCurrentTreeFlag_GJetsHTException() const{ return isGJets_HT_currentTree; }
  bool getPTGExceptionRange(float& vlow, float& vhigh) const{ vlow = pTG_true_exception_range[0]; vhigh = pTG_true_exception_range[1]; return (vlow!=vhigh); }
  SystematicsHelpers::SystematicVariationTypes const& getSystematic() const{ return registeredSyst; }
  std::vector<IvyBase*> const& getObjectHandlers() const{ return registeredHandlers; }
  std::vector<ScaleFactorHandlerBase*> const& getSFHandlers() const{ return registeredSFHandlers; }
  std::unordered_map<TString, BulkReweightingBuilder*> const& getRegisteredRewgtBuilders() const{ return registeredRewgtBuilders; }
  std::unordered_map<TString, std::vector< std::string > > const& getHLTMenus() const{ return registeredHLTMenus; }
  std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > > const& getHLTMenuProperties() const{ return registeredHLTMenuProperties; }
  ParticleDisambiguator& getParticleDisambiguator(){ return particleDisambiguator; }
  ParticleDisambiguator const& getParticleDisambiguator() const{ return particleDisambiguator; }
  DileptonHandler& getDileptonHandler(){ return dileptonHandler; }
  DileptonHandler const& getDileptonHandler() const{ return dileptonHandler; }
  CMS3MELAHelpers::GMECBlock& getMEblock(){ return MEblock; }
  CMS3MELAHelpers::GMECBlock const& getMEblock() const{ return MEblock; }

  bool hasSimpleHLTMenus() const{ return registeredHLTMenus.size()>0; }
  bool hasHLTMenuProperties() const{ return registeredHLTMenuProperties.size()>0; }

  bool hasLinkedOutputTrees() const{ return !productTreeList.empty(); }
  bool hasGenMEs() const{ return !lheMElist.empty(); }
  bool hasRecoMEs() const{ return !recoMElist.empty(); }

  // Function to add and increment a selection type to count
  void incrementSelection(TString const& strsel, unsigned int inc=1);

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
