#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "BaseTree.h"
#include "ExtendedBinning.h"
#include "ReweightingFunctions.h"


struct ReweightingSpecifications{
  bool allowNegativeWeights;
  float weightThresholdReference;

  std::vector<TString> strReweightingWeights;
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights;
  std::unordered_map<BaseTree*, std::vector<float*>> componentRefs_reweightingweights;

  std::unordered_map<BaseTree*, std::vector<float>> weightThresholds;
  std::unordered_map<BaseTree*, std::vector<float>> sumPostThrWeights;
  std::unordered_map<BaseTree*, std::vector<float>> sumPostThrSqWeights;
  std::unordered_map<BaseTree*, std::vector<unsigned int>> sumNonZeroWgtEvents;
  std::unordered_map<BaseTree*, std::vector<float>> sumNonZeroWgtNominalWeights;

  std::unordered_map<BaseTree*, std::vector<float>> cachedSampleNormalizationsPerBin;
  std::vector<float> cachedNormComponentsPerBin;

  std::vector<std::vector<SimpleEntry>> indexList;

  ReweightingSpecifications();
  ReweightingSpecifications(
    bool allowNegativeWeights_,
    std::vector<TString> const& strReweightingWeights_,
    ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_
  );
  ReweightingSpecifications(ReweightingSpecifications const& other);

  void swap(ReweightingSpecifications& other);
  ReweightingSpecifications& operator=(const ReweightingSpecifications& other);

  void initialize(BaseTree* theTree, unsigned int ns, bool doReweighting);
  void setupThresholds(BaseTree* theTree, float fractionRequirement, unsigned int minimumNevents);

  float eval_reweightingweights(BaseTree* theTree) const;

};


class BulkReweightingBuilder{
public:
  constexpr static bool useNeffInNormComponent=true;

protected:
  bool allowNegativeWeights;

  ExtendedBinning weightBinning;
  std::unordered_map<BaseTree*, float*> binningVarRefs;

  std::vector<TString> const strNominalWeights;
  ReweightingFunctions::ReweightingFunction_t rule_nominalweights;
  std::unordered_map<BaseTree*, std::vector<float*>> componentRefs_nominalweights;
  std::unordered_map<BaseTree*, float> sumNominalWeights;

  std::vector<TString> const strCrossSectionWeights;
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights;
  std::unordered_map<BaseTree*, std::vector<float*>> componentRefs_xsecweights;
  std::unordered_map<BaseTree*, float> xsecVals;

  std::unordered_map<TString, ReweightingSpecifications> reweightingSpecs;
  ReweightingSpecifications* currentReweightingSpecs;

  int findBin(BaseTree* theTree) const;
  bool setCurrentReweightingSpecs(TString strScheme);
  ReweightingSpecifications* getReweightingSpecs(TString strScheme);
  ReweightingSpecifications const* getReweightingSpecs(TString strScheme) const;

public:
  BulkReweightingBuilder(
    std::vector<TString> const& strNominalWeights_,
    std::vector<TString> const& strCrossSectionWeights_,
    ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
    ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
  );
  virtual ~BulkReweightingBuilder(){}

  void addReweightingWeights(
    TString strScheme,
    std::vector<TString> const& strReweightingWeights_,
    ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_
  );

  virtual float eval_nominalweights(BaseTree* theTree) const;
  virtual float eval_reweightingweights(BaseTree* theTree, TString strScheme) const;
  virtual float eval_xsecweights(BaseTree* theTree) const;

  std::vector<float> getWeightThresholds(BaseTree* theTree, TString strScheme) const;
  float getPostThresholdWeight(BaseTree* theTree, TString strScheme) const;
  float getFinalEventWeight(BaseTree* theTree, TString strScheme) const;

  void rejectNegativeWeights(const bool flag);

  void setWeightBinning(const ExtendedBinning& binning);
  void setWeightThresholdReference(const float& weightThresholdReference_, TString strScheme);

  void setupWeightVariables(BaseTree* theTree, float fractionRequirement=0.999, unsigned int minimumNevents=0);
  void setupCaches();

  std::vector<BaseTree*> getRegisteredTrees() const;

  float getNormComponent(BaseTree* theTree, TString strScheme) const; // Tree is passed here to find the bin
  float getNormComponent(int bin, TString strScheme) const;

  float getBareSumNominalWeights(BaseTree* theTree) const;

};


#endif
