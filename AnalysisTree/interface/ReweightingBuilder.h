#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "BaseTree.h"
#include "ExtendedBinning.h"
#include "ReweightingFunctions.h"


class ReweightingBuilder{
public:
  constexpr static bool useNeffInNormComponent=true;

protected:
  bool allowNegativeWeights;
  float weightThresholdReference;

  std::vector<TString> const strNominalWeights;
  std::vector<TString> const strReweightingWeights;
  std::vector<TString> const strCrossSectionWeights;

  ReweightingFunctions::ReweightingFunction_t rule_nominalweights;
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights;
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights;

  std::unordered_map<BaseTree*, std::vector<float*>> componentRefs_nominalweights;
  std::unordered_map<BaseTree*, std::vector<float*>> componentRefs_reweightingweights;
  std::unordered_map<BaseTree*, std::vector<float*>> componentRefs_xsecweights;

  ExtendedBinning weightBinning;
  std::unordered_map<BaseTree*, float*> binningVarRefs;

  std::unordered_map<BaseTree*, std::vector<float>> weightThresholds;
  std::unordered_map<BaseTree*, std::vector<float>> sumPostThrWeights;
  std::unordered_map<BaseTree*, std::vector<float>> sumPostThrSqWeights;
  std::unordered_map<BaseTree*, std::vector<unsigned int>> sumNonZeroWgtEvents;
  std::unordered_map<BaseTree*, std::vector<float>> sumNonZeroWgtNominalWeights;

  std::unordered_map<BaseTree*, std::vector<float>> cachedSampleNormalizationsPerBin;
  std::vector<float> cachedNormComponentsPerBin;

  std::unordered_map<BaseTree*, float> xsecVals;
  std::unordered_map<BaseTree*, float> sumNominalWeights;

  int findBin(BaseTree* theTree) const;

public:
  ReweightingBuilder(
    std::vector<TString> const& strNominalWeights_,
    std::vector<TString> const& strReweightingWeights_,
    std::vector<TString> const& strCrossSectionWeights_,
    ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
    ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_,
    ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
  );
  
  virtual ~ReweightingBuilder(){}

  virtual float eval_nominalweights(BaseTree* theTree) const;
  virtual float eval_reweightingweights(BaseTree* theTree) const;
  virtual float eval_xsecweights(BaseTree* theTree) const;

  std::vector<float> getWeightThresholds(BaseTree* theTree) const;
  float getPostThresholdWeight(BaseTree* theTree) const;
  float getFinalEventWeight(BaseTree* theTree) const;

  void rejectNegativeWeights(const bool flag);

  void setWeightBinning(const ExtendedBinning& binning);
  void setWeightThresholdReference(const float& weightThresholdReference_);
  void setupWeightVariables(BaseTree* theTree, float fractionRequirement=0.999, unsigned int minimumNevents=0);
  void setupCaches();

  std::vector<BaseTree*> getRegisteredTrees() const;

  float getNormComponent(BaseTree* theTree) const; // Tree is passed here to find the bin
  float getNormComponent(int bin) const;

  float getBareSumNominalWeights(BaseTree* theTree) const;

};


#endif
