#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "BaseTree.h"
#include "ExtendedBinning.h"
#include "ReweightingFunctions.h"


class BulkReweightingBuilder{
protected:
  std::vector<BaseTree*> registeredTrees;

  ExtendedBinning binning;
  std::vector<TString> const strBinningVars;
  ReweightingFunctions::ReweightingVariableBinFunction_t rule_binningVar;
  std::unordered_map< BaseTree*, std::vector<float*> > binningVarRefs;

  std::vector<TString> const strNominalWeights;
  ReweightingFunctions::ReweightingFunction_t rule_nominalweights;
  std::unordered_map< BaseTree*, std::vector<float*> > componentRefs_nominalweights;

  std::vector<TString> const strCrossSectionWeights;
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights;
  std::unordered_map< BaseTree*, std::vector<float*> > componentRefs_xsecweights;

  std::vector<std::vector<TString>> strReweightingWeightsList;
  std::vector<ReweightingFunctions::ReweightingFunction_t> rule_reweightingweights_list;
  std::unordered_map< BaseTree*, std::vector<std::vector<float*> > > componentRefsList_reweightingweights;

  std::unordered_map<BaseTree*, float> normWeights;

  // Derived variables
  std::unordered_map< BaseTree*, std::vector< std::vector<float> > > absWeightThresholdsPerBinList;
  std::unordered_map< BaseTree*, std::vector<double> > sum_normwgts_all;
  std::unordered_map< BaseTree*, std::vector<double> > sum_normwgts_nonzerorewgt;
  std::unordered_map< BaseTree*, std::vector< std::vector< std::pair<double, double> > > > sum_wgts_withrewgt;
  std::unordered_map< BaseTree*, std::vector<double> > NeffsPerBin;
  std::unordered_map< BaseTree*, std::vector<double> > sampleNormalization;

public:
  BulkReweightingBuilder(
    ExtendedBinning const& binning_,
    std::vector<TString> const& strBinningVars_,
    std::vector<TString> const& strNominalWeights_,
    std::vector<TString> const& strCrossSectionWeights_,
    ReweightingFunctions::ReweightingVariableBinFunction_t rule_binningVar_,
    ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
    ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
  );
  virtual ~BulkReweightingBuilder(){}

  void addReweightingWeights(
    std::vector<TString> const& strReweightingWeights_,
    ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_
  );
  void registerTree(BaseTree* tree, float extNorm=1);

  void setup(unsigned int ihypo_Neff, std::vector<std::pair<BaseTree*, BaseTree*>> const* tree_normTree_pairs=nullptr);

  double getOverallReweightingNormalization(BaseTree* tree) const;

  bool checkWeightsBelowThreshold(BaseTree* tree) const;

};


#endif
