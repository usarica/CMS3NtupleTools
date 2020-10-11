#include "BulkReweightingBuilder.h"
#include "SimpleEntry.h"
#include "MELAAccumulators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


BulkReweightingBuilder::BulkReweightingBuilder(
  ExtendedBinning const& binning_,
  std::vector<TString> const& strBinningVars_,
  std::vector<TString> const& strNominalWeights_,
  std::vector<TString> const& strCrossSectionWeights_,
  ReweightingFunctions::ReweightingVariableBinFunction_t rule_binningVar_,
  ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
) :
  binning(binning_),
  strBinningVars(strBinningVars_),
  rule_binningVar(rule_binningVar_),

  strNominalWeights(strNominalWeights_),
  rule_nominalweights(rule_nominalweights_),

  strCrossSectionWeights(strCrossSectionWeights_),
  rule_xsecweights(rule_xsecweights_)
{}

void BulkReweightingBuilder::addReweightingWeights(
  std::vector<TString> const& strReweightingWeights_,
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_
){
  strReweightingWeightsList.push_back(strReweightingWeights_);
  rule_reweightingweights_list.push_back(rule_reweightingweights_);
}

void BulkReweightingBuilder::registerTree(BaseTree* tree, float extNorm){
  for (auto const& s:strBinningVars) tree->bookBranch<float>(s, 0.);
  binningVarRefs[tree] = ReweightingFunctions::getWeightRefs(tree, strBinningVars);

  for (auto const& s:strNominalWeights) tree->bookBranch<float>(s, 0.);
  componentRefs_nominalweights[tree] = ReweightingFunctions::getWeightRefs(tree, strNominalWeights);

  for (auto const& s:strCrossSectionWeights) tree->bookBranch<float>(s, 0.);
  componentRefs_xsecweights[tree] = ReweightingFunctions::getWeightRefs(tree, strCrossSectionWeights);

  componentRefsList_reweightingweights[tree] = std::vector<std::vector<float*>>(strReweightingWeightsList.size(), std::vector<float*>());
  for (unsigned int ihypo=0; ihypo<strReweightingWeightsList.size(); ihypo++){
    for (auto const& s:strReweightingWeightsList.at(ihypo)) tree->bookBranch<float>(s, 0.);
    componentRefsList_reweightingweights[tree].at(ihypo) = ReweightingFunctions::getWeightRefs(tree, strReweightingWeightsList.at(ihypo));
  }

  normWeights[tree] = extNorm;
  registeredTrees.push_back(tree);
}

void BulkReweightingBuilder::setup(unsigned int ihypo_Neff, std::vector<std::pair<BaseTree*, BaseTree*>> const* tree_normTree_pairs){
  unsigned int const nhypos = strReweightingWeightsList.size();
  assert(ihypo_Neff<nhypos);
  unsigned int const nbins = (binning.isValid() ? binning.getNbins() : static_cast<unsigned int>(1));
  for (auto const& tree:registeredTrees){
    MELAout << "BulkReweightingBuilder::setup: Processing " << tree->sampleIdentifier << "..." << endl;

    MELAout << "\t- Obtaining Neff thresholds..." << endl;
    NeffThrsPerBin[tree] = ReweightingFunctions::getSimpleNeffThrsPerBin(tree, binning, binningVarRefs[tree], rule_binningVar, -1);

    MELAout << "\t- Obtaining weight thresholds..." << endl;
    absWeightThresholdsPerBinList[tree] = std::vector<std::vector<float>>(nhypos, std::vector<float>());
    for (unsigned int ihypo=0; ihypo<nhypos; ihypo++) absWeightThresholdsPerBinList[tree].at(ihypo) = ReweightingFunctions::getAbsWeightThresholdsPerBinByNeff(
      tree,
      componentRefsList_reweightingweights[tree].at(ihypo), rule_reweightingweights_list.at(ihypo),
      binning, binningVarRefs[tree], rule_binningVar,
      NeffThrsPerBin[tree]
    );

    NeffsPerBin[tree] = std::vector<double>(nbins, 0);

    // Get sum of weights
    int nEntries = tree->getNEvents();
    std::vector<double> sum_normwgts_all_tree(nbins, 0);
    std::vector<double> sum_normwgts_nonzerorewgt_tree(nbins, 0);
    std::vector<std::vector<std::pair<double, double>>> sum_wgts_withrewgt_tree(nhypos, std::vector<std::pair<double, double>>(nbins, std::pair<double, double>(0, 0)));
    MELAout << "\t- Obtaining weight sums..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      tree->getEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);

      int ibin = rule_binningVar(tree, binning, binningVarRefs[tree]);
      if (ibin<0) ibin=0;
      else if (ibin>=(int) nbins) ibin = nbins-1;

      double wgt_nominal = rule_nominalweights(tree, componentRefs_nominalweights[tree]) * normWeights[tree];
      double wgt_xsec = rule_xsecweights(tree, componentRefs_xsecweights[tree]);
      double wgt_nominal_xsec = wgt_nominal*wgt_xsec;
      sum_normwgts_all_tree.at(ibin) += wgt_nominal_xsec;
      bool allHyposFine = true;
      std::vector<double> wgt_rewgt_list; wgt_rewgt_list.reserve(nhypos);
      for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
        float const& wgt_thr = absWeightThresholdsPerBinList[tree].at(ihypo).at(ibin);
        float wgt_rewgt = rule_reweightingweights_list.at(ihypo)(tree, componentRefsList_reweightingweights[tree].at(ihypo));
        if (wgt_thr>0.f && std::abs(wgt_rewgt)>wgt_thr) wgt_rewgt = 0.f;
        if (wgt_rewgt==0.f) allHyposFine = false;
        wgt_rewgt_list.push_back(wgt_rewgt);
      }
      if (allHyposFine){
        for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
          sum_normwgts_nonzerorewgt_tree.at(ibin) += wgt_nominal_xsec;
          auto& sum_wgts_pair = sum_wgts_withrewgt_tree.at(ihypo).at(ibin);
          double wgt_product = wgt_nominal_xsec*wgt_rewgt_list.at(ihypo);
          sum_wgts_pair.first += wgt_product;
          sum_wgts_pair.second += std::pow(wgt_product, 2);
        }
      }
    }

    // Assign the vectors to the maps
    sum_normwgts_all[tree] = sum_normwgts_all_tree;
    sum_normwgts_nonzerorewgt[tree] = sum_normwgts_nonzerorewgt_tree;
    sum_wgts_withrewgt[tree] = sum_wgts_withrewgt_tree;

    // Will be filled after all trees are done, but there is no harm in assigning here.
    sampleNormalization[tree] = std::vector<double>(nbins, 1);

    // Compute Neff to be able to combine the trees ultimately
    for (unsigned int ibin=0; ibin<nbins; ibin++){
      auto& sum_wgts_pair = sum_wgts_withrewgt_tree.at(ihypo_Neff).at(ibin);
      if (sum_wgts_pair.second==0.) continue;
      NeffsPerBin[tree].at(ibin) = std::pow(sum_wgts_pair.first, 2)/sum_wgts_pair.second;
    }
  }

  // Once all trees are complete, compute final normalization factors
  for (auto const& tree:registeredTrees){
    for (unsigned int ibin=0; ibin<nbins; ibin++){
      double& sampleNormalization_bin = sampleNormalization[tree].at(ibin);
      double sum_Neff_bin = 0;
      for (auto const& it:NeffsPerBin) sum_Neff_bin += it.second.at(ibin);
      if (sum_Neff_bin>0.) sampleNormalization_bin = NeffsPerBin[tree].at(ibin) / sum_Neff_bin;
    }
    if (tree_normTree_pairs){
      for (auto const& tree_normTree_pair:(*tree_normTree_pairs)){
        if (tree_normTree_pair.first == tree){
          auto const& normTree = tree_normTree_pair.second;
          auto const& sum_wgts_withrewgt_basetree = sum_wgts_withrewgt[tree].at(ihypo_Neff);
          auto const& sum_wgts_withrewgt_normtree = sum_wgts_withrewgt[normTree].at(ihypo_Neff);

          double sum_wgts_common_basetree = 0;
          double sum_wgts_common_normtree = 0;
          for (unsigned int ibin=0; ibin<nbins; ibin++){
            double const& sum_wgts_common_basetree_bin = sum_wgts_withrewgt_basetree.at(ibin).first;
            double const& sum_wgts_common_normtree_bin = sum_wgts_withrewgt_normtree.at(ibin).first;
            if (sum_wgts_common_basetree_bin!=0. && sum_wgts_common_normtree_bin!=0.){
              sum_wgts_common_basetree += sum_wgts_common_basetree_bin;
              sum_wgts_common_normtree += sum_wgts_common_normtree_bin;
            }
          }
          if (sum_wgts_common_basetree==0.) MELAerr << "BulkReweightingBuilder::setup: Base tree " << tree->sampleIdentifier << " has no overlap with its norm tree " << normTree->sampleIdentifier << endl;
          if (sum_wgts_common_normtree==0.) MELAerr << "BulkReweightingBuilder::setup: Norm tree " << tree->sampleIdentifier << " has no overlap with its base tree " << normTree->sampleIdentifier << endl;

          double extra_norm = 1;
          if (sum_wgts_common_basetree!=0.) extra_norm = sum_wgts_common_normtree / sum_wgts_common_basetree;
          MELAout << "\t- Base tree " << tree->sampleIdentifier << " has an extra normalization of " << extra_norm << endl;
          for (unsigned int ibin=0; ibin<nbins; ibin++) sampleNormalization[tree].at(ibin) *= extra_norm;
        }
      }
    }
  }
}
double BulkReweightingBuilder::getOverallReweightingNormalization(BaseTree* tree) const{
  auto it_binningVarRefs = binningVarRefs.find(tree);
  if (it_binningVarRefs == binningVarRefs.cend()){
    MELAerr << "BulkReweightingBuilder::getOverallReweightingNormalization: Tree " << tree->sampleIdentifier << " is not registered properly." << endl;
    return 1;
  }

  unsigned int const nbins = (binning.isValid() ? binning.getNbins() : static_cast<unsigned int>(1));
  int ibin = rule_binningVar(tree, binning, it_binningVarRefs->second);
  if (ibin<0) ibin=0;
  else if (ibin>=(int) nbins) ibin = nbins-1;

  auto it = sampleNormalization.find(tree);
  if (it!=sampleNormalization.cend()) return it->second.at(ibin);
  else{
    MELAerr << "BulkReweightingBuilder::getOverallReweightingNormalization: No normalization factor is found for tree " << tree->sampleIdentifier << "." << endl;
    return 1;
  }
}

