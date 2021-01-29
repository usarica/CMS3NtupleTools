#include "BulkReweightingBuilder.h"
#include "HostHelpersCore.h"
#include "SimpleEntry.h"
#include "ExtendedHistogram_1D.h"
#include "ExtendedHistogram_2D.h"
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
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_,
  double thr_wgt, double tolerance
){
  strReweightingWeightsList.push_back(strReweightingWeights_);
  rule_reweightingweights_list.push_back(rule_reweightingweights_);
  reweightingweights_frac_tolerance_pair_list.emplace_back(thr_wgt, tolerance);
}

void BulkReweightingBuilder::registerTree(BaseTree* tree, double extNorm){
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

void BulkReweightingBuilder::setup(
  int ihypo_Neff, std::vector<std::pair<BaseTree*, BaseTree*>> const* tree_normTree_pairs,
  float thr_frac_Neff
){
  unsigned int const nhypos = strReweightingWeightsList.size();
  assert(ihypo_Neff<(int) nhypos);
  assert(thr_frac_Neff<=1.f);
  unsigned int const nbins = (binning.isValid() ? binning.getNbins() : static_cast<unsigned int>(1));
  for (auto const& tree:registeredTrees){
    MELAout << "BulkReweightingBuilder::setup: Processing " << tree->sampleIdentifier << "..." << endl;

    // Mute all branches except the ones of interest in order to speed up the process
    {
      std::vector<TString> bnames_preserved;
      for (auto const& strReweightingWeights:strReweightingWeightsList){
        for (auto const& bname:strReweightingWeights){
          if (!HelperFunctions::checkListVariable(bnames_preserved, bname)) bnames_preserved.push_back(bname);
        }
      }
      for (auto const& bname:strBinningVars){
        if (!HelperFunctions::checkListVariable(bnames_preserved, bname)) bnames_preserved.push_back(bname);
      }
      for (auto const& bname:strNominalWeights){
        if (!HelperFunctions::checkListVariable(bnames_preserved, bname)) bnames_preserved.push_back(bname);
      }
      for (auto const& bname:strCrossSectionWeights){
        if (!HelperFunctions::checkListVariable(bnames_preserved, bname)) bnames_preserved.push_back(bname);
      }
      tree->muteAllBranchesExcept(bnames_preserved);
    }

    MELAout << "\t- Obtaining weight thresholds..." << endl;
    absWeightThresholdsPerBinList[tree] = ReweightingFunctions::getAbsWeightThresholdsPerBinByFixedFractionalThreshold(
      tree,
      componentRefsList_reweightingweights[tree], rule_reweightingweights_list, reweightingweights_frac_tolerance_pair_list,
      binning, binningVarRefs[tree], rule_binningVar
    );

    // Initialize normalization variables
    sampleNormalization[tree] = std::vector<double>(nbins, 1);
    samplePairwiseNormalization[tree] = 1;

    // Get sum of weights
    int nEntries = tree->getNEvents();
    std::vector<double> sum_normwgts_all_tree(nbins, 0);
    std::vector<double> sum_normwgts_nonzerorewgt_tree(nbins, 0);
    std::vector<double> NeffsPerBin_tree(nbins, 0);
    std::vector<double> sampleZeroMECompensation_tree(nbins, 1);
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
        if (wgt_rewgt==0.f){
          allHyposFine = false;
          break; // We can break because the following statement only proceeds if allHyposFine==true
        }
        wgt_rewgt_list.push_back(wgt_rewgt);
      }
      if (allHyposFine){
        sum_normwgts_nonzerorewgt_tree.at(ibin) += wgt_nominal_xsec;

        for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
          auto& sum_wgts_pair = sum_wgts_withrewgt_tree.at(ihypo).at(ibin);
          double wgt_product = wgt_nominal_xsec*wgt_rewgt_list.at(ihypo);
          sum_wgts_pair.first += wgt_product;
          sum_wgts_pair.second += std::pow(wgt_product, 2);
        }
      }
    }

    // Compute Neff to be able to combine the trees ultimately
    for (unsigned int ibin=0; ibin<nbins; ibin++){
      if (sum_normwgts_all_tree.at(ibin)<0. || sum_normwgts_nonzerorewgt_tree.at(ibin)<0.){
        MELAout << "\t\t- Bin " << ibin << " has a negative sum of norm. weights before or after threshold applications. The contribution of this sample from this bin will be set to 0." << endl;
        sum_normwgts_all_tree.at(ibin) = 0;
        sum_normwgts_nonzerorewgt_tree.at(ibin) = 0;
        for (auto& sum_wgts_withrewgt_tree_hypo:sum_wgts_withrewgt_tree){ sum_wgts_withrewgt_tree_hypo.at(ibin).first = sum_wgts_withrewgt_tree_hypo.at(ibin).second = 0; }
        sampleZeroMECompensation_tree.at(ibin) = 0;
      }
      else if (sum_normwgts_nonzerorewgt_tree.at(ibin)>0.) sampleZeroMECompensation_tree.at(ibin) = sum_normwgts_all_tree.at(ibin) / sum_normwgts_nonzerorewgt_tree.at(ibin);
      if (ihypo_Neff>=0){
        auto& sum_wgts_pair = sum_wgts_withrewgt_tree.at(ihypo_Neff).at(ibin);
        if (sum_wgts_pair.second==0.) continue;
        NeffsPerBin_tree.at(ibin) = std::pow(sum_wgts_pair.first, 2)/sum_wgts_pair.second;
      }
      else{
        double Neff_bin_worst = -1;
        for (auto const& sum_wgts_withrewgt_tree_hypo:sum_wgts_withrewgt_tree){
          auto& sum_wgts_pair = sum_wgts_withrewgt_tree_hypo.at(ibin);
          if (sum_wgts_pair.second==0.) continue;
          double Neff_hypo_bin = std::pow(sum_wgts_pair.first, 2)/sum_wgts_pair.second;
          if (Neff_bin_worst<0. || Neff_bin_worst>Neff_hypo_bin) Neff_bin_worst = Neff_hypo_bin;
        }
        NeffsPerBin_tree.at(ibin) = std::max(0., Neff_bin_worst);
      }
    }

    // Assign the vectors to the maps
    NeffsPerBin[tree] = NeffsPerBin_tree;
    sum_normwgts_all[tree] = sum_normwgts_all_tree;
    sum_normwgts_nonzerorewgt[tree] = sum_normwgts_nonzerorewgt_tree;
    sum_wgts_withrewgt[tree] = sum_wgts_withrewgt_tree;
    sampleZeroMECompensation[tree] = sampleZeroMECompensation_tree;

    tree->recoverMutedBranches();
  }

  // Once all trees are complete, compute final normalization factors based on the chosen values of Neff.

  // First pass on relative normalizations
  if (thr_frac_Neff>0.){
    MELAout << "\t- Adjusting Neff values for Neff threshold " << thr_frac_Neff << ":" << endl;
    std::unordered_map<BaseTree*, unsigned int> tree_bestBin_map;
    bool useBestBins = (strBinningVars.size()==1);
    if (useBestBins){
      for (auto const& tree:registeredTrees){
        auto const& NeffsPerBin_tree = NeffsPerBin.find(tree)->second;
        double best_Neff_val = -1;
        int best_bin_index = -1;
        for (unsigned int ibin=0; ibin<nbins; ibin++){
          double const& NeffsPerBin_bin = NeffsPerBin_tree.at(ibin);
          if (NeffsPerBin_bin==0.) continue;
          if (best_Neff_val<0. || best_Neff_val<NeffsPerBin_bin){
            best_Neff_val = NeffsPerBin_bin;
            best_bin_index = ibin;
          }
        }
        if (best_bin_index>=0) tree_bestBin_map[tree] = best_bin_index;
      }
    }

    for (unsigned int ibin=0; ibin<nbins; ibin++){
      double Neff_total = 0;
      for (auto const& tree:registeredTrees){
        auto const& NeffsPerBin_tree = NeffsPerBin.find(tree)->second;
        double const& NeffsPerBin_bin = NeffsPerBin_tree.at(ibin);
        Neff_total += NeffsPerBin_bin;
      }
      for (auto const& tree:registeredTrees){
        auto& NeffsPerBin_tree = NeffsPerBin.find(tree)->second;
        auto& sum_normwgts_all_tree = sum_normwgts_all.find(tree)->second;
        auto& sum_normwgts_nonzerorewgt_tree = sum_normwgts_nonzerorewgt.find(tree)->second;
        auto& sum_wgts_withrewgt_tree = sum_wgts_withrewgt.find(tree)->second;
        auto& sampleZeroMECompensation_tree = sampleZeroMECompensation.find(tree)->second;
        int tree_bin_best = -1;
        auto it_tree_bestBin = tree_bestBin_map.find(tree);
        if (it_tree_bestBin!=tree_bestBin_map.cend()) tree_bin_best = it_tree_bestBin->second;

        double& NeffsPerBin_bin = NeffsPerBin_tree.at(ibin);
        double const frac_Neff = NeffsPerBin_bin/Neff_total;
        if (frac_Neff<thr_frac_Neff && (tree_bin_best<0 || std::abs(tree_bin_best-static_cast<int>(ibin))>1)){
          MELAout << "\t\t- Bin " << ibin << " in sample " << tree->sampleIdentifier << " has a fractional Neff=" << frac_Neff << " < " << thr_frac_Neff << ". Removing this bin from reweighting considerations." << endl;
          NeffsPerBin_bin = 0;
          sum_normwgts_all_tree.at(ibin) = 0;
          sum_normwgts_nonzerorewgt_tree.at(ibin) = 0;
          for (auto& sum_wgts_withrewgt_tree_hypo:sum_wgts_withrewgt_tree){ sum_wgts_withrewgt_tree_hypo.at(ibin).first = sum_wgts_withrewgt_tree_hypo.at(ibin).second = 0; }
          sampleZeroMECompensation_tree.at(ibin) = 0;
        }
      }
    }

    // If best bins are used, ensure that the sample sets bins to 0 beyond the first encounter of 0.
    if (useBestBins){
      for (auto const& tree:registeredTrees){
        auto& NeffsPerBin_tree = NeffsPerBin.find(tree)->second;
        auto& sum_normwgts_all_tree = sum_normwgts_all.find(tree)->second;
        auto& sum_normwgts_nonzerorewgt_tree = sum_normwgts_nonzerorewgt.find(tree)->second;
        auto& sum_wgts_withrewgt_tree = sum_wgts_withrewgt.find(tree)->second;
        auto& sampleZeroMECompensation_tree = sampleZeroMECompensation.find(tree)->second;
        int tree_bin_best = -1;
        auto it_tree_bestBin = tree_bestBin_map.find(tree);
        if (it_tree_bestBin!=tree_bestBin_map.cend()) tree_bin_best = it_tree_bestBin->second;

        // Move to the right
        {
          bool seeZero = false;
          for (unsigned int ibin = tree_bin_best; ibin<nbins; ibin++){
            double& NeffsPerBin_bin = NeffsPerBin_tree.at(ibin);
            if (NeffsPerBin_bin==0.) seeZero = true;
            if (seeZero){
              MELAout << "\t\t- Bin " << ibin << " in sample " << tree->sampleIdentifier << " has acceptable Neff, but it is beyond the first occurence of 0s from the right. Removing this bin from reweighting considerations." << endl;
              NeffsPerBin_bin = 0;
              sum_normwgts_all_tree.at(ibin) = 0;
              sum_normwgts_nonzerorewgt_tree.at(ibin) = 0;
              for (auto& sum_wgts_withrewgt_tree_hypo:sum_wgts_withrewgt_tree){ sum_wgts_withrewgt_tree_hypo.at(ibin).first = sum_wgts_withrewgt_tree_hypo.at(ibin).second = 0; }
              sampleZeroMECompensation_tree.at(ibin) = 0;
            }
          }
        }
        // Move to the left
        {
          bool seeZero = false;
          for (int ibin = tree_bin_best; ibin>=0; ibin--){
            double& NeffsPerBin_bin = NeffsPerBin_tree.at(ibin);
            if (NeffsPerBin_bin==0.) seeZero = true;
            if (seeZero){
              MELAout << "\t\t- Bin " << ibin << " in sample " << tree->sampleIdentifier << " has acceptable Neff, but it is beyond the first occurence of 0s from the left. Removing this bin from reweighting considerations." << endl;
              NeffsPerBin_bin = 0;
              sum_normwgts_all_tree.at(ibin) = 0;
              sum_normwgts_nonzerorewgt_tree.at(ibin) = 0;
              for (auto& sum_wgts_withrewgt_tree_hypo:sum_wgts_withrewgt_tree){ sum_wgts_withrewgt_tree_hypo.at(ibin).first = sum_wgts_withrewgt_tree_hypo.at(ibin).second = 0; }
              sampleZeroMECompensation_tree.at(ibin) = 0;
            }
          }
        }
      }
    }
  }

  // Second pass on relative normalizations
  MELAout << "\t- Sample normalizations before extra normalizations:" << endl;
  for (auto const& tree:registeredTrees){
    auto& sampleNormalization_tree = sampleNormalization.find(tree)->second;
    auto const& NeffsPerBin_tree = NeffsPerBin.find(tree)->second;
    for (unsigned int ibin=0; ibin<nbins; ibin++){
      double& sampleNormalization_bin = sampleNormalization_tree.at(ibin);
      double const& NeffsPerBin_bin = NeffsPerBin_tree.at(ibin);
      double sum_Neff_bin = 0;
      for (auto const& it:NeffsPerBin) sum_Neff_bin += it.second.at(ibin);
      if (sum_Neff_bin>0.) sampleNormalization_bin = NeffsPerBin_bin / sum_Neff_bin;
    }
    MELAout << "[" << tree->sampleIdentifier << "]: " << sampleNormalization_tree << endl;
  }

  // Pairwise normalizations
  if (tree_normTree_pairs && ihypo_Neff>=0){
    MELAout << "\t- Computing the pairwise normalization factors..." << endl;
    std::unordered_map<BaseTree*, double> extraNormFactors;
    // Determine the pairwise normalizations first
    for (auto const& tree_normTree_pair:(*tree_normTree_pairs)){
      auto const& tree = tree_normTree_pair.first;
      auto const& normTree = tree_normTree_pair.second;
      auto const& sum_wgts_withrewgt_basetree = sum_wgts_withrewgt[tree].at(ihypo_Neff);
      auto const& sum_wgts_withrewgt_normtree = sum_wgts_withrewgt[normTree].at(ihypo_Neff);

      // 0-MEs compensation factors need to be accounted for in the 
      auto const& sampleZeroMECompensation_basetree = sampleZeroMECompensation.find(tree)->second;
      auto const& sampleZeroMECompensation_normtree = sampleZeroMECompensation.find(normTree)->second;

      double sum_wgts_common_basetree = 0;
      double sum_wgts_common_normtree = 0;
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        double const& sum_wgts_common_basetree_bin = sum_wgts_withrewgt_basetree.at(ibin).first;
        double const& sum_wgts_common_normtree_bin = sum_wgts_withrewgt_normtree.at(ibin).first;
        if (sum_wgts_common_basetree_bin!=0. && sum_wgts_common_normtree_bin!=0.){
          double const& nonzeroNorm_basetree = sampleZeroMECompensation_basetree.at(ibin);
          double const& nonzeroNorm_normtree = sampleZeroMECompensation_normtree.at(ibin);
          sum_wgts_common_basetree += sum_wgts_common_basetree_bin * nonzeroNorm_basetree;
          sum_wgts_common_normtree += sum_wgts_common_normtree_bin * nonzeroNorm_normtree;
        }
      }
      if (sum_wgts_common_basetree==0.) MELAerr << "\t\t- Base tree " << tree->sampleIdentifier << " has no overlap with its norm tree " << normTree->sampleIdentifier << endl;
      if (sum_wgts_common_normtree==0.) MELAerr << "\t\t- Norm tree " << tree->sampleIdentifier << " has no overlap with its base tree " << normTree->sampleIdentifier << endl;

      double extra_norm = 1;
      if (sum_wgts_common_basetree!=0.) extra_norm = sum_wgts_common_normtree / sum_wgts_common_basetree;
      MELAout << "\t\t- Base tree " << tree->sampleIdentifier << " has an extra normalization of " << extra_norm << endl;
      extraNormFactors[tree] = extra_norm;
    }
    // Multiply all relevant ones
    MELAout << "\t- Applying the computed pairwise normalization factors..." << endl;
    for (auto const& tree_normTree_pair:(*tree_normTree_pairs)){
      BaseTree* tree = tree_normTree_pair.first;
      BaseTree* baseTree = tree;
      BaseTree* normTree = tree_normTree_pair.second;
      while (normTree!=nullptr){
        MELAout << "\t\t- Normalizing " << tree->sampleIdentifier << " by the normalization factor (=" << extraNormFactors[baseTree] << ") of " << baseTree->sampleIdentifier << endl;
        for (auto& v:sampleNormalization[tree]) v *= extraNormFactors[baseTree];
        samplePairwiseNormalization[tree] *= extraNormFactors[baseTree];
        baseTree = normTree; normTree = nullptr;
        for (auto const& pp:(*tree_normTree_pairs)){
          if (pp.first == baseTree){
            normTree = pp.second;
            break;
          }
        }
      }
    }
    MELAout << "\t- Sample normalizations after pairwise normalizations:" << endl;
    for (auto const& tree:registeredTrees) MELAout << "[" << tree->sampleIdentifier << "]: " << sampleNormalization[tree] << endl;
  }

  // Finally, add the 0-ME compensation factors
  MELAout << "\t- Sample normalizations after 0-ME compensation factors:" << endl;
  for (auto const& tree:registeredTrees){
    auto& sampleNormalization_tree = sampleNormalization.find(tree)->second;
    auto const& sampleZeroMECompensation_tree = sampleZeroMECompensation.find(tree)->second;
    for (unsigned int ibin=0; ibin<nbins; ibin++) sampleNormalization_tree.at(ibin) *= sampleZeroMECompensation_tree.at(ibin);
    MELAout << "[" << tree->sampleIdentifier << "]: " << sampleNormalization_tree << endl;
  }

  MELAout << "BulkReweightingBuilder::setup has completed successfully." << endl;
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
double BulkReweightingBuilder::getSamplePairwiseNormalization(BaseTree* tree) const{
  auto it = samplePairwiseNormalization.find(tree);
  if (it!=samplePairwiseNormalization.cend()) return it->second;
  else{
    MELAerr << "BulkReweightingBuilder::getSamplePairwiseNormalization: No normalization factor is found for tree " << tree->sampleIdentifier << "." << endl;
    return 1;
  }
}
bool BulkReweightingBuilder::checkWeightsBelowThreshold(BaseTree* tree) const{
  auto it_binningVarRefs = binningVarRefs.find(tree);
  if (it_binningVarRefs == binningVarRefs.cend()){
    MELAerr << "BulkReweightingBuilder::checkWeightsBelowThreshold: Tree " << tree->sampleIdentifier << " is not registered properly." << endl;
    return false;
  }

  unsigned int const nbins = (binning.isValid() ? binning.getNbins() : static_cast<unsigned int>(1));
  int ibin = rule_binningVar(tree, binning, it_binningVarRefs->second);
  if (ibin<0) ibin=0;
  else if (ibin>=(int) nbins) ibin = nbins-1;

  bool res = true;
  auto const& absWeightThresholdsPerBin = absWeightThresholdsPerBinList.find(tree)->second;
  auto const& componentRefs_reweightingweights = componentRefsList_reweightingweights.find(tree)->second;
  unsigned int const nhypos = strReweightingWeightsList.size();
  for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
    float const& wgt_thr = absWeightThresholdsPerBin.at(ihypo).at(ibin);
    float wgt_rewgt = rule_reweightingweights_list.at(ihypo)(tree, componentRefs_reweightingweights.at(ihypo));
    if (wgt_thr>0.f && std::abs(wgt_rewgt)>wgt_thr){
      res = false;
      break;
    }
  }

  return res;
}

void BulkReweightingBuilder::writeToFile(TFile* foutput) const{
  unsigned int const nhypos = strReweightingWeightsList.size();
  unsigned int const nbins = (binning.isValid() ? binning.getNbins() : static_cast<unsigned int>(1));

  TDirectory* curdir = gDirectory;
  foutput->cd();

  for (auto const& tree:registeredTrees){
    TString const sid = tree->sampleIdentifier;
    TDirectory* dir_tree = foutput->mkdir(sid);
    dir_tree->cd();

    TH1D* h_normWeights = new TH1D("normWeight", sid, 1, 0, 1);
    h_normWeights->SetBinContent(1, normWeights.find(tree)->second);
    dir_tree->WriteTObject(h_normWeights); delete h_normWeights;

    TH1D* h_samplePairwiseNormalization = new TH1D("samplePairwiseNormalization", sid, 1, 0, 1);
    h_samplePairwiseNormalization->SetBinContent(1, samplePairwiseNormalization.find(tree)->second);
    dir_tree->WriteTObject(h_samplePairwiseNormalization); delete h_samplePairwiseNormalization;

#define HIST_1D_COMMANDS \
HIST_1D_COMMAND(sum_normwgts_all) \
HIST_1D_COMMAND(sum_normwgts_nonzerorewgt) \
HIST_1D_COMMAND(NeffsPerBin) \
HIST_1D_COMMAND(sampleNormalization) \
HIST_1D_COMMAND(sampleZeroMECompensation)
#define HIST_1D_COMMAND(NAME) \
    auto const& NAME##_tree = NAME.find(tree)->second; \
    TH1D* h_##NAME = new TH1D(#NAME, sid, nbins, 0, nbins); \
    for (unsigned int ibin=0; ibin<nbins; ibin++) h_##NAME->SetBinContent(ibin+1, NAME##_tree.at(ibin)); \
    dir_tree->WriteTObject(h_##NAME); delete h_##NAME;
    HIST_1D_COMMANDS;
#undef HIST_1D_COMMAND
#undef HIST_1D_COMMANDS

    auto const& absWeightThresholdsPerBinList_tree = absWeightThresholdsPerBinList.find(tree)->second;
    TH2D* h_absWeightThresholdsPerBinList = new TH2D("absWeightThresholdsPerBinList", sid, nhypos, 0, nhypos, nbins, 0, nbins);
    for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        h_absWeightThresholdsPerBinList->SetBinContent(ihypo+1, ibin+1, absWeightThresholdsPerBinList_tree.at(ihypo).at(ibin));
      }
    }
    dir_tree->WriteTObject(h_absWeightThresholdsPerBinList); delete h_absWeightThresholdsPerBinList;

    auto const& sum_wgts_withrewgt_tree = sum_wgts_withrewgt.find(tree)->second;
    TH2D* h_sum_wgts_withrewgt = new TH2D("sum_wgts_withrewgt", sid, nhypos, 0, nhypos, nbins, 0, nbins);
    for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        h_sum_wgts_withrewgt->SetBinContent(ihypo+1, ibin+1, sum_wgts_withrewgt_tree.at(ihypo).at(ibin).first);
        h_sum_wgts_withrewgt->SetBinError(ihypo+1, ibin+1, sum_wgts_withrewgt_tree.at(ihypo).at(ibin).second);
      }
    }
    dir_tree->WriteTObject(h_sum_wgts_withrewgt); delete h_sum_wgts_withrewgt;

    dir_tree->Close();
    foutput->cd();
  }

  // Common stuff
  TH2D* h_reweightingweights_frac_tolerance_pair_list = new TH2D("reweightingweights_frac_tolerance_pair_list", "", nhypos, 0, nhypos, 2, 0, 2);
  for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
    h_reweightingweights_frac_tolerance_pair_list->SetBinContent(ihypo+1, 1, reweightingweights_frac_tolerance_pair_list.at(ihypo).first);
    h_reweightingweights_frac_tolerance_pair_list->SetBinContent(ihypo+1, 2, reweightingweights_frac_tolerance_pair_list.at(ihypo).second);
  }
  foutput->WriteTObject(h_reweightingweights_frac_tolerance_pair_list); delete h_reweightingweights_frac_tolerance_pair_list;

  TTree* t_common = new TTree("CommonVariables", "");
  TString binning_name = binning.getName(); t_common->Branch("binning_name", &binning_name);
  std::vector<double> bin_thresholds = binning.getBinningVector(); t_common->Branch("bin_thresholds", &bin_thresholds);
  auto v_strBinningVars = strBinningVars; t_common->Branch("strBinningVars", &v_strBinningVars);
  auto v_strNominalWeights = strNominalWeights; t_common->Branch("strNominalWeights", &v_strNominalWeights);
  auto v_strCrossSectionWeights = strCrossSectionWeights; t_common->Branch("strCrossSectionWeights", &v_strCrossSectionWeights);
  auto v_strReweightingWeightsList = strReweightingWeightsList; t_common->Branch("strReweightingWeightsList", &v_strReweightingWeightsList);
  t_common->Fill();
  foutput->WriteTObject(t_common); delete t_common;

  curdir->cd();
}
void BulkReweightingBuilder::setupFromFile(TString cinput){
  HostHelpers::ExpandEnvironmentVariables(cinput);

  if (!HostHelpers::FileReadable(cinput.Data())){
    MELAerr << "BulkReweightingBuilder::setupFromFile: File " << cinput << " is not readable." << endl;
    assert(0);
  }
  MELAout << "BulkReweightingBuilder::setupFromFile: setting up reweighting from file " << cinput << "." << endl;
  TDirectory* curdir = gDirectory;
  TFile* finput = TFile::Open(cinput, "read");
  finput->cd();

  // Cross-check the requested weighting scheme
  TTree* t_common = (TTree*) finput->Get("CommonVariables");
  TString binning_name; t_common->SetBranchAddress("binning_name", &binning_name);
  std::vector<double>* bin_thresholds = nullptr; t_common->SetBranchAddress("bin_thresholds", &bin_thresholds);
  std::vector<TString>* v_strBinningVars = nullptr; t_common->SetBranchAddress("strBinningVars", &v_strBinningVars);
  std::vector<TString>* v_strNominalWeights = nullptr; t_common->SetBranchAddress("strNominalWeights", &v_strNominalWeights);
  std::vector<TString>* v_strCrossSectionWeights = nullptr; t_common->SetBranchAddress("strCrossSectionWeights", &v_strCrossSectionWeights);
  std::vector<std::vector<TString>>* v_strReweightingWeightsList = nullptr; t_common->SetBranchAddress("strReweightingWeightsList", &v_strReweightingWeightsList);
  t_common->GetEntry(0);

  // Check binning variables
  {
    bool hasSameBinning = (v_strBinningVars->size() == strBinningVars.size());
    if (hasSameBinning){
      for (auto const& v:strBinningVars){
        if (!HelperFunctions::checkListVariable(*v_strBinningVars, v)){
          hasSameBinning = false;
          break;
        }
      }
    }
    if (!hasSameBinning){
      MELAerr << "\t- Binning variables were specified as " << strBinningVars << ", but the input file has " << *v_strBinningVars << "." << endl;
      assert(0);
    }
    // Reset binning
    binning = ExtendedBinning(*bin_thresholds, binning_name, "");
  }

  // Check nominal weight variables
  {
    bool hasSameWeights = (v_strNominalWeights->size() == strNominalWeights.size());
    if (hasSameWeights){
      for (auto const& v:strNominalWeights){
        if (!HelperFunctions::checkListVariable(*v_strNominalWeights, v)){
          hasSameWeights = false;
          break;
        }
      }
    }
    if (!hasSameWeights){
      MELAerr << "\t- Nominal weight variables were specified as " << strNominalWeights << ", but the input file has " << *v_strNominalWeights << "." << endl;
      assert(0);
    }
  }

  // Check nominal weight variables
  {
    bool hasSameWeights = (v_strCrossSectionWeights->size() == strCrossSectionWeights.size());
    if (hasSameWeights){
      for (auto const& v:strCrossSectionWeights){
        if (!HelperFunctions::checkListVariable(*v_strCrossSectionWeights, v)){
          hasSameWeights = false;
          break;
        }
      }
    }
    if (!hasSameWeights){
      MELAerr << "\t- Cross section weight variables were specified as " << strCrossSectionWeights << ", but the input file has " << *v_strCrossSectionWeights << "." << endl;
      assert(0);
    }
  }

  unsigned int const nhypos = strReweightingWeightsList.size();
  unsigned int const nbins = (binning.isValid() ? binning.getNbins() : static_cast<unsigned int>(1));

  std::vector<unsigned int> hypo_order;
  for (auto const& strReweightingWeights:strReweightingWeightsList){
    int jhypo = -1;
    for (unsigned int ihypo=0; ihypo<v_strReweightingWeightsList->size(); ihypo++){
      auto const& v_strReweightingWeights = v_strReweightingWeightsList->at(ihypo);
      bool hasSameWeights = (v_strReweightingWeights.size() == strReweightingWeights.size());
      if (hasSameWeights){
        for (auto const& v:strReweightingWeights){
          if (!HelperFunctions::checkListVariable(v_strReweightingWeights, v)){
            hasSameWeights = false;
            break;
          }
        }
      }
      if (hasSameWeights){ jhypo = ihypo; break; }
    }
    if (jhypo<0){
      MELAout << "\t- Reweighting hypothesis with weights " << strReweightingWeights << " is not specified in the input file." << endl;
      assert(0);
    }
    hypo_order.push_back(jhypo);
  }
  assert(nhypos==hypo_order.size());

  // Set reweightingweights_frac_tolerance_pair_list
  reweightingweights_frac_tolerance_pair_list.assign(nhypos, std::pair<double, double>(0, 0));
  TH2D* h_reweightingweights_frac_tolerance_pair_list = (TH2D*) finput->Get("reweightingweights_frac_tolerance_pair_list");
  for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
    unsigned int const& jhypo = hypo_order.at(ihypo);
    reweightingweights_frac_tolerance_pair_list.at(ihypo).first = h_reweightingweights_frac_tolerance_pair_list->GetBinContent(jhypo+1, 1);
    reweightingweights_frac_tolerance_pair_list.at(ihypo).second = h_reweightingweights_frac_tolerance_pair_list->GetBinContent(jhypo+1, 2);
  }

  // Set unordered maps for the trees
  for (auto const& tree:registeredTrees){
    TString const sid = tree->sampleIdentifier;
    TDirectory* dir_tree = (TDirectory*) finput->Get(sid);
    if (!dir_tree){
      MELAout << "\t- No directory named " << sid << endl;
      assert(0);
    }
    dir_tree->cd();

    TH1D* h_normWeights = (TH1D*) dir_tree->Get("normWeight");
    normWeights[tree] = h_normWeights->GetBinContent(1);

    TH1D* h_samplePairwiseNormalization = (TH1D*) dir_tree->Get("samplePairwiseNormalization");
    samplePairwiseNormalization[tree] = h_samplePairwiseNormalization->GetBinContent(1);

#define HIST_1D_COMMANDS \
HIST_1D_COMMAND(sum_normwgts_all) \
HIST_1D_COMMAND(sum_normwgts_nonzerorewgt) \
HIST_1D_COMMAND(NeffsPerBin) \
HIST_1D_COMMAND(sampleNormalization) \
HIST_1D_COMMAND(sampleZeroMECompensation)
#define HIST_1D_COMMAND(NAME) \
    auto const& NAME##_tree = NAME.find(tree)->second; \
    NAME[tree] = std::vector<double>(nbins, 0); \
    TH1D* h_##NAME = (TH1D*) dir_tree->Get(#NAME); \
    for (unsigned int ibin=0; ibin<nbins; ibin++) NAME[tree].at(ibin) = h_##NAME->GetBinContent(ibin+1); \
    HIST_1D_COMMANDS;
#undef HIST_1D_COMMAND
#undef HIST_1D_COMMANDS

    absWeightThresholdsPerBinList[tree] = std::vector< std::vector<float> >(nhypos, std::vector<float>(nbins, -1));
    TH2D* h_absWeightThresholdsPerBinList = (TH2D*) dir_tree->Get("absWeightThresholdsPerBinList");
    for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
      unsigned int const& jhypo = hypo_order.at(ihypo);
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        absWeightThresholdsPerBinList[tree].at(ihypo).at(ibin) = h_absWeightThresholdsPerBinList->GetBinContent(jhypo+1, ibin+1);
      }
    }

    sum_wgts_withrewgt[tree] = std::vector< std::vector< std::pair<double, double> > >(nhypos, std::vector< std::pair<double, double> >(nbins, std::pair<double, double>(0, 0)));
    TH2D* h_sum_wgts_withrewgt = (TH2D*) dir_tree->Get("sum_wgts_withrewgt");
    for (unsigned int ihypo=0; ihypo<nhypos; ihypo++){
      unsigned int const& jhypo = hypo_order.at(ihypo);
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        sum_wgts_withrewgt[tree].at(ihypo).at(ibin).first = h_sum_wgts_withrewgt->GetBinContent(jhypo+1, ibin+1);
        sum_wgts_withrewgt[tree].at(ihypo).at(ibin).second = h_sum_wgts_withrewgt->GetBinError(jhypo+1, ibin+1);
      }
    }

    finput->cd();
  }

  finput->Close();
  curdir->cd();
}
