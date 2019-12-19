#include "ReweightingBuilder.h"
#include "SimpleEntry.h"
#include "MELAAccumulators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


ReweightingBuilder::ReweightingBuilder(
  std::vector<TString> const& strNominalWeights_,
  std::vector<TString> const& strReweightingWeights_,
  std::vector<TString> const& strCrossSectionWeights_,
  ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_,
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
) :
  strNominalWeights(strNominalWeights_),
  strReweightingWeights(strReweightingWeights_),
  strCrossSectionWeights(strCrossSectionWeights_),
  rule_nominalweights(rule_nominalweights_),
  rule_reweightingweights(rule_reweightingweights_),
  rule_xsecweights(rule_xsecweights_)
{}


float ReweightingBuilder::eval_nominalweights(BaseTree* theTree) const{
  std::unordered_map<BaseTree*, std::vector<float*>>::const_iterator it = componentRefs_nominalweights.find(theTree);
  if (it==componentRefs_nominalweights.cend()){
    MELAerr << "ReweightingBuilder::eval_nominalweights: Could not find the weights " << strNominalWeights << ". Call ReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }

  float weight=0;
  if (rule_nominalweights) weight=rule_nominalweights(theTree, it->second);
  if (!checkVarNanInf(weight)){
    MELAerr << "ReweightingBuilder::eval_nominalweights: Weight " << weight << " is nan/inf!" << endl;
    weight=0;
  }
  if (!allowNegativeWeights && weight<0.){
    weight=0;
    MELAerr << "ReweightingBuilder::eval_nominalweights: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
  return weight;
}
float ReweightingBuilder::eval_reweightingweights(BaseTree* theTree) const{
  std::unordered_map<BaseTree*, std::vector<float*>>::const_iterator it = componentRefs_reweightingweights.find(theTree);
  if (it==componentRefs_reweightingweights.cend()){
    MELAerr << "ReweightingBuilder::eval_reweightingweights: Could not find the weights " << strReweightingWeights << ". Call ReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }

  float weight=0;
  if (rule_reweightingweights) weight=rule_reweightingweights(theTree, it->second);
  if (!checkVarNanInf(weight)){
    MELAerr << "ReweightingBuilder::eval_reweightingweights: Weight " << weight << " is nan/inf!" << endl;
    weight=0;
  }
  if (!allowNegativeWeights && weight<0.){
    weight=0;
    MELAerr << "ReweightingBuilder::eval_reweightingweights: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
  return weight;
}
float ReweightingBuilder::eval_xsecweights(BaseTree* theTree) const{
  std::unordered_map<BaseTree*, std::vector<float*>>::const_iterator it = componentRefs_xsecweights.find(theTree);
  if (it==componentRefs_xsecweights.cend()){
    MELAerr << "ReweightingBuilder::eval_xsecweights: Could not find the weights " << strCrossSectionWeights << ". Call ReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }

  float weight=0;
  if (rule_xsecweights) weight=rule_xsecweights(theTree, it->second);
  if (!checkVarNanInf(weight)){
    MELAerr << "ReweightingBuilder::eval_xsecweights: Weight " << weight << " is nan/inf!" << endl;
    weight=0;
  }
  if (!allowNegativeWeights && weight<0.){
    weight=0;
    MELAerr << "ReweightingBuilder::eval_xsecweights: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
  return weight;
}

int ReweightingBuilder::findBin(BaseTree* theTree) const{
  const ExtendedBinning& binning = weightBinning;
  int bin=0;
  if (binning.isValid()){
    auto binningVarRefsIt = binningVarRefs.find(theTree);
    if (binningVarRefsIt!=binningVarRefs.cend()){
      float const& orderingVal=*(binningVarRefsIt->second);
      bin = binning.getBin(orderingVal);
      if (bin<0 || bin>=(int) binning.getNbins()) bin=-1;
    }
  }
  return bin;
}
void ReweightingBuilder::rejectNegativeWeights(const bool flag){ allowNegativeWeights = !flag; }
void ReweightingBuilder::setWeightBinning(const ExtendedBinning& binning){ weightBinning = binning; }
void ReweightingBuilder::setWeightThresholdReference(const float& weightThresholdReference_){ weightThresholdReference = weightThresholdReference_; }

void ReweightingBuilder::setupWeightVariables(BaseTree* theTree, float fractionRequirement, unsigned int minimumNevents){
  MELAout << "ReweightingBuilder[" << "Nominal weights: " << strNominalWeights << ", reweighting weights: " << strReweightingWeights << ", xsec weights: " << strCrossSectionWeights << "]::setupWeightVariables is called for tree " << theTree << "." << endl;

  if (!theTree) return;
  const int nevents = theTree->getNEvents();

  const ExtendedBinning& binning = this->weightBinning;
  TString strOrderingVal = binning.getLabel();
  if (strOrderingVal=="") return;
  float* orderingValRef = ReweightingFunctions::getWeightRef(theTree, strOrderingVal);
  if (!orderingValRef) return;

  // First link components
  componentRefs_nominalweights[theTree] = ReweightingFunctions::getWeightRefs(theTree, strNominalWeights);
  componentRefs_reweightingweights[theTree] = ReweightingFunctions::getWeightRefs(theTree, strReweightingWeights);
  componentRefs_xsecweights[theTree] = ReweightingFunctions::getWeightRefs(theTree, strCrossSectionWeights);
  binningVarRefs[theTree] = orderingValRef;

  // Now compute other weight-related variables
  const bool noBoundaries = !binning.isValid();
  const unsigned int ns = (!noBoundaries ? binning.getNbins() : 1);

  // Initialize
  bool doReweighting = (fractionRequirement>=0.);
  sumPostThrWeights[theTree]=std::vector<float>();
  sumPostThrSqWeights[theTree]=std::vector<float>();
  weightThresholds[theTree]=std::vector<float>();
  sumNonZeroWgtEvents[theTree]=std::vector<unsigned int>();
  sumNonZeroWgtNominalWeights[theTree]=std::vector<float>();
  sumPostThrWeights[theTree].assign(ns, 0);
  sumPostThrSqWeights[theTree].assign(ns, 0);
  sumNominalWeights[theTree] = 0;
  xsecVals[theTree]=1;
  if (doReweighting){
    weightThresholds[theTree].assign(ns, 0);
    sumNonZeroWgtEvents[theTree].assign(ns, 0);
    sumNonZeroWgtNominalWeights[theTree].assign(ns, 0);
  }
  else{
    weightThresholds[theTree].assign(ns, -1);
    sumNonZeroWgtEvents[theTree].assign(ns, -1);
    sumNonZeroWgtNominalWeights[theTree].assign(ns, 0);
  }

  vector<vector<SimpleEntry>> indexList;
  indexList.assign(ns, vector<SimpleEntry>());
  for (auto& index:indexList) index.reserve(theTree->getNEvents()/ns+1);

  MELAout << "\t- Ordering the " << nevents << " events";
  if (!noBoundaries) MELAout << " over the " << ns << " bins: [ " << binning.getBinningVector() << " ]";
  MELAout << " and computing sum of nominal weights without xsec." << endl;
  {
    bool firstTreeEvent=true;
    auto it_sumNominalWeights = sumNominalWeights.find(theTree);
    for (int ev=0; ev<theTree->getNEvents(); ev++){
      HelperFunctions::progressbar(ev, nevents);

      bool hasEvent = theTree->getEvent(ev);
      if (!hasEvent) continue;

      if (firstTreeEvent) xsecVals[theTree] = this->eval_xsecweights(theTree);

      float wgt_nominal = this->eval_nominalweights(theTree);
      it_sumNominalWeights->second += wgt_nominal;

      if (doReweighting){
        float wgt = wgt_nominal * this->eval_reweightingweights(theTree);

        float const& orderingVal = *orderingValRef;
        SimpleEntry theEntry(0, fabs(wgt-weightThresholdReference), wgt);

        int bin=0;
        if (!noBoundaries) bin = binning.getBin(orderingVal);
        if (bin>=0 && bin<(int) ns){
          // Accumulate the events
          indexList.at(bin).push_back(theEntry);
          if (theEntry.trackingval!=0.){
            sumNonZeroWgtEvents[theTree].at(bin)++;
            sumNonZeroWgtNominalWeights[theTree].at(bin) += wgt_nominal;
          }
        }
      }

      if (firstTreeEvent) firstTreeEvent=false;
    }
    MELAout << "\t- Sum of nominal weights is computed as " << it_sumNominalWeights->second << "." << endl;
  }

  if (!doReweighting) return;

  for (unsigned int ibin=0; ibin<ns; ibin++){
    vector<SimpleEntry>& index = indexList.at(ibin);
    if (minimumNevents>sumNonZeroWgtEvents[theTree].at(ibin)){
      MELAout << "\t- Bin " << ibin << " has less number of events with non-zero weight than the requested number " << minimumNevents << ". Resetting the bin..." << endl;
      index.clear();

      sumNonZeroWgtEvents[theTree].at(ibin)=0;
      sumNonZeroWgtNominalWeights[theTree].at(ibin)=0;
    }
    const unsigned int nTotalPerBin = index.size();
    //MELAout << "\t- Looping over the " << nTotalPerBin << " events to find the threshold in bin " << ibin << endl;
    float threshold=0;
    if (nTotalPerBin>2){
      const unsigned int maxPrunedSize = std::max((unsigned int) (fractionRequirement>=0. ? std::ceil(float(nTotalPerBin)*(1.-fractionRequirement)) : nTotalPerBin), (unsigned int) 3);
      vector<SimpleEntry> indexPruned; indexPruned.reserve(maxPrunedSize);

      //unsigned int itrk=0;
      for (auto& theEntry:index){
        //HelperFunctions::progressbar(itrk, nTotalPerBin);
        if (indexPruned.size()<maxPrunedSize) addByHighest(indexPruned, theEntry, false);
        else if (indexPruned.at(indexPruned.size()-1).trackingval<theEntry.trackingval){
          addByHighest(indexPruned, theEntry, false);
          indexPruned.pop_back();
        }
        //itrk++;
      }

      // Find the threshold
      unsigned int index_entry = maxPrunedSize-2;
      unsigned int index_entry_prev=index_entry+1;
      threshold = (indexPruned.at(index_entry_prev).trackingval + indexPruned.at(index_entry).trackingval)*0.5;
      MELAout << "Threshold raw: " << threshold << endl;
      if ((threshold>0. && indexPruned.front().trackingval<threshold*5.) || fractionRequirement>=1.) threshold = indexPruned.front().trackingval; // Prevent false-positives
    }
    else if (nTotalPerBin==2) threshold = std::max(index.at(0).trackingval, index.at(1).trackingval);
    else if (nTotalPerBin==1) threshold = index.at(0).trackingval;
    weightThresholds[theTree].at(ibin)=threshold;

    //MELAout << "\t- Looping over the " << nTotalPerBin << " events to find the sum of weights after threshold " << threshold << " in bin " << ibin << endl;
    // Do a precise summation with the Kahan method
    KahanAccumulator<float> sum;
    KahanAccumulator<float> sumsq;
    for (auto& theEntry:index){
      float weight = theEntry.weight;
      float weightRel = weight-weightThresholdReference;
      if (fabs(weightRel)>threshold){
        weightRel = pow(threshold, 2)/weightRel;
        weight = weightRel + weightThresholdReference;
      }
      sum += weight;
      sumsq += pow(weight, 2);
    }
    // Assign the sum
    sumPostThrWeights[theTree].at(ibin)=sum;
    sumPostThrSqWeights[theTree].at(ibin)=sumsq;
    if (sumNonZeroWgtEvents[theTree].at(ibin)>0) MELAout << "\t- Threshold at bin " << ibin << ": " << weightThresholdReference << " +- " << threshold
      << ", sum of post-threshold weights: " << sumPostThrWeights[theTree].at(ibin) << " +- " << sqrt(sumPostThrSqWeights[theTree].at(ibin))
      << ", number of events with non-zero weight: " << sumNonZeroWgtEvents[theTree].at(ibin)
      << ", sum of nominal weights for non-zero - weight events: " << sumNonZeroWgtNominalWeights[theTree].at(ibin)
      << endl;
  }
}
void ReweightingBuilder::setupCaches(){
  const ExtendedBinning& binning = weightBinning;
  int nbins;
  if (binning.isValid()) nbins = binning.getNbins();
  else nbins=1;

  std::vector<BaseTree*> const trees = this->getRegisteredTrees();

  // Clear cached objects
  cachedNormComponentsPerBin.clear();
  cachedSampleNormalizationsPerBin.clear();

  // Setup default values
  cachedNormComponentsPerBin.assign(nbins, 0);
  for (auto& tree:trees){
    cachedSampleNormalizationsPerBin[tree] = std::vector<float>(nbins, 0);
  }

  for (int bin=0; bin<nbins; bin++){
    KahanAccumulator<float> numerator;
    KahanAccumulator<float> denominator;

    float& cachedNormComponentPerBin = cachedNormComponentsPerBin.at(bin);

    for (auto& tree:trees){
      auto itSumWeights = this->sumPostThrWeights.find(tree); if (itSumWeights==this->sumPostThrWeights.cend()) continue;
      auto itSumSqWeights = this->sumPostThrSqWeights.find(tree); if (itSumSqWeights==this->sumPostThrSqWeights.cend()) continue;
      auto itSumNonZeroWgtNominalWeights = this->sumNonZeroWgtNominalWeights.find(tree); if (itSumNonZeroWgtNominalWeights==this->sumNonZeroWgtNominalWeights.cend()) continue;
      auto itSumNonZeroWgtEvents = this->sumNonZeroWgtEvents.find(tree); if (itSumNonZeroWgtEvents==this->sumNonZeroWgtEvents.cend()) continue;
      auto itSumNominalWeights = this->sumNominalWeights.find(tree); if (itSumNominalWeights==this->sumNominalWeights.cend()) continue;
      auto itXsec = this->xsecVals.find(tree); if (itXsec==this->xsecVals.cend()) continue;
      auto itCachedSampleNormalizationsPerBin = this->cachedSampleNormalizationsPerBin.find(tree); if (itCachedSampleNormalizationsPerBin==this->cachedSampleNormalizationsPerBin.end()) continue;

      float const& sumwgts = itSumWeights->second.at(bin); // S_tj = sum{j} (w_jt * r_jt)
      float const& sumsqwgts = itSumSqWeights->second.at(bin); // V_tj = sum{j} ((w_jt * r_jt)^2)
      int const& nevtsnonzerowgt = itSumNonZeroWgtEvents->second.at(bin); // Indicator to show whether this bin is usable, should not enter into the calculation
      //float const& sumAllNominalWeights = itSumNominalWeights->second; // W_t from the sample
      float const& xsec = itXsec->second; // sigma_t
      float& sampleNormalizationsPerBinCache = itCachedSampleNormalizationsPerBin->second.at(bin);

      // W_t = sum{j} (W_tj)
      float sumAllNonZeroWgtNominalWeights = 0;
      for (float const& s:itSumNonZeroWgtNominalWeights->second) sumAllNonZeroWgtNominalWeights += s; // Summation over the W_tj

      float numerator_pertree = 0;
      float denominator_pertree = 0;
      if (nevtsnonzerowgt>0){
        // Sample weight
        if (ReweightingBuilder::useNeffInNormComponent) denominator_pertree = pow(sumwgts, 2)/sumsqwgts;
        else denominator_pertree = 1./sumsqwgts;

        // S_tj / W_t * sigma_t * sample weight
        sampleNormalizationsPerBinCache = xsec/sumAllNonZeroWgtNominalWeights*denominator_pertree;
        numerator_pertree = sumwgts * sampleNormalizationsPerBinCache;
      }

      // Sum numerator and denominator over samples t
      numerator += numerator_pertree;
      denominator += denominator_pertree;
    } // End loop over trees


    if (denominator!=0.){
      // Cache norm component for each bin
      cachedNormComponentPerBin = numerator/denominator;
      // Cache sample weights
      for (auto& tree:trees){
        auto itCachedSampleNormalizationsPerBin = this->cachedSampleNormalizationsPerBin.find(tree); if (itCachedSampleNormalizationsPerBin!=this->cachedSampleNormalizationsPerBin.end()) continue;
        float& sampleNormalizationsPerBinCache = itCachedSampleNormalizationsPerBin->second.at(bin);
        sampleNormalizationsPerBinCache /= denominator;
      }
    }
  } // End loop over bins
}


std::vector<BaseTree*> ReweightingBuilder::getRegisteredTrees() const{
  std::vector<BaseTree*> res;
  for (auto it=componentRefs_nominalweights.cbegin(); it!=componentRefs_nominalweights.cend(); it++) res.push_back(it->first);
  return res;
}


std::vector<float> ReweightingBuilder::getWeightThresholds(BaseTree* theTree) const{
  if (!theTree) return vector<float>();
  auto it = weightThresholds.find(theTree);

  if (it!=weightThresholds.cend()) return it->second;
  else return vector<float>();
}
float ReweightingBuilder::getPostThresholdWeight(BaseTree* theTree) const{
  float weight = this->eval_nominalweights(theTree) * this->eval_reweightingweights(theTree);
  float weightRel = weight-weightThresholdReference;
  int bin=this->findBin(theTree);
  if (bin>=0){
    const float& threshold=weightThresholds.find(theTree)->second.at(bin);
    if (threshold>=0. && fabs(weightRel)>threshold){
      weightRel = pow(threshold, 2)/weightRel;
      weight = weightRel + weightThresholdReference;
    }
  }
  return weight;
}
float ReweightingBuilder::getFinalEventWeight(BaseTree* theTree) const{
  int bin=this->findBin(theTree);
  if (bin<0) return 0;

  float weight = this->getPostThresholdWeight(theTree);

  // Get cached sample normalization factor
  auto itCachedSampleNormalizationsPerBin = this->cachedSampleNormalizationsPerBin.find(theTree);
  if (itCachedSampleNormalizationsPerBin==this->cachedSampleNormalizationsPerBin.cend()){
    MELAerr << "ReweightingBuilder::getFinalEventWeight: You must set up the caches first!" << endl;
    assert(0);
  }
  weight *= itCachedSampleNormalizationsPerBin->second.at(bin);

  return weight;
}

float ReweightingBuilder::getNormComponent(BaseTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->ReweightingBuilder::getNormComponent(bin);
}
float ReweightingBuilder::getNormComponent(int bin) const{
  if (bin<0) return 0;
  if (cachedNormComponentsPerBin.empty()){
    MELAerr << "ReweightingBuilder::getNormComponent: You must set up the caches first!" << endl;
    assert(0);
    return 0;
  }
  return cachedNormComponentsPerBin.at(bin);
}

float ReweightingBuilder::getBareSumNominalWeights(BaseTree* theTree) const{
  auto itSumNominalWeights = this->sumNominalWeights.find(theTree);
  if (itSumNominalWeights==this->sumNominalWeights.cend()) return 0;
  else return itSumNominalWeights->second;
}
