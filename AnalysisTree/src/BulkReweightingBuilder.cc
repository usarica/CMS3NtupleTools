#include "BulkReweightingBuilder.h"
#include "SimpleEntry.h"
#include "MELAAccumulators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


ReweightingSpecifications::ReweightingSpecifications() :
  allowNegativeWeights(true),
  weightThresholdReference(0),
  rule_reweightingweights(nullptr)
{}
ReweightingSpecifications::ReweightingSpecifications(
  bool allowNegativeWeights_,
  std::vector<TString> const& strReweightingWeights_,
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_
) :
  allowNegativeWeights(allowNegativeWeights_),
  weightThresholdReference(0),
  strReweightingWeights(strReweightingWeights_),
  rule_reweightingweights(rule_reweightingweights_)
{}
ReweightingSpecifications::ReweightingSpecifications(ReweightingSpecifications const& other) :
  allowNegativeWeights(other.allowNegativeWeights),
  weightThresholdReference(other.weightThresholdReference),

  strReweightingWeights(other.strReweightingWeights),
  rule_reweightingweights(other.rule_reweightingweights),
  componentRefs_reweightingweights(other.componentRefs_reweightingweights),

  weightThresholds(other.weightThresholds),
  sumPostThrWeights(other.sumPostThrWeights),
  sumPostThrSqWeights(other.sumPostThrSqWeights),
  sumNonZeroWgtEvents(other.sumNonZeroWgtEvents),
  sumNonZeroWgtNominalWeights(other.sumNonZeroWgtNominalWeights),

  cachedSampleNormalizationsPerBin(other.cachedSampleNormalizationsPerBin),
  cachedNormComponentsPerBin(other.cachedNormComponentsPerBin),

  indexList(other.indexList)
{}
void ReweightingSpecifications::swap(ReweightingSpecifications& other){
  std::swap(allowNegativeWeights, other.allowNegativeWeights);
  std::swap(weightThresholdReference, other.weightThresholdReference);

  std::swap(strReweightingWeights, other.strReweightingWeights);
  std::swap(rule_reweightingweights, other.rule_reweightingweights);
  std::swap(componentRefs_reweightingweights, other.componentRefs_reweightingweights);

  std::swap(weightThresholds, other.weightThresholds);
  std::swap(sumPostThrWeights, other.sumPostThrWeights);
  std::swap(sumPostThrSqWeights, other.sumPostThrSqWeights);
  std::swap(sumNonZeroWgtEvents, other.sumNonZeroWgtEvents);
  std::swap(sumNonZeroWgtNominalWeights, other.sumNonZeroWgtNominalWeights);

  std::swap(cachedSampleNormalizationsPerBin, other.cachedSampleNormalizationsPerBin);
  std::swap(cachedNormComponentsPerBin, other.cachedNormComponentsPerBin);

  std::swap(indexList, other.indexList);
}
ReweightingSpecifications& ReweightingSpecifications::operator=(const ReweightingSpecifications& other){
  ReweightingSpecifications tmp(other);
  swap(tmp);
  return *this;
}


void ReweightingSpecifications::initialize(BaseTree* theTree, unsigned int ns, bool doReweighting){
  componentRefs_reweightingweights[theTree] = ReweightingFunctions::getWeightRefs(theTree, strReweightingWeights);

  sumPostThrWeights[theTree]=std::vector<float>(ns, 0);
  sumPostThrSqWeights[theTree]=std::vector<float>(ns, 0);
  weightThresholds[theTree]=std::vector<float>();
  sumNonZeroWgtEvents[theTree]=std::vector<unsigned int>();
  sumNonZeroWgtNominalWeights[theTree]=std::vector<float>();
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

  indexList.assign(ns, std::vector<SimpleEntry>());
  for (auto& index:indexList) index.reserve(theTree->getNEvents()/ns+1);
}
void ReweightingSpecifications::setupThresholds(BaseTree* theTree, float fractionRequirement, unsigned int minimumNevents){
  unsigned int ns = indexList.size();

  for (unsigned int ibin=0; ibin<ns; ibin++){
    std::vector<SimpleEntry>& index = indexList.at(ibin);
    if (minimumNevents>sumNonZeroWgtEvents[theTree].at(ibin)){
      MELAout << "\t- Bin " << ibin << " has " << sumNonZeroWgtEvents[theTree].at(ibin) << " events with non-zero weight, which is less than the requested number " << minimumNevents << ". Resetting the bin..." << endl;
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

  indexList.clear();
}
float ReweightingSpecifications::eval_reweightingweights(BaseTree* theTree) const{
  std::unordered_map<BaseTree*, std::vector<float*>>::const_iterator it = this->componentRefs_reweightingweights.find(theTree);
  if (it==this->componentRefs_reweightingweights.cend()){
    MELAerr << "ReweightingSpecifications::eval_reweightingweights: Could not find the weights " << this->strReweightingWeights << ". Call BulkReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }
  /*
  if (it->second.empty()){
    MELAerr << "ReweightingSpecifications::eval_reweightingweights: There are no component references!" << endl;
    return 0;
  }
  */

  float weight=0;
  if (this->rule_reweightingweights) weight = this->rule_reweightingweights(theTree, it->second);
  else MELAerr << "ReweightingSpecifications::eval_reweightingweights: Rule is not set!" << endl;
  if (!checkVarNanInf(weight)){
    MELAerr << "ReweightingSpecifications::eval_reweightingweights: Weight " << weight << " is nan/inf!" << endl;
    weight=0;
  }
  if (!allowNegativeWeights && weight<0.f){
    weight=0;
    MELAerr << "ReweightingSpecifications::eval_reweightingweights: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
  return weight;
}


BulkReweightingBuilder::BulkReweightingBuilder(
  std::vector<TString> const& strNominalWeights_,
  std::vector<TString> const& strCrossSectionWeights_,
  ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
) :
  allowNegativeWeights(true),

  strNominalWeights(strNominalWeights_),
  rule_nominalweights(rule_nominalweights_),

  strCrossSectionWeights(strCrossSectionWeights_),
  rule_xsecweights(rule_xsecweights_),

  currentReweightingSpecs(nullptr)
{}


float BulkReweightingBuilder::eval_nominalweights(BaseTree* theTree) const{
  std::unordered_map<BaseTree*, std::vector<float*>>::const_iterator it = componentRefs_nominalweights.find(theTree);
  if (it==componentRefs_nominalweights.cend()){
    MELAerr << "BulkReweightingBuilder::eval_nominalweights: Could not find the weights " << strNominalWeights << ". Call BulkReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }

  float weight=0;
  if (rule_nominalweights) weight=rule_nominalweights(theTree, it->second);
  if (!checkVarNanInf(weight)){
    MELAerr << "BulkReweightingBuilder::eval_nominalweights: Weight " << weight << " is nan/inf!" << endl;
    weight=0;
  }
  if (!allowNegativeWeights && weight<0.){
    weight=0;
    MELAerr << "BulkReweightingBuilder::eval_nominalweights: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
  return weight;
}
float BulkReweightingBuilder::eval_reweightingweights(BaseTree* theTree, TString strScheme) const{
  auto scheme = getReweightingSpecs(strScheme);
  if (!scheme) return 0;
  return scheme->eval_reweightingweights(theTree);
}
float BulkReweightingBuilder::eval_xsecweights(BaseTree* theTree) const{
  std::unordered_map<BaseTree*, std::vector<float*>>::const_iterator it = componentRefs_xsecweights.find(theTree);
  if (it==componentRefs_xsecweights.cend()){
    MELAerr << "BulkReweightingBuilder::eval_xsecweights: Could not find the weights " << strCrossSectionWeights << ". Call BulkReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }

  float weight=0;
  if (rule_xsecweights) weight=rule_xsecweights(theTree, it->second);
  if (!checkVarNanInf(weight)){
    MELAerr << "BulkReweightingBuilder::eval_xsecweights: Weight " << weight << " is nan/inf!" << endl;
    weight=0;
  }
  if (!allowNegativeWeights && weight<0.){
    weight=0;
    MELAerr << "BulkReweightingBuilder::eval_xsecweights: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
  return weight;
}

int BulkReweightingBuilder::findBin(BaseTree* theTree) const{
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
ReweightingSpecifications* BulkReweightingBuilder::getReweightingSpecs(TString strScheme){
  auto it = reweightingSpecs.find(strScheme);
  if (it==reweightingSpecs.end()){
    MELAerr << "BulkReweightingBuilder::getReweightingSpecs: Scheme " << strScheme << " is not set." << endl;
    return nullptr;
  }
  else return &(it->second);
}
ReweightingSpecifications const* BulkReweightingBuilder::getReweightingSpecs(TString strScheme) const{
  auto it = reweightingSpecs.find(strScheme);
  if (it==reweightingSpecs.cend()){
    MELAerr << "BulkReweightingBuilder::getReweightingSpecs: Scheme " << strScheme << " is not set." << endl;
    return nullptr;
  }
  else return &(it->second);
}
bool BulkReweightingBuilder::setCurrentReweightingSpecs(TString strScheme){ currentReweightingSpecs = getReweightingSpecs(strScheme); return (currentReweightingSpecs!=nullptr); }
void BulkReweightingBuilder::rejectNegativeWeights(const bool flag){
  allowNegativeWeights = !flag;
  for (auto& it:reweightingSpecs) it.second.allowNegativeWeights = allowNegativeWeights;
}
void BulkReweightingBuilder::setWeightBinning(const ExtendedBinning& binning){ weightBinning = binning; }
void BulkReweightingBuilder::setWeightThresholdReference(const float& weightThresholdReference_, TString strScheme){
  auto scheme = getReweightingSpecs(strScheme);
  if (scheme) scheme->weightThresholdReference = weightThresholdReference_;
}
void BulkReweightingBuilder::addReweightingWeights(
  TString strScheme,
  std::vector<TString> const& strReweightingWeights_,
  ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_
){
  auto hasScheme = getReweightingSpecs(strScheme);
  if (hasScheme){
    MELAerr << "BulkReweightingBuilder::addReweightingWeights: Scheme " << strScheme << " already exists!" << endl;
    return;
  }
  else reweightingSpecs[strScheme] = ReweightingSpecifications(allowNegativeWeights, strReweightingWeights_, rule_reweightingweights_);
}

void BulkReweightingBuilder::setupWeightVariables(BaseTree* theTree, float fractionRequirement, unsigned int minimumNevents){
  MELAout << "BulkReweightingBuilder::setupWeightVariables(\n"
    << "\t- Nominal weights: " << strNominalWeights << "\n\t- Reweighting weights: {" << endl;
  for (auto const& it:reweightingSpecs) MELAout << "\t\t[ " << it.first << ": " << it.second.strReweightingWeights << " ]" << endl;
  MELAout << "\t}\n\t- Xsec weights: " << strCrossSectionWeights << "\n) is called for tree " << theTree << "." << endl;

  if (!theTree) return;
  const int nevents = theTree->getNEvents();

  const ExtendedBinning& binning = this->weightBinning;
  const bool noBoundaries = !binning.isValid();
  const unsigned int ns = (!noBoundaries ? binning.getNbins() : 1);

  TString strOrderingVal = binning.getLabel();
  if (strOrderingVal=="") return;
  float* orderingValRef = ReweightingFunctions::getWeightRef(theTree, strOrderingVal);
  if (!orderingValRef) return;
  binningVarRefs[theTree] = orderingValRef;

  bool doReweighting = (fractionRequirement>=0.);

  // Link components and initialize
  componentRefs_nominalweights[theTree] = ReweightingFunctions::getWeightRefs(theTree, strNominalWeights);
  componentRefs_xsecweights[theTree] = ReweightingFunctions::getWeightRefs(theTree, strCrossSectionWeights);
  sumNominalWeights[theTree] = 0;
  xsecVals[theTree]=1;
  currentReweightingSpecs = nullptr;
  for (auto& it:reweightingSpecs) it.second.initialize(theTree, ns, doReweighting);

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
        float const& orderingVal = *orderingValRef;
        int bin=0;
        if (!noBoundaries) bin = binning.getBin(orderingVal);

        for (auto& it:reweightingSpecs){
          currentReweightingSpecs = &(it.second);

          float wgt = wgt_nominal * currentReweightingSpecs->eval_reweightingweights(theTree);
          SimpleEntry theEntry(0, fabs(wgt-currentReweightingSpecs->weightThresholdReference), wgt);

          if (bin>=0 && bin<(int) ns){
            // Accumulate the events
            currentReweightingSpecs->indexList.at(bin).push_back(theEntry);
            if (wgt!=0.f){
              currentReweightingSpecs->sumNonZeroWgtEvents[theTree].at(bin)++;
              currentReweightingSpecs->sumNonZeroWgtNominalWeights[theTree].at(bin) += wgt_nominal;
            }
          }
        }
        currentReweightingSpecs = nullptr; // Reset
      }

      if (firstTreeEvent) firstTreeEvent=false;
    }
    MELAout << "\t- Sum of nominal weights is computed as " << it_sumNominalWeights->second << "." << endl;
  }

  if (!doReweighting) return;

  for (auto& it:reweightingSpecs){
    MELAout << "=> Setting up thresholds for scheme " << it.first << ": " << endl;
    it.second.setupThresholds(theTree, fractionRequirement, minimumNevents);
  }
}
void BulkReweightingBuilder::setupCaches(){
  const ExtendedBinning& binning = weightBinning;
  int nbins;
  if (binning.isValid()) nbins = binning.getNbins();
  else nbins=1;

  std::vector<BaseTree*> const trees = this->getRegisteredTrees();

  currentReweightingSpecs = nullptr;
  for (auto& it_reweightingSpecs:reweightingSpecs){
    currentReweightingSpecs = &(it_reweightingSpecs.second);

    // Clear cached objects
    currentReweightingSpecs->cachedNormComponentsPerBin.clear();
    currentReweightingSpecs->cachedSampleNormalizationsPerBin.clear();

    // Setup default values
    currentReweightingSpecs->cachedNormComponentsPerBin.assign(nbins, 0);
    for (auto& tree:trees){
      currentReweightingSpecs->cachedSampleNormalizationsPerBin[tree] = std::vector<float>(nbins, 0);
    }

    for (int bin=0; bin<nbins; bin++){
      KahanAccumulator<float> numerator;
      KahanAccumulator<float> denominator;

      float& cachedNormComponentPerBin = currentReweightingSpecs->cachedNormComponentsPerBin.at(bin);

      for (auto& tree:trees){
        auto itSumWeights = currentReweightingSpecs->sumPostThrWeights.find(tree); if (itSumWeights==currentReweightingSpecs->sumPostThrWeights.cend()) continue;
        auto itSumSqWeights = currentReweightingSpecs->sumPostThrSqWeights.find(tree); if (itSumSqWeights==currentReweightingSpecs->sumPostThrSqWeights.cend()) continue;
        auto itSumNonZeroWgtNominalWeights = currentReweightingSpecs->sumNonZeroWgtNominalWeights.find(tree); if (itSumNonZeroWgtNominalWeights==currentReweightingSpecs->sumNonZeroWgtNominalWeights.cend()) continue;
        auto itSumNonZeroWgtEvents = currentReweightingSpecs->sumNonZeroWgtEvents.find(tree); if (itSumNonZeroWgtEvents==currentReweightingSpecs->sumNonZeroWgtEvents.cend()) continue;
        auto itSumNominalWeights = this->sumNominalWeights.find(tree); if (itSumNominalWeights==this->sumNominalWeights.cend()) continue;
        auto itXsec = this->xsecVals.find(tree); if (itXsec==this->xsecVals.cend()) continue;
        auto itCachedSampleNormalizationsPerBin = currentReweightingSpecs->cachedSampleNormalizationsPerBin.find(tree); if (itCachedSampleNormalizationsPerBin==currentReweightingSpecs->cachedSampleNormalizationsPerBin.end()) continue;

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
          if (BulkReweightingBuilder::useNeffInNormComponent) denominator_pertree = pow(sumwgts, 2)/sumsqwgts;
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
          auto itCachedSampleNormalizationsPerBin = currentReweightingSpecs->cachedSampleNormalizationsPerBin.find(tree); if (itCachedSampleNormalizationsPerBin!=currentReweightingSpecs->cachedSampleNormalizationsPerBin.end()) continue;
          float& sampleNormalizationsPerBinCache = itCachedSampleNormalizationsPerBin->second.at(bin);
          sampleNormalizationsPerBinCache /= denominator;
        }
      }
    } // End loop over bins
  } // End loop over schemes

  currentReweightingSpecs = nullptr;
}


std::vector<BaseTree*> BulkReweightingBuilder::getRegisteredTrees() const{
  std::vector<BaseTree*> res;
  for (auto it=componentRefs_nominalweights.cbegin(); it!=componentRefs_nominalweights.cend(); it++) res.push_back(it->first);
  return res;
}


std::vector<float> BulkReweightingBuilder::getWeightThresholds(BaseTree* theTree, TString strScheme) const{
  if (!theTree) return vector<float>();

  auto scheme = getReweightingSpecs(strScheme);
  if (!scheme) return vector<float>();

  auto it = scheme->weightThresholds.find(theTree);

  if (it!=scheme->weightThresholds.cend()) return it->second;
  else return vector<float>();
}
float BulkReweightingBuilder::getPostThresholdWeight(BaseTree* theTree, TString strScheme) const{
  auto scheme = getReweightingSpecs(strScheme);
  if (!scheme) return 0;

  float weight = this->eval_nominalweights(theTree) * scheme->eval_reweightingweights(theTree);
  float weightRel = weight-scheme->weightThresholdReference;
  int bin=this->findBin(theTree);
  if (bin>=0){
    const float& threshold=scheme->weightThresholds.find(theTree)->second.at(bin);
    if (threshold>=0. && fabs(weightRel)>threshold){
      weightRel = pow(threshold, 2)/weightRel;
      weight = weightRel + scheme->weightThresholdReference;
    }
  }
  return weight;
}
float BulkReweightingBuilder::getFinalEventWeight(BaseTree* theTree, TString strScheme) const{
  int bin=this->findBin(theTree);
  if (bin<0) return 0;

  auto scheme = getReweightingSpecs(strScheme);
  if (!scheme) return 0;

  float weight = this->getPostThresholdWeight(theTree, strScheme);

  // Get cached sample normalization factor
  auto itCachedSampleNormalizationsPerBin = scheme->cachedSampleNormalizationsPerBin.find(theTree);
  if (itCachedSampleNormalizationsPerBin==scheme->cachedSampleNormalizationsPerBin.cend()){
    MELAerr << "BulkReweightingBuilder::getFinalEventWeight: You must set up the caches first!" << endl;
    assert(0);
  }
  weight *= itCachedSampleNormalizationsPerBin->second.at(bin);

  return weight;
}

float BulkReweightingBuilder::getNormComponent(BaseTree* theTree, TString strScheme) const{
  int bin=this->findBin(theTree);
  return this->BulkReweightingBuilder::getNormComponent(bin, strScheme);
}
float BulkReweightingBuilder::getNormComponent(int bin, TString strScheme) const{
  if (bin<0) return 0;

  auto scheme = getReweightingSpecs(strScheme);
  if (!scheme) return 0;

  if (scheme->cachedNormComponentsPerBin.empty()){
    MELAerr << "BulkReweightingBuilder::getNormComponent: You must set up the caches first!" << endl;
    assert(0);
    return 0;
  }
  return scheme->cachedNormComponentsPerBin.at(bin);
}

float BulkReweightingBuilder::getBareSumNominalWeights(BaseTree* theTree) const{
  auto itSumNominalWeights = this->sumNominalWeights.find(theTree);
  if (itSumNominalWeights==this->sumNominalWeights.cend()) return 0;
  else return itSumNominalWeights->second;
}
