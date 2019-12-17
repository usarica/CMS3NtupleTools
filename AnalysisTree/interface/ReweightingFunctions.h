#ifndef REWEIGHTINGFUNCTIONS_H
#define REWEIGHTINGFUNCTIONS_H

#include "BaseTree.h"
#include <vector>
#include "TString.h"
#include "HelperFunctions.h"


namespace ReweightingFunctions{
  std::vector<float*> getWeightRefs(BaseTree* tree, const std::vector<TString>& strWeights);
  float* getWeightRef(BaseTree* tree, const TString& strWeight);

  typedef float(*ReweightingFunction_t)(BaseTree*, const std::vector<float*>&);
  float getSimpleWeight(BaseTree* tree, const std::vector<float*>& vals);
  float getA1PlusB1Weight(BaseTree* tree, const std::vector<float*>& vals);
  float getA1MinusB1Weight(BaseTree* tree, const std::vector<float*>& vals);
  float getOnePlusB1OverA1Weight(BaseTree* tree, const std::vector<float*>& vals);
  float getOneMinusB1OverA1Weight(BaseTree* tree, const std::vector<float*>& vals);
  float getA1OverB1Weight(BaseTree* tree, const std::vector<float*>& vals);

}


#endif
