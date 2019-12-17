#include "ReweightingFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


std::vector<float*> ReweightingFunctions::getWeightRefs(BaseTree* tree, const std::vector<TString>& strWeights){
  std::vector<float*> res;
  if (!tree || strWeights.empty()) return res;
  res.reserve(strWeights.size());
  for (auto const& s:strWeights){
    float* v = nullptr;
    tree->getValRef<float>(s, v);
    if (!v) MELAerr << "ReweightingFunctions::getWeightRefs: Could not get the reference for weight " << s << endl;
    res.push_back(v);
  }
  return res;
}
float* ReweightingFunctions::getWeightRef(BaseTree* tree, const TString& strWeight){
  float* res=nullptr;
  if (!tree || strWeight=="") return res;
  tree->getValRef<float>(strWeight, res);
  if (!res) MELAerr << "ReweightingFunctions::getWeightRef: Could not get the reference for weight " << strWeight << endl;
  return res;
}

float ReweightingFunctions::getSimpleWeight(BaseTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  res=1;
  for (auto const& v:vals){ if (v) res *= *v; }
  return res;
}
float ReweightingFunctions::getA1PlusB1Weight(BaseTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2);
  res = *(vals.at(0)) + *(vals.at(1));
  return res;
}
float ReweightingFunctions::getA1MinusB1Weight(BaseTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2);
  res = *(vals.at(0)) - *(vals.at(1));
  return res;
}
float ReweightingFunctions::getOnePlusB1OverA1Weight(BaseTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2);
  res = 1. + *(vals.at(1)) / *(vals.at(0));
  return res;
}
float ReweightingFunctions::getOneMinusB1OverA1Weight(BaseTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  assert(vals.size()==2);
  res = 1. - *(vals.at(1)) / *(vals.at(0));
  return res;
}
float ReweightingFunctions::getA1OverB1Weight(BaseTree* tree, const std::vector<float*>& vals){
  float res=0;
  if (!tree) return res;
  res=1;
  if (vals.size()!=2 || !(vals.at(0) && vals.at(1))){
    MELAerr << "ReweightingFunctions::getA1OverB1Weight: (vals.size()=" << vals.size() << "!=2) || !(vals.at(0)=" << vals.at(0) << " && vals.at(1)=" << vals.at(1) << "). Returning " << res << "..." << endl;
    return res;
  }
  res = (*(vals.at(0)))/(*(vals.at(1)));
  return res;
}
