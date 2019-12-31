#include "VHProdDecACGsDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


VHProdDecACGsDiscriminant::VHProdDecACGsDiscriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void VHProdDecACGsDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  enum{
    iZHSM=0,
    iWHSM,
    iDecSM,
    iZHBSM,
    iDecBSM,
    iZHConst,
    iWHConst,
    iZHMJJ,
    iZHMJJtrue,
    iWHMJJ,
    iWHMJJtrue,
    nvarsreq
  };
  enum{
    iGVH, // There ought to be a dedicated g for WH+ZH
    iGDec,
    nGs
  };
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==(unsigned int) nvarsreq && theG.size()==(unsigned int)nGs){
    float pZHMJJwgt = vars[iZHMJJ]/vars[iZHMJJtrue];
    float pWHMJJwgt = vars[iWHMJJ]/vars[iWHMJJtrue];
    float pZHSM = vars[iZHSM]/vars[iZHConst]*pZHMJJwgt;
    float pWHSM = vars[iWHSM]/vars[iWHConst]*pWHMJJwgt;
    float pZHBSM = vars[iZHBSM]/vars[iZHConst]*pZHMJJwgt;
    constexpr float pWHBSM = 0;
    float pVHdecSM = (pZHSM + pWHSM)*vars[iDecSM];
    float pVHdecBSM = (pZHBSM + pWHBSM)*vars[iDecBSM];

    float gCommon = pow((theG.at(iGVH).second->Eval(valReco))*(theG.at(iGDec).second->Eval(valReco))*gscale, 2);

    val = pVHdecSM / (pVHdecSM + gCommon*pVHdecBSM);
  }
}
