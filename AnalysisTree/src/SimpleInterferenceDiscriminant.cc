#include "SimpleInterferenceDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


SimpleInterferenceDiscriminant::SimpleInterferenceDiscriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void SimpleInterferenceDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=3;
  const unsigned int nvarsreq_withAvgME=5;
  this->resetVal();
  bool isWithAvgME=(vars.size()==nvarsreq_withAvgME);
  if (
    checkNanInf(vars) && (
    (isWithAvgME && checkNonNegative(vars))
    ||
    (vars.size()==nvarsreq && checkNonNegative(vars, -1, 2))
    )
    ){
    float constant = getCval(valReco);
    if (!isWithAvgME) val = vars[2]*sqrt(constant)/(vars[0]+constant*vars[1]);
    else val = (vars[2]*(1./vars[3]+1./vars[4])-vars[0]/vars[3]-vars[1]/vars[4])*sqrt(vars[3]*vars[4]*constant)/(vars[0]+constant*vars[1]);
  }
}
