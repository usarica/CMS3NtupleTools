#include "SimpleDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


SimpleDiscriminant::SimpleDiscriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void SimpleDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=2;
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==nvarsreq){
    float constant = getCval(valReco);
    val = vars[0]/(vars[0]+constant*vars[1]);
  }
}
