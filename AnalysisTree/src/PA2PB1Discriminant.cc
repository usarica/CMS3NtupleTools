#include "PA2PB1Discriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


PA2PB1Discriminant::PA2PB1Discriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void PA2PB1Discriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=3;
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==nvarsreq){
    float constant = getCval(valReco);
    val = vars[0]*vars[1]/(vars[0]*vars[1]+constant*vars[2]);
  }
}
