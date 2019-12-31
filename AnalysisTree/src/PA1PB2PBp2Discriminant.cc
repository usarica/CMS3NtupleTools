#include "PA1PB2PBp2Discriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


PA1PB2PBp2Discriminant::PA1PB2PBp2Discriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void PA1PB2PBp2Discriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=5;
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==nvarsreq){
    float constant = getCval(valReco);
    val = vars[0]/(vars[0]+constant*(vars[1]/vars[3] + vars[2]/vars[4])/(1./vars[3] + 1./vars[4]));
  }
}
