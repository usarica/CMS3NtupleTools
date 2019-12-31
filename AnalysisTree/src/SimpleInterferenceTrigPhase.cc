#include "SimpleInterferenceTrigPhase.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


SimpleInterferenceTrigPhase::SimpleInterferenceTrigPhase() : Discriminant("", "", "", "", 1){}

void SimpleInterferenceTrigPhase::eval(const std::vector<float>& vars, const float& /*valReco*/){
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
    if (!isWithAvgME){
      val=0;
      float d=vars[0]*vars[1];
      if (d>0.) val += vars[2]/(2.*sqrt(d));
    }
    else{
      val=0;
      float pA = vars[0]/vars[3];
      float pB = vars[1]/vars[4];
      float pAB = vars[2]*(1./vars[3]+1./vars[4]);
      float pABInt = pAB-pA-pB;
      float d=pA*pB;
      if (d>0.) val += pABInt/(2.*sqrt(d));
    }
  }
}
