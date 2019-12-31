#include "SimpleAverageInterferenceTrigPhase.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


SimpleAverageInterferenceTrigPhase::SimpleAverageInterferenceTrigPhase() : Discriminant("", "", "", "", 1){}

void SimpleAverageInterferenceTrigPhase::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=6;
  const unsigned int nvarsreq_withAvgME=10;
  this->resetVal();
  bool isWithAvgME=(vars.size()==nvarsreq_withAvgME);
  if (
    checkNanInf(vars) && (
      (isWithAvgME && checkNonNegative(vars))
      ||
      (vars.size()==nvarsreq && checkNonNegative(vars, -1, 2) && checkNonNegative(vars, 3, 5))
      )
    ){
    float constant = getCval(valReco);
    val=0;
    float d1=0, d2=0;
    if (!isWithAvgME){
      d1=vars[0]*vars[1];
      d2=vars[3]*vars[4];
      if (d1>0.) val += vars[2]/(2.*sqrt(d1));
      if (d2>0.) val += constant*vars[5]/(2.*sqrt(d2));
    }
    else{
      float pA = vars[0]/vars[3];
      float pB = vars[1]/vars[4];
      float pAB = vars[2]*(1./vars[3]+1./vars[4]);
      float pABInt = pAB-pA-pB;
      d1=pA*pB;
      if (d1>0.) val += pABInt/(2.*sqrt(d1));

      float pC = vars[5]/vars[8];
      float pD = vars[6]/vars[9];
      float pCD = vars[7]*(1./vars[8]+1./vars[9]);
      float pCDInt = pCD-pC-pD;
      d2=pC*pD;
      if (d2>0.) val += pCDInt/(2.*sqrt(d2));
    }
    if (d1>0. && d2>0.) val /= 2.;
  }
}
