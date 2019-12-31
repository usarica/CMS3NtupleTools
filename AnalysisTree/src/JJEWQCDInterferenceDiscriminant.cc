#include "JJEWQCDInterferenceDiscriminant.h"
#include "HelperFunctions.h"
#include <cassert>


using namespace std;
using namespace HelperFunctions;


JJEWQCDInterferenceDiscriminant::JJEWQCDInterferenceDiscriminant(
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename,
  const float gscale_
) : Discriminant(cfilename, splinename, gfilename, gsplinename, gscale_){}

void JJEWQCDInterferenceDiscriminant::eval(const std::vector<float>& vars, const float& valReco){
  const unsigned int nvarsreq=21;
  this->resetVal();
  if (checkNanInf(vars) && checkNonNegative(vars) && vars.size()==nvarsreq){
    float constant = getCval(valReco);

    float vbf = vars[0]/vars[11];
    float zh = vars[1]/vars[12];
    float wh = vars[2]/vars[13];
    float constA = 1./(1./vars[11]+1./vars[12]+1./vars[13]);

    float vbs = vars[3]/vars[14];
    float zzz = vars[4]/vars[15];
    float wzz = vars[5]/vars[16];
    float qcdzz = vars[6]/vars[17];
    float constB = 1./(1./vars[14]+1./vars[15]+1./vars[16]+1./vars[17]);

    float vbf_vbs_int = vars[18]*(1./vars[11]+1./vars[14]) - vbf - vbs;
    float zh_zzz_int = vars[19]*(1./vars[12]+1./vars[15]) - zh - zzz;
    float wh_wzz_int = vars[20]*(1./vars[13]+1./vars[16]) - wh - wzz;

    const float scale_Pmjj_vb=1;
    float scale_Pmjj_z = vars[7]/vars[8];
    float scale_Pmjj_w = vars[9]/vars[10];

    vbf *= scale_Pmjj_vb;
    vbs *= scale_Pmjj_vb;
    vbf_vbs_int *= scale_Pmjj_vb;

    zh *= scale_Pmjj_z;
    zzz *= scale_Pmjj_z;
    zh_zzz_int *= scale_Pmjj_z;

    wh *= scale_Pmjj_w;
    wzz *= scale_Pmjj_w;
    wh_wzz_int *= scale_Pmjj_w;

    float PA = (vbf + zh + wh)*constA;
    float PB = (vbs + zzz + wzz + qcdzz)*constB;
    float Pint = (vbf_vbs_int + zh_zzz_int + wh_wzz_int)*sqrt(constA*constB);
    val = Pint*sqrt(constant)/(PA+constant*PB);
  }
}
