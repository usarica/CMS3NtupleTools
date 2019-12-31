#ifndef JJEWQCDINTERFERENCEDISCRIMINANT_H
#define JJEWQCDINTERFERENCEDISCRIMINANT_H

#include "Discriminant.h"


class JJEWQCDInterferenceDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  JJEWQCDInterferenceDiscriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef JJEWQCDInterferenceDiscriminant DintjjEWQCD_t;

#endif
