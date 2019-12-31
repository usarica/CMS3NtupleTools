#ifndef SIMPLEINTERFERENCEDISCRIMINANT_H
#define SIMPLEINTERFERENCEDISCRIMINANT_H

#include "Discriminant.h"


class SimpleInterferenceDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  SimpleInterferenceDiscriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef SimpleInterferenceDiscriminant Dintkin_t;

#endif
