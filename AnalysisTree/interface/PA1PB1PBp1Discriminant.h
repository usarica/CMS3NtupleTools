#ifndef PA1PB1PBP1DISCRIMINANT_H
#define PA1PB1PBP1DISCRIMINANT_H

#include "Discriminant.h"


class PA1PB1PBp1Discriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  PA1PB1PBp1Discriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};


#endif
