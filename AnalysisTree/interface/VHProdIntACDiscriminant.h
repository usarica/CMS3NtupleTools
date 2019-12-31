#ifndef VHPRODINTACDISCRIMINANT_H
#define VHPRODINTACDISCRIMINANT_H

#include "Discriminant.h"


class VHProdIntACDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  VHProdIntACDiscriminant(
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef VHProdIntACDiscriminant DaiVHint_t;

#endif
