#ifndef VHPRODDECACDISCRIMINANT_H
#define VHPRODDECACDISCRIMINANT_H

#include "Discriminant.h"


class VHProdDecACDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  VHProdDecACDiscriminant(
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef VHProdDecACDiscriminant DaiVHdec_t;

#endif
