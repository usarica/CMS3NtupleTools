#ifndef VHPRODDECACGSDISCRIMINANT_H
#define VHPRODDECACGSDISCRIMINANT_H

#include "Discriminant.h"


class VHProdDecACGsDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  VHProdDecACGsDiscriminant(
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef VHProdDecACGsDiscriminant DaiGsVHdec_t;

#endif
