#ifndef JJEWQCDBKGM4LDISCRIMINANT_H
#define JJEWQCDBKGM4LDISCRIMINANT_H

#include "Discriminant.h"


class JJEWQCDBkgM4LDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  JJEWQCDBkgM4LDiscriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef JJEWQCDBkgM4LDiscriminant Dbkgm4ljjEWQCD_t;

#endif
