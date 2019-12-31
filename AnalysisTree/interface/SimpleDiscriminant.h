#ifndef SIMPLEDISCRIMINANT_H
#define SIMPLEDISCRIMINANT_H

#include "Discriminant.h"


class SimpleDiscriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  SimpleDiscriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef SimpleDiscriminant Dbkgkin_t;
typedef SimpleDiscriminant Dggbkgkin_t;
typedef SimpleDiscriminant DjjVBF_t;

#endif
