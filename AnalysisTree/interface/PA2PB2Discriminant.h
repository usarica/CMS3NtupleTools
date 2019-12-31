#ifndef PA2PB2DISCRIMINANT_H
#define PA2PB2DISCRIMINANT_H

#include "Discriminant.h"


class PA2PB2Discriminant : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  PA2PB2Discriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
};

typedef PA2PB2Discriminant Dbkgm4l_t;
typedef PA2PB2Discriminant DjjVH_t;
typedef PA2PB2Discriminant DaiVBFdec_t;

#endif
