#ifndef DISCRIMINANT_H
#define DISCRIMINANT_H

#include <vector>
#include "TFile.h"
#include "TString.h"
#include "TSpline.h"


class Discriminant{
protected:
  std::vector<std::pair<TFile*, TSpline3*>> theC;
  std::vector<std::pair<TFile*, TSpline3*>> theG;

  float WPCshift;
  float gscale;
  bool invertG;

  float val;

  void resetVal();
  virtual void eval(const std::vector<float>& vars, const float& valReco)=0;

public:
  Discriminant(
    const TString cfilename="", const TString splinename="sp_gr_varReco_Constant_Smooth",
    const TString gfilename="", const TString gsplinename="",
    const float gscale_=1
  );
  virtual ~Discriminant();

  operator float() const;
  operator float&();
  operator float*();

  bool operator<(const float& other) const;
  bool operator>(const float& other) const;
  bool operator<=(const float& other) const;
  bool operator>=(const float& other) const;
  bool operator==(const float& other) const;
  bool operator!=(const float& other) const;

  virtual float getCval(const float valReco) const;
  float update(const std::vector<float>& vars, const float valReco);
  float applyAdditionalC(const float cval);

  void setWP(float inval=0.5);
  void setGScale(float inval=1);
  void setInvertG(bool flag);

  bool addAdditionalC(const TString filename, const TString splinename);
  bool addAdditionalG(const TString filename, const TString splinename);

};


#endif
