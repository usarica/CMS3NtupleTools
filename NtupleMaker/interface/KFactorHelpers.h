#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include "TFile.h"
#include "TSpline.h"


namespace KFactorHelpers{
  enum KFactorType{
    kf_QCD_NNLO_GGZZ_SIG,
    kf_QCD_NLO_GGZZ_SIG,

    kf_QCD_NNLO_QQZZ_BKG,
    kf_QCD_NNLO_QQWZ_BKG,
    kf_QCD_NNLO_QQWW_BKG,

    kf_EW_NLO_QQZZ_BKG,

    nKFactorTypes
  };

  class KFactorHandlerBase{
  protected:
    virtual void setup() = 0;

  public:
    KFactorHandlerBase(){}
    KFactorHandlerBase(KFactorHandlerBase const& other){}
    virtual ~KFactorHandlerBase(){}

    virtual void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) = 0;
    virtual void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map) = 0;
  };

  class KFactorHandler_QCD_ggZZ_Sig : public KFactorHandlerBase{
  public:
    static const std::string KFactorArgName;

  protected:
    std::shared_ptr<TFile> kfFile_NNLO;
    std::shared_ptr<TFile> kfFile_NLO;

    std::vector<TSpline3*> sp_NNLO;
    std::vector<TSpline3*> sp_NLO;

    void setup();

  public:
    KFactorHandler_QCD_ggZZ_Sig(int year);
    KFactorHandler_QCD_ggZZ_Sig(KFactorHandler_QCD_ggZZ_Sig const& other);
    ~KFactorHandler_QCD_ggZZ_Sig(){}

    void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map);
    void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map){ this->eval(type, KFactorHelpers::nKFactorTypes, kfactors_map); }

  };

}
