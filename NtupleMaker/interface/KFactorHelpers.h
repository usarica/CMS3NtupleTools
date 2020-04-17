#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

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
    kf_EW_NLO_QQWZ_BKG,
    kf_EW_NLO_QQWW_BKG,

    nKFactorTypes
  };

  class KFactorHandlerBase{
  protected:
    virtual void setup() = 0;

  public:
    KFactorHandlerBase(){}
    KFactorHandlerBase(KFactorHandlerBase const& other){}
    virtual ~KFactorHandlerBase(){}

    virtual void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) const = 0;
    virtual void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map) const = 0;
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

    void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map) const;
    void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) const{ this->eval(type, KFactorHelpers::nKFactorTypes, kfactors_map); }

  };

  class KFactorHandler_EW_qqVV_Bkg : public KFactorHandlerBase{
  public:
    //static const std::string KFactorArgName;

  protected:
    // Raw table
    std::vector< std::vector<double> > table_ZZ;
    std::vector< std::vector<double> > table_WZ;
    std::vector< std::vector<double> > table_WW;

    // Locations of the first entries with a different shat
    std::vector< std::vector<std::vector<double>>::const_iterator > table_sqrtShatBegin_ZZ;
    std::vector< std::vector<std::vector<double>>::const_iterator > table_sqrtShatBegin_WZ;
    std::vector< std::vector<std::vector<double>>::const_iterator > table_sqrtShatBegin_WW;

    static void readTableFromFile(TString const& fname, std::vector< std::vector<double> >& table);
    static void findTableEntry(double const& shat, double const& that, std::vector< std::vector<double> > const& table, std::vector< std::vector<std::vector<double>>::const_iterator > const& table_sqrtShatBegin, double& val_uc, double& val_ds, double& val_b);

    void setup();

    void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) const{}
    void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType /*denominator*/, std::unordered_map<std::string, float>& kfactors_map) const{ return this->eval(type, kfactors_map); }

  public:
    KFactorHandler_EW_qqVV_Bkg();
    KFactorHandler_EW_qqVV_Bkg(KFactorHandler_EW_qqVV_Bkg const& other);
    ~KFactorHandler_EW_qqVV_Bkg(){}

    void eval(KFactorHelpers::KFactorType type, pat::PackedGenParticleCollection const& genparticles, std::unordered_map<std::string, float>& kfactors_map) const;

  };

}
