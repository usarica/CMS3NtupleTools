#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "TFile.h"
#include "TSpline.h"


namespace KFactorHelpers{
  // Various enumerators
  enum KFactorType{
    kf_QCD_NNLO_GGVV_SIG,
    kf_QCD_NLO_GGVV_SIG,

    kf_QCD_NNLO_QQZZ_BKG,
    kf_QCD_NNLO_QQWZ_BKG,
    kf_QCD_NNLO_QQWW_BKG,

    kf_EW_NLO_QQZZ_BKG,
    kf_EW_NLO_QQWZ_BKG,
    kf_EW_NLO_QQWW_BKG,

    nKFactorTypes
  };
  enum VVFinalStateType{
    kZZ,
    kWZ,
    kWW,
    nVVFinalStateTypes
  };

  typedef double(*KFactorFcn_qqZZ_QCD_t)(double, bool, unsigned char);
  typedef double(*KFactorFcn_qqWW_QCD_t)(double, unsigned char);
  typedef double(*KFactorFcn_qqWZ_QCD_t)(bool, unsigned char);


  // Functions that can be useful for different purposes
  float getBeamEnergy(int const& year);

  void getVVTopology(
    VVFinalStateType const& type,
    std::vector<reco::GenParticle> const& genparticles,
    std::vector<reco::GenParticle const*>& incomingQuarks,
    std::vector<reco::GenParticle const*>& incomingGluons,
    std::vector<reco::GenParticle const*>& outgoingQuarks,
    std::vector<reco::GenParticle const*>& outgoingGluons,
    std::pair<reco::GenParticle const*, reco::GenParticle const*>& V1pair,
    std::pair<reco::GenParticle const*, reco::GenParticle const*>& V2pair
  );


  // K factor classes
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

  class KFactorHandler_QCD_ggVV_Sig : public KFactorHandlerBase{
  public:
    static const std::string KFactorArgName;

  protected:
    std::shared_ptr<TFile> kfFile_NNLO;
    std::shared_ptr<TFile> kfFile_NLO;

    std::vector<TSpline3*> sp_NNLO;
    std::vector<TSpline3*> sp_NLO;

    void setup();

  public:
    KFactorHandler_QCD_ggVV_Sig() : KFactorHandlerBase(){}
    KFactorHandler_QCD_ggVV_Sig(int const& year);
    KFactorHandler_QCD_ggVV_Sig(KFactorHandler_QCD_ggVV_Sig const& other);
    ~KFactorHandler_QCD_ggVV_Sig(){}

    void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map) const;
    void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) const{ this->eval(type, KFactorHelpers::nKFactorTypes, kfactors_map); }

  };

  class KFactorHandler_EW_qqVV_Bkg : public KFactorHandlerBase{
  protected:
    KFactorHelpers::KFactorType type;
    double beamEnergy;

    KFactorFcn_qqZZ_QCD_t fcn_Kfactor_QCD_qqZZ;
    KFactorFcn_qqWZ_QCD_t fcn_Kfactor_QCD_qqWZ;
    KFactorFcn_qqWW_QCD_t fcn_Kfactor_QCD_qqWW;

    // Raw table
    std::vector< std::vector<double> > table_VV;

    // Locations of the first entries with a different shat
    std::vector< std::vector<std::vector<double>>::const_iterator > table_sqrtShatBegin_VV;

    static void readTableFromFile(TString const& fname, std::vector< std::vector<double> >& table);

    std::vector<double> findTableEntry(double const& shat, double const& that) const;

    void setup();

    void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) const{}
    void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType /*denominator*/, std::unordered_map<std::string, float>& kfactors_map) const{ return this->eval(type, kfactors_map); }

  public:
    KFactorHandler_EW_qqVV_Bkg();
    KFactorHandler_EW_qqVV_Bkg(int const& year, KFactorHelpers::KFactorType const& type_);
    KFactorHandler_EW_qqVV_Bkg(KFactorHandler_EW_qqVV_Bkg const& other);
    ~KFactorHandler_EW_qqVV_Bkg(){}

    KFactorHelpers::KFactorType const& getType() const{ return type; }

    void eval(
      GenEventInfoProduct const& eventInfo,
      std::vector<reco::GenParticle const*> const& incomingQuarks,
      std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV1pair,
      std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV2pair,
      std::unordered_map<std::string, float>& kfactors_map
    ) const;

  };

  class KFactorHandler_QCD_qqVV_Bkg : public KFactorHandlerBase{
  protected:
    KFactorFcn_qqZZ_QCD_t fcn_Kfactor_QCD_qqZZ;
    KFactorFcn_qqWZ_QCD_t fcn_Kfactor_QCD_qqWZ;
    KFactorFcn_qqWW_QCD_t fcn_Kfactor_QCD_qqWW;

    void setup(){}
    void eval(KFactorHelpers::KFactorType type, std::unordered_map<std::string, float>& kfactors_map) const{}
    void eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType /*denominator*/, std::unordered_map<std::string, float>& kfactors_map) const{ return this->eval(type, kfactors_map); }

  public:
    KFactorHandler_QCD_qqVV_Bkg();
    KFactorHandler_QCD_qqVV_Bkg(int const& year);
    KFactorHandler_QCD_qqVV_Bkg(KFactorHandler_QCD_qqVV_Bkg const& other);
    ~KFactorHandler_QCD_qqVV_Bkg(){}

    void eval(
      KFactorHelpers::KFactorType type,
      std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV1pair,
      std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV2pair,
      std::unordered_map<std::string, float>& kfactors_map
    ) const;

  };

}
