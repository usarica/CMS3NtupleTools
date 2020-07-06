/*
Adapted from CJLST/ZZAnalysis
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/PatCandidates/interface/Muon.h>

#include <CMS3/NtupleMaker/interface/RoccoR.h>

#include "TLorentzVector.h"
#include "TRandom3.h"

#include <vector>
#include <string>
#include <sstream>
#include <memory>


using namespace edm;
using namespace std;
using namespace reco;


class RochesterPATMuonCorrector : public edm::EDProducer{
public:
  explicit RochesterPATMuonCorrector(const edm::ParameterSet&);
  ~RochesterPATMuonCorrector(){}

private:
  virtual void beginJob(){}
  virtual void endJob(){}
  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  string identifier_;
  bool isMC_;
  edm::EDGetTokenT< edm::View<pat::Muon> > muonToken_;

  std::shared_ptr<RoccoR> calibrator;

};


RochesterPATMuonCorrector::RochesterPATMuonCorrector(const edm::ParameterSet& iConfig) :
  identifier_(iConfig.getParameter<string>("identifier")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  muonToken_(consumes< edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("src")))
{
  std::stringstream ss; ss << "CMS3/NtupleMaker/data/RochesterMuonCorrections/" << identifier_ << ".txt";
  edm::FileInPath corrPath(ss.str());

  calibrator = std::make_shared<RoccoR>(corrPath.fullPath());

  produces<pat::MuonCollection>();
}


void RochesterPATMuonCorrector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Input collection
  edm::Handle< edm::View<pat::Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  // Output collection
  auto result = std::make_unique<pat::MuonCollection>(); result->reserve(muonHandle->size());

  TRandom3 rand;
  for (View<pat::Muon>::const_iterator muon = muonHandle->begin(); muon != muonHandle->end(); muon++){
    pat::Muon mu(*muon); // Clone the muon. This is the single muon to be put into the resultant collection

    double oldpt = mu.pt();

    auto gen_particle = mu.genParticle();
    double scale_factor = 1;
    double scale_error = 0.;
    double smear_error = 0.;

    // Deterministic scale/smear
    rand.SetSeed(std::abs(static_cast<int>(sin(mu.phi())*100000)));
    double u = rand.Rndm();

    if (calibrator  && mu.muonBestTrackType() == 1 && oldpt <= 200.){
      int nl = mu.track()->hitPattern().trackerLayersWithMeasurement();

      ////Protection against muons with low number of layers, they are not used in the analysis anyway as we apply tight muon ID
      //if (isMC_ && nl > 5){
      if (isMC_){
        //cout << "pt = " << mu.pt() << " eta = " << mu.eta() << " phi = " << mu.phi() << "u = " << u << endl;
        /// ====== ON MC (correction plus smearing) =====
        if (gen_particle){
          scale_factor = calibrator->kSpreadMC(mu.charge(), oldpt, mu.eta(), mu.phi(), gen_particle->pt());
          smear_error = calibrator->kSpreadMCerror(mu.charge(), oldpt, mu.eta(), mu.phi(), gen_particle->pt());
        }
        else{
          scale_factor = calibrator->kSmearMC(mu.charge(), oldpt, mu.eta(), mu.phi(), nl, u);
          smear_error = calibrator->kSmearMCerror(mu.charge(), oldpt, mu.eta(), mu.phi(), nl, u);
        }

        scale_error = calibrator->kScaleDTerror(mu.charge(), oldpt, mu.eta(), mu.phi());//there is no scale for MC so calculate it pretending it is data
      }
      //else if (!isMC_ && nl > 5){
      else if (!isMC_){
        /// ====== ON DATA (correction only) =====
        if (oldpt>2. && std::abs(mu.eta())<2.4){
          scale_factor = calibrator->kScaleDT(mu.charge(), oldpt, mu.eta(), mu.phi());
          scale_error = calibrator->kScaleDTerror(mu.charge(), oldpt, mu.eta(), mu.phi());
          smear_error = calibrator->kSmearMCerror(mu.charge(), oldpt, mu.eta(), mu.phi(), nl, u);//there is no smear in data so calculate it pretending it is mc
        }
        else{
          // keep old values
          scale_factor = 1.;
          scale_error = 0.;
          smear_error = 0.;
        }
      }
    }

    mu.addUserFloat("scale_smear_pt_corr", scale_factor);
    mu.addUserFloat("scale_smear_pt_corr_scale_totalUp", scale_factor*(1. + scale_error));
    mu.addUserFloat("scale_smear_pt_corr_smear_totalUp", scale_factor*(1. + smear_error));
    mu.addUserFloat("scale_smear_pt_corr_scale_totalDn", scale_factor*(1. - scale_error));
    mu.addUserFloat("scale_smear_pt_corr_smear_totalDn", scale_factor*(1. - smear_error));

    result->push_back(mu);
  }

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(RochesterPATMuonCorrector);
