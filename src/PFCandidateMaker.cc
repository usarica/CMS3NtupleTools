#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CMS3/NtupleMaker/interface/PFCandidateMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

PFCandidateMaker::PFCandidateMaker(const edm::ParameterSet& iConfig){

  pfCandidatesTag_  = iConfig.getParameter<InputTag> ("pfCandidatesTag");

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<vector<LorentzVector> > (branchprefix+"p4"                   ).setBranchAlias(aliasprefix_+"_p4"                   );
  produces<vector<float> >         (branchprefix+"mass"                 ).setBranchAlias(aliasprefix_+"_mass"                 );
  produces<vector<float> >         (branchprefix+"dz"                   ).setBranchAlias(aliasprefix_+"_dz"                   );
  produces<vector<int> >           (branchprefix+"charge"               ).setBranchAlias(aliasprefix_+"_charge"               );
  produces<vector<int> >           (branchprefix+"particleId"           ).setBranchAlias(aliasprefix_+"_particleId"           );
  produces<vector<uint8_t> >       (branchprefix+"fromPV"               ).setBranchAlias(aliasprefix_+"_fromPV"               );
  produces<vector<bool> >          (branchprefix+"isStandAloneMuon"     ).setBranchAlias(aliasprefix_+"_isStandAloneMuon"     );
  produces<vector<bool> >          (branchprefix+"isGlobalMuon"         ).setBranchAlias(aliasprefix_+"_isGlobalMuon"         );
  produces<vector<uint8_t> >       (branchprefix+"pvAssociationQuality" ).setBranchAlias(aliasprefix_+"_pvAssociationQuality" );
  produces<vector<int> >           (branchprefix+"IdAssociatedPV"       ).setBranchAlias(aliasprefix_+"_IdAssociatedPV"       );
  produces<vector<float> >         (branchprefix+"dzAssociatedPV"       ).setBranchAlias(aliasprefix_+"_dzAssociatedPV"       );
  produces<vector<float> >         (branchprefix+"puppiWeight"          ).setBranchAlias(aliasprefix_+"_puppiWeight"          );

  produces<float>                  ("evtfixgridrhoctr" ).setBranchAlias("evt_fixgrid_rho_ctr");
  produces<float>                  ("evtfixgridrhofwd" ).setBranchAlias("evt_fixgrid_rho_fwd");
  produces<float>                  ("evtfixgridrhoall" ).setBranchAlias("evt_fixgrid_rho_all");
}

PFCandidateMaker::~PFCandidateMaker(){}
void  PFCandidateMaker::beginRun(const edm::Run&, const edm::EventSetup& es){}
void PFCandidateMaker::beginJob() {}
void PFCandidateMaker::endJob()   {}

void PFCandidateMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector <LorentzVector> > pfcands_p4                  (new vector<LorentzVector> );
  auto_ptr<vector <float> >         pfcands_mass                (new vector<float>         );
  auto_ptr<vector <float> >         pfcands_dz                  (new vector<float>         );
  auto_ptr<vector <int> >           pfcands_charge              (new vector<int>           );
  auto_ptr<vector <int> >           pfcands_particleId          (new vector<int>           );
  auto_ptr<vector <bool> >          pfcands_isStandAloneMuon    (new vector<bool>          );
  auto_ptr<vector <bool> >          pfcands_isGlobalMuon        (new vector<bool>          );
  auto_ptr<vector <uint8_t> >       pfcands_fromPV              (new vector<uint8_t>       );
  auto_ptr<vector<uint8_t> >        pfcands_pvAssociationQuality(new vector<uint8_t>       );
  auto_ptr<vector<int> >            pfcands_IdAssociatedPV      (new vector<int>           );
  auto_ptr<vector<float> >          pfcands_dzAssociatedPV      (new vector<float>         );
  auto_ptr<vector<float> >          pfcands_puppiWeight         (new vector<float>         );

  auto_ptr<float >                  evt_fixgrid_rho_ctr      (new float               );
  auto_ptr<float >                  evt_fixgrid_rho_fwd      (new float               );
  auto_ptr<float >                  evt_fixgrid_rho_all      (new float               );
    
  //get pfcandidates
  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  pfCandidates  = pfCandidatesHandle.product();

  for( pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {
	pfcands_p4   -> push_back( LorentzVector(pf_it->p4()) );
    pfcands_mass -> push_back( pf_it->mass()              );

	if (!pf_it->vertexRef().isNull()){
	  pfcands_dz                   -> push_back( pf_it->dz()                   );
	  pfcands_pvAssociationQuality -> push_back( pf_it->pvAssociationQuality() );
	  pfcands_dzAssociatedPV       -> push_back( pf_it->dzAssociatedPV()       );
	  pfcands_IdAssociatedPV       -> push_back( pf_it->vertexRef().key()      );
	}
	else {
	  pfcands_dz                   -> push_back( -9999.                        );
	  pfcands_pvAssociationQuality -> push_back( 0                             );
	  pfcands_dzAssociatedPV       -> push_back( -9999.                        );
	  pfcands_IdAssociatedPV       -> push_back( -9999                         );
	}

	pfcands_charge           -> push_back( pf_it->charge()           );
    pfcands_particleId       -> push_back( pf_it->pdgId()            );
    pfcands_fromPV           -> push_back( pf_it->fromPV()           );
	pfcands_puppiWeight      -> push_back( pf_it->puppiWeight()      );
    pfcands_isStandAloneMuon -> push_back( pf_it->isStandAloneMuon() );
    pfcands_isGlobalMuon     -> push_back( pf_it->isGlobalMuon()     );
  }//loop over candidate collection

  //define the phi bins
  vector<float> phibins;
  for (int i=0;i<10;i++) phibins.push_back(-TMath::Pi()+(2*i+1)*TMath::TwoPi()/20.);

  //define the eta bins
  vector<float> etabins_ctr;
  for (int i=0;i<8;++i) etabins_ctr.push_back(-2.1+0.6*i);
  vector<float> etabins_fwd;
  for (int i=0;i<10;++i) {
	if (i<5) etabins_fwd.push_back(-5.1+0.6*i);
	else etabins_fwd.push_back(2.7+0.6*(i-5));
  }
  vector<float> etabins_all;
  for (int i=0;i<18;++i) etabins_all.push_back(-5.1+0.6*i);

  //compute it
  *evt_fixgrid_rho_ctr = getFixGridRho(etabins_ctr,phibins);
  *evt_fixgrid_rho_fwd = getFixGridRho(etabins_fwd,phibins);
  *evt_fixgrid_rho_all = getFixGridRho(etabins_all,phibins);
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  //Keep it
  iEvent.put(pfcands_p4                   , branchprefix+"p4"                  );
  iEvent.put(pfcands_mass                 , branchprefix+"mass"                );
  iEvent.put(pfcands_dz                   , branchprefix+"dz"                  );
  iEvent.put(pfcands_charge               , branchprefix+"charge"              );
  iEvent.put(pfcands_particleId           , branchprefix+"particleId"          );
  iEvent.put(pfcands_isGlobalMuon         , branchprefix+"isGlobalMuon"        );
  iEvent.put(pfcands_isStandAloneMuon     , branchprefix+"isStandAloneMuon"    );
  iEvent.put(pfcands_fromPV               , branchprefix+"fromPV"              );
  iEvent.put(pfcands_pvAssociationQuality , branchprefix+"pvAssociationQuality");
  iEvent.put(pfcands_IdAssociatedPV       , branchprefix+"IdAssociatedPV"      );
  iEvent.put(pfcands_dzAssociatedPV       , branchprefix+"dzAssociatedPV"      );
  iEvent.put(pfcands_puppiWeight          , branchprefix+"puppiWeight"         );

  iEvent.put(evt_fixgrid_rho_ctr      , "evtfixgridrhoctr");
  iEvent.put(evt_fixgrid_rho_fwd      , "evtfixgridrhofwd");
  iEvent.put(evt_fixgrid_rho_all      , "evtfixgridrhoall");
 
}

float PFCandidateMaker::getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins) {
  float etadist = etabins[1]-etabins[0];
  float phidist = phibins[1]-phibins[0];
  float etahalfdist = (etabins[1]-etabins[0])/2.;
  float phihalfdist = (phibins[1]-phibins[0])/2.;
  vector<float> sumPFNallSMDQ;
  sumPFNallSMDQ.reserve(etabins.size()*phibins.size());
  for (unsigned int ieta=0;ieta<etabins.size();++ieta) {
	for (unsigned int iphi=0;iphi<phibins.size();++iphi) {
	  float pfniso_ieta_iphi = 0;
	  for(pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++) {
		if (fabs(etabins[ieta]-pf_it->eta())>etahalfdist) continue;
		if (fabs(reco::deltaPhi(phibins[iphi],pf_it->phi()))>phihalfdist) continue;
		pfniso_ieta_iphi+=pf_it->pt();
	  }
	  sumPFNallSMDQ.push_back(pfniso_ieta_iphi);
	}
  }
  float evt_smdq = 0;
  sort(sumPFNallSMDQ.begin(),sumPFNallSMDQ.end());
  if (sumPFNallSMDQ.size()%2) evt_smdq = sumPFNallSMDQ[(sumPFNallSMDQ.size()-1)/2];
  else evt_smdq = (sumPFNallSMDQ[sumPFNallSMDQ.size()/2]+sumPFNallSMDQ[(sumPFNallSMDQ.size()-2)/2])/2.;
  return evt_smdq/(etadist*phidist);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFCandidateMaker);
