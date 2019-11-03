// -*- C++ -*-
//
// Package:    MuonMaker
// Class:      MuonMaker
// 
/**\class MuonMaker MuonMaker.cc CMS2/MuonMaker/src/MuonMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonMaker.cc,v 1.68 2012/07/20 01:19:39 dbarge Exp $
//
//


// system include files
#include <memory>
#include <sstream>

// user include files
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CMS3/NtupleMaker/interface/plugins/MuonMaker.h"
#include "CMS3/NtupleMaker/interface/plugins/MatchUtilities.h"

#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"


//////////////
// typedefs //
//////////////

typedef math::XYZPoint Point;


////////////////
// namespaces //
////////////////

using namespace std;
using namespace reco;
using namespace edm;


/////////////////
// Constructor //
/////////////////

MuonMaker::MuonMaker( const ParameterSet& iConfig ) {

    /////////////////////////////
    // Branch & Alias prefixes //
    /////////////////////////////

    aliasprefix_        = iConfig.getUntrackedParameter<string>("aliasPrefix");
    branchprefix_       = aliasprefix_; if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


    //////////////////////
    // Input Parameters //
    //////////////////////

    muonsToken    = consumes<View<pat::Muon> >(iConfig.getParameter<InputTag> ("muonsInputTag"   ));
    vtxToken         = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));

    produces<pat::MuonCollection>().setBranchAlias(aliasprefix_);
} // end Constructor

void MuonMaker::beginJob () {}  // method called once each job just before starting event loop
void MuonMaker::endJob   () {}  // method called once each job just after ending the event loop


//////////////
// Producer //
//////////////

void MuonMaker::produce(Event& iEvent, const EventSetup& iSetup) {
  auto result = std::make_unique<pat::MuonCollection>();

    Handle<View<pat::Muon> > muon_h;
    iEvent.getByToken( muonsToken , muon_h );
    iEvent.getByToken( vtxToken , vertexHandle );

    size_t muonIndex = 0;
    for ( View<pat::Muon>::const_iterator muon = muon_h->begin(); muon != muon_h->end(); muon++, muonIndex++ ) {
    pat::Muon muon_result(*muon); // Clone the muon. This is the single muon to be put into the resultant collection

        // References
        const RefToBase<pat::Muon>    muonRef                 = muon_h->refAt(muonIndex); 
        const TrackRef                globalTrack             = muon->globalTrack();
        const TrackRef                siTrack                 = muon->innerTrack();
        const TrackRef                staTrack                = muon->outerTrack();
        const TrackRef                bestTrack               = muon->muonBestTrack();
        const MuonQuality             quality                 = muon->combinedQuality();
        const VertexCollection*       vertexCollection        = vertexHandle.product();

        // Iterators
        VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
        int firstGoodVertexIdx = 0;
        for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++firstGoodVertexIdx) {
            if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
                firstGoodVertex = vtx;
                break;
            }
        }

        muon_result.addUserInt("selectors", muon->selectors());
        muon_result.addUserInt("simType", muon->simType());
        muon_result.addUserInt("simExtType", muon->simExtType());

        ////////////
        // Global //
        ////////////

        muon_result.addUserFloat("gfit_chi2", globalTrack.isNonnull() ? globalTrack->chi2() : -9999. );
        muon_result.addUserInt("gfit_ndof", globalTrack.isNonnull() ? globalTrack->ndof() : -9999  );
        muon_result.addUserInt("gfit_validSTAHits", globalTrack.isNonnull() ? globalTrack->hitPattern().numberOfValidMuonHits() : -9999  );
        muon_result.addUserFloat("gfit_pt", globalTrack.isNonnull() ? globalTrack->pt() : 0.0 );
        muon_result.addUserFloat("gfit_eta", globalTrack.isNonnull() ? globalTrack->eta() : 0.0 );
        muon_result.addUserFloat("gfit_phi", globalTrack.isNonnull() ? globalTrack->phi() : 0.0 );
        muon_result.addUserInt("gfit_algo", globalTrack.isNonnull() ? globalTrack->algo() : -9999  );
        muon_result.addUserFloat("gfit_ptErr", globalTrack.isNonnull() ? globalTrack->ptError() : -9999. );

        //////////
        // Best //
        //////////

        muon_result.addUserFloat("bfit_pt", bestTrack.isNonnull() ? bestTrack->pt() : 0.0 );
        muon_result.addUserFloat("bfit_eta", bestTrack.isNonnull() ? bestTrack->eta() : 0.0 );
        muon_result.addUserFloat("bfit_phi", bestTrack.isNonnull() ? bestTrack->phi() : 0.0 );
        muon_result.addUserInt("bfit_algo", bestTrack.isNonnull() ? bestTrack->algo() : -9999  );
        muon_result.addUserFloat("bfit_ptErr", bestTrack.isNonnull() ? bestTrack->ptError() : -9999. );

        //////////////////
        // Muon Quality //
        //////////////////

        muon_result.addUserFloat("trkKink", quality.trkKink);
        muon_result.addUserFloat("chi2LocalPosition", quality.chi2LocalPosition);
        muon_result.addUserFloat("chi2LocalMomentum", quality.chi2LocalMomentum);


        //////////
        // Muon //
        //////////
    
        /////////////////////
        // Muon Quantities //
        /////////////////////


        muon_result.addUserInt("type", muon->type());
        muon_result.addUserInt("charge", muon->charge());
        muon_result.addUserFloat("caloCompatibility", muon->caloCompatibility());
        muon_result.addUserFloat("segmentCompatibility", muon::segmentCompatibility(*muon));
        muon_result.addUserInt("numberOfMatchedStations", muon::segmentCompatibility(*muon));
        muon_result.addUserFloat("pt", muon->pt());
        muon_result.addUserFloat("eta", muon->eta());
        muon_result.addUserFloat("phi", muon->phi());
        muon_result.addUserFloat("mass", muon->mass());

        ////////
        // ID //
        ////////

        bool matchIsValid = muon->isMatchesValid();

        muon_result.addUserInt("pid_TMLastStationLoose", matchIsValid ? muon::isGoodMuon( *muon, muon::TMLastStationLoose     ) : -9999  );
        muon_result.addUserInt("pid_TMLastStationTight", matchIsValid ? muon::isGoodMuon( *muon, muon::TMLastStationTight     ) : -9999  );
        muon_result.addUserInt("pid_TM2DCompatibilityLoose", matchIsValid ? muon::isGoodMuon( *muon, muon::TM2DCompatibilityLoose ) : -9999  );
        muon_result.addUserInt("pid_TM2DCompatibilityTight", matchIsValid ? muon::isGoodMuon( *muon, muon::TM2DCompatibilityTight ) : -9999  );
        muon_result.addUserInt("pid_TMOneStationTight", matchIsValid ? muon::isGoodMuon( *muon, muon::TMOneStationTight      ) : -9999  );
        muon_result.addUserInt("pid_PFMuon", muon->isPFMuon() );

        ////////////
        // Energy //
        ////////////

        bool energyIsValid = muon->isEnergyValid();

        muon_result.addUserFloat("ecal_time",  energyIsValid ? muon->calEnergy().ecal_time : -9999. );
        muon_result.addUserFloat("hcal_time",  energyIsValid ? muon->calEnergy().hcal_time : -9999. );

        ///////////////
        // Isolation //
        ///////////////

        muon_result.addUserFloat("iso_trckvetoDep", muon->isEnergyValid()    ? muon->isolationR03().trackerVetoPt  : -9999.        );
        muon_result.addUserFloat("iso_ecalvetoDep", muon->isEnergyValid()    ? muon->isolationR03().emVetoEt       : -9999.        );
        muon_result.addUserFloat("iso_hcalvetoDep", muon->isEnergyValid()    ? muon->isolationR03().hadVetoEt      : -9999.        );
        muon_result.addUserFloat("iso_hovetoDep", muon->isEnergyValid()    ? muon->isolationR03().hoVetoEt       : -9999.        );
        muon_result.addUserFloat("iso03_sumPt", muon->isIsolationValid() ? muon->isolationR03().sumPt          : -9999.        );
        muon_result.addUserFloat("iso03_emEt", muon->isIsolationValid() ? muon->isolationR03().emEt           : -9999.        );
        muon_result.addUserFloat("iso03_hadEt", muon->isIsolationValid() ? muon->isolationR03().hadEt          : -9999.        );
        muon_result.addUserInt("iso03_ntrk", muon->isIsolationValid() ? muon->isolationR03().nTracks        : -9999         );

        ////////////
        // Tracks //
        ////////////

        // 0imuon_result.addUserFloat("kjldf_.Ea",kjldt(x
        // "wp -- to paste macro
        // "wy -- to yank macro into w after selecting it

        muon_result.addUserFloat("trk_pt", siTrack.isNonnull()     ? siTrack.get()->pt() : 0.);
        muon_result.addUserFloat("trk_eta", siTrack.isNonnull()     ? siTrack.get()->eta() : 0.);
        muon_result.addUserFloat("trk_phi", siTrack.isNonnull()     ? siTrack.get()->phi() : 0.);
        muon_result.addUserInt("validHits", siTrack.isNonnull()     ? siTrack->numberOfValidHits()                         : -9999         );
        muon_result.addUserInt("lostHits", siTrack.isNonnull()     ? siTrack->numberOfLostHits()                          : -9999         );
        muon_result.addUserFloat("d0Err", siTrack.isNonnull()     ? siTrack->d0Error()                                   :  -9999.       );
        muon_result.addUserFloat("z0Err", siTrack.isNonnull()     ? siTrack->dzError()                                   :  -9999.       );
        muon_result.addUserFloat("ptErr", siTrack.isNonnull()     ? siTrack->ptError()                                   :  -9999.       );
        muon_result.addUserInt("algo", siTrack.isNonnull()     ? siTrack->algo       ()                               : -9999.        );
        muon_result.addUserInt("algoOrig", siTrack.isNonnull()     ? siTrack->originalAlgo       ()                               : -9999.        );
        muon_result.addUserInt("nlayers", siTrack.isNonnull()     ? siTrack->hitPattern().trackerLayersWithMeasurement() :  -9999        );
        muon_result.addUserInt("validPixelHits", siTrack.isNonnull()     ? siTrack->hitPattern().numberOfValidPixelHits()       :  -9999        );
        muon_result.addUserInt("exp_innerlayers", siTrack.isNonnull()     ? siTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)   :  -9999        );
        muon_result.addUserInt("exp_outerlayers", siTrack.isNonnull()     ? siTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS)   :  -9999        );

        if (firstGoodVertex!=vertexCollection->end()) { 
            muon_result.addUserFloat("dxyPV", siTrack.isNonnull()     ? siTrack->dxy( firstGoodVertex->position() )           : -9999.        );
            muon_result.addUserFloat("dzPV", siTrack.isNonnull()     ? siTrack->dz(  firstGoodVertex->position() )           : -9999.        );
        }
        else {
            muon_result.addUserFloat("dxyPV", -999. );
            muon_result.addUserFloat("dzPV", -999. );
        }

        muon_result.addUserFloat("dz_firstPV",siTrack.isNonnull() ? siTrack->dz((vertexCollection->begin())->position()) : -999. );
        muon_result.addUserFloat("dxy_firstPV",siTrack.isNonnull() ? siTrack->dxy((vertexCollection->begin())->position()) : -999. );

        ////////
        // PF //
        ////////

        // PF Isolation
        MuonPFIsolation pfStructR03 = muon->pfIsolationR03();
        MuonPFIsolation pfStructR04 = muon->pfIsolationR04();

        muon_result.addUserFloat("isoR03_pf_ChargedHadronPt", pfStructR03.sumChargedHadronPt              );
        muon_result.addUserFloat("isoR03_pf_ChargedParticlePt", pfStructR03.sumChargedParticlePt            );
        muon_result.addUserFloat("isoR03_pf_NeutralHadronEt", pfStructR03.sumNeutralHadronEt              );
        muon_result.addUserFloat("isoR03_pf_PhotonEt", pfStructR03.sumPhotonEt                     );
        muon_result.addUserFloat("isoR03_pf_sumNeutralHadronEtHighThreshold", pfStructR03.sumNeutralHadronEtHighThreshold );
        muon_result.addUserFloat("isoR03_pf_sumPhotonEtHighThreshold", pfStructR03.sumPhotonEtHighThreshold        );
        muon_result.addUserFloat("isoR03_pf_PUPt", pfStructR03.sumPUPt                         );

        muon_result.addUserFloat("isoR04_pf_ChargedHadronPt", pfStructR04.sumChargedHadronPt              );
        muon_result.addUserFloat("isoR04_pf_ChargedParticlePt", pfStructR04.sumChargedParticlePt            );
        muon_result.addUserFloat("isoR04_pf_NeutralHadronEt", pfStructR04.sumNeutralHadronEt              );
        muon_result.addUserFloat("isoR04_pf_PhotonEt", pfStructR04.sumPhotonEt                     );
        muon_result.addUserFloat("isoR04_pf_sumNeutralHadronEtHighThreshold", pfStructR04.sumNeutralHadronEtHighThreshold );
        muon_result.addUserFloat("isoR04_pf_sumPhotonEtHighThreshold", pfStructR04.sumPhotonEtHighThreshold        );
        muon_result.addUserFloat("isoR04_pf_PUPt", pfStructR04.sumPUPt                         );

        // Other PF
        reco::CandidatePtr pfCandRef = muon->sourceCandidatePtr(0);

        if(pfCandRef.isNonnull()){
            
            muon_result.addUserFloat("pfpt", pfCandRef->p4().pt());
            muon_result.addUserFloat("pfeta", pfCandRef->p4().eta());
            muon_result.addUserFloat("pfphi", pfCandRef->p4().phi());
            muon_result.addUserFloat("pfmass", pfCandRef->p4().mass());
            muon_result.addUserInt("pfcharge", pfCandRef->charge()                                                     );
            muon_result.addUserInt("pfparticleId", pfCandRef->pdgId()                                                      );
            muon_result.addUserInt("pfidx", pfCandRef.key()                                                         );
        }
        else {
            
            muon_result.addUserFloat("pfpt", 0.0);
            muon_result.addUserFloat("pfeta", 0.0);
            muon_result.addUserFloat("pfphi", 0.0);
            muon_result.addUserFloat("pfmass", 0.0);
            muon_result.addUserInt("pfcharge", -9999);
            muon_result.addUserInt("pfparticleId", -9999);
            muon_result.addUserInt("pfidx", -9999);

        } //

        ///////////
        // IP 3D //
        ///////////

        muon_result.addUserFloat("ip3d", muon->dB(pat::Muon::PV3D) ); 
        muon_result.addUserFloat("ip3derr", muon->edB(pat::Muon::PV3D) );
        muon_result.addUserFloat("ip2d", muon->dB(pat::Muon::PV2D) ); 
        muon_result.addUserFloat("ip2derr", muon->edB(pat::Muon::PV2D) );
 
        //////////////////////
        // genMatch miniAOD //
        //////////////////////
    
        LorentzVector mc_p4(0,0,0,0);	 
        const reco::GenParticle * gen = muon->genParticle();
        if (gen != 0) {
            mc_p4 = gen->p4();
            muon_result.addUserInt("mc_patMatch_id", gen->pdgId()  );
            muon_result.addUserFloat("mc_patMatch_pt", mc_p4.pt()         );
            muon_result.addUserFloat("mc_patMatch_eta", mc_p4.eta()         );
            muon_result.addUserFloat("mc_patMatch_phi", mc_p4.phi()         );
            muon_result.addUserFloat("mc_patMatch_mass", mc_p4.mass()         );
            muon_result.addUserFloat("mc_patMatch_dr", ROOT::Math::VectorUtil::DeltaR(gen->p4(), muon->p4())  );

        }
        else {
            muon_result.addUserInt("mc_patMatch_id", -999);
            muon_result.addUserFloat("mc_patMatch_pt", 0.);
            muon_result.addUserFloat("mc_patMatch_eta", 0.);
            muon_result.addUserFloat("mc_patMatch_phi", 0.);
            muon_result.addUserFloat("mc_patMatch_mass", 0.);
            muon_result.addUserFloat("mc_patMatch_dr", -999.);

        }


        //////////////////////
        // mini-isolation   //
        //////////////////////
    
        auto mu2 = muon->clone();
        auto miniiso = mu2->miniPFIsolation();
        muon_result.addUserFloat("miniIso_uncor",miniiso.chargedHadronIso() + miniiso.neutralHadronIso() + miniiso.photonIso());
        muon_result.addUserFloat("miniIso_ch",miniiso.chargedHadronIso()); 
        muon_result.addUserFloat("miniIso_nh",miniiso.neutralHadronIso()); 
        muon_result.addUserFloat("miniIso_em",miniiso.photonIso()); 
        muon_result.addUserFloat("miniIso_db",miniiso.puChargedHadronIso()); 
        delete mu2;
    
        result->emplace_back(muon_result);

    } // end loop on muons

    iEvent.put(std::move(result));

} //


double MuonMaker::muonIsoValuePF(const Muon& mu, const Vertex& vtx, float coner, float minptn, float dzcut, int filterId){

    float pfciso = 0;
    float pfniso = 0;

    const TrackRef siTrack  = mu.innerTrack();

    float mudz = siTrack.isNonnull() ? siTrack->dz(vtx.position()) : mu.standAloneMuon()->dz(vtx.position());

    for (PFCandidateCollection::const_iterator pf=pfCand_h->begin(); pf<pfCand_h->end(); ++pf){

        float dR = deltaR(pf->eta(), pf->phi(), mu.eta(), mu.phi());
        if (dR>coner) continue;

        int pfid = abs(pf->pdgId());
        if (filterId!=0 && filterId!=pfid) continue;

        float pfpt = pf->pt();
        if (pf->charge()==0) {
            //neutrals
            if (pfpt>minptn) pfniso+=pfpt;
        } else {
            //charged
            //avoid double counting of muon itself
            const TrackRef pfTrack  = pf->trackRef();
            if (siTrack.isNonnull()  && pfTrack.isNonnull() && siTrack.key()==pfTrack.key()) continue;
            //first check electrons with gsf track
            if (abs(pf->pdgId())==11 && pf->gsfTrackRef().isNonnull()) {
                if(fabs(pf->gsfTrackRef()->dz(vtx.position()) - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
                continue;//and avoid double counting
            }
            //then check anything that has a ctf track
            if (pfTrack.isNonnull()) {//charged (with a ctf track)
                if(fabs( pfTrack->dz(vtx.position()) - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
            } 
        }
    } 
    return pfciso+pfniso;
}
void MuonMaker::muIsoCustomCone( edm::View<pat::Muon>::const_iterator& mu, float dr, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso){
    chiso     = 0.;
    nhiso     = 0.;
    emiso     = 0.;
    dbiso     = 0.;
    float deadcone_ch = 0.0001;
    float deadcone_pu = 0.01;
    float deadcone_ph = 0.01;
    float deadcone_nh = 0.01;

  double phi = mu->p4().Phi();
  double eta = mu->p4().Eta();
  double pi = M_PI;

    for( pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {
        float id = pf_it->pdgId();
        if (fabs(id) != 211 && fabs(id) != 130 && fabs(id) != 22) continue;

        double deltaPhi = phi-pf_it->p4().Phi();
        if ( deltaPhi > pi ) deltaPhi -= 2.0*pi;
        else if ( deltaPhi <= -pi ) deltaPhi += 2.0*pi;
        deltaPhi = fabs(deltaPhi);
        if (deltaPhi > dr) continue;
        double deltaEta = fabs(pf_it->p4().Eta()-eta);
        if (deltaEta > dr) continue;
        double thisDR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

        if ( thisDR>dr ) continue;  
        float pt = pf_it->p4().pt();
        if ( fabs(id)==211 ) {
            if (pf_it->fromPV() > 1 && (!useVetoCones || thisDR > deadcone_ch) ) chiso+=pt;
            else if ((pf_it->fromPV() <= 1) && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) dbiso+=pt;
        }
        if ( fabs(id)==130 && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_nh) ) nhiso+=pt;
        if ( fabs(id)==22 && (pt > ptthresh) && (!useVetoCones || thisDR > deadcone_ph) ) emiso+=pt;
    }
    return;
}

void MuonMaker::muMiniIso( edm::View<pat::Muon>::const_iterator& mu, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso){

    float pt = mu->p4().pt();
    float dr = 0.2;
    if (pt>50) dr = 10./pt;
    if (pt>200) dr = 0.05;
    muIsoCustomCone(mu,dr,useVetoCones,ptthresh, chiso, nhiso, emiso, dbiso);
    return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
