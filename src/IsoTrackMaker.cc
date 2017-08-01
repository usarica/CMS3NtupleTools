#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "CMS3/NtupleMaker/interface/IsoTrackMaker.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

IsoTrackMaker::IsoTrackMaker(const edm::ParameterSet& iConfig){

    pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));
    isoTracksToken = consumes<pat::IsolatedTrackCollection>(iConfig.getParameter<edm::InputTag>("isoTracksTag"));
    pt_cut_        = iConfig.getParameter<double>("pT_cut");
    pt_cut_noIso_  = iConfig.getParameter<double>("pT_cut_noIso");

    produces<vector<LorentzVector> > ("isotracksp4"         ).setBranchAlias("isotracks_p4"         );
    produces<vector<float> >         ("isotracksdz"         ).setBranchAlias("isotracks_dz"         );
    produces<vector<float> >         ("isotracksdxy"        ).setBranchAlias("isotracks_dxy"         );
    produces<vector<float> >         ("isotracksdzError"    ).setBranchAlias("isotracks_dzError"         );
    produces<vector<float> >         ("isotracksdxyError"   ).setBranchAlias("isotracks_dxyError"         );
    produces<vector<int> >           ("isotrackscharge"     ).setBranchAlias("isotracks_charge"     );
    produces<vector<int> >           ("isotracksparticleId" ).setBranchAlias("isotracks_particleId" );
    produces<vector<int> >           ("isotracksfromPV"     ).setBranchAlias("isotracks_fromPV"	    );
    produces<vector<bool> >          ("isotracksisPFCand"   ).setBranchAlias("isotracks_isPFCand"	    );
    produces<vector<float> >         ("isotrackspfisoch"    ).setBranchAlias("isotracks_pfIso_ch"         );
    produces<vector<float> >         ("isotrackspfisonh"    ).setBranchAlias("isotracks_pfIso_nh"         );
    produces<vector<float> >         ("isotrackspfisoem"    ).setBranchAlias("isotracks_pfIso_em"         );
    produces<vector<float> >         ("isotrackspfisodb"    ).setBranchAlias("isotracks_pfIso_db"         );
    produces<vector<float> >         ("isotracksminiisoch"  ).setBranchAlias("isotracks_miniIso_ch"         );
    produces<vector<float> >         ("isotracksminiisonh"  ).setBranchAlias("isotracks_miniIso_nh"         );
    produces<vector<float> >         ("isotracksminiisoem"  ).setBranchAlias("isotracks_miniIso_em"         );
    produces<vector<float> >         ("isotracksminiisodb"  ).setBranchAlias("isotracks_miniIso_db"         );
    produces<vector<bool> >          ("isotracksisHighPurityTrack").setBranchAlias("isotracks_isHighPurityTrack");
    produces<vector<bool> >          ("isotracksisTightTrack").setBranchAlias("isotracks_isTightTrack");
    produces<vector<float> >         ("isotracksmatchedcalojetemenergy").setBranchAlias("isotracks_matchedCaloJetEmEnergy" );
    produces<vector<float> >         ("isotracksmatchedcalojethadenergy").setBranchAlias("isotracks_matchedCaloJetHadEnergy" );
    produces<vector<float> >         ("isotracksdedxstrip").setBranchAlias("isotracks_dEdxStrip" );
    produces<vector<float> >         ("isotracksdedxpixel").setBranchAlias("isotracks_dEdxPixel" );
    produces<vector<float> >         ("isotracksdeltaeta").setBranchAlias("isotracks_deltaEta" );
    produces<vector<float> >         ("isotracksdeltaphi").setBranchAlias("isotracks_deltaPhi" );
    produces<vector<vector<uint16_t> > >("isotrackscrossedecalstatus").setBranchAlias("isotracks_crossedEcalStatus" );
    produces<vector<vector<uint32_t> > >("isotrackscrossedhcalstatus").setBranchAlias("isotracks_crossedHcalStatus" );
    produces<vector<int> >           ("isotrackstrackerLayersWithMeasurement" ).setBranchAlias("isotracks_trackerLayersWithMeasurement"  );
    produces<vector<int> >           ("isotrackspixelLayersWithMeasurement" ).setBranchAlias("isotracks_pixelLayersWithMeasurement"  );
    produces<vector<int> >           ("isotracksnumberOfValidPixelHits" ).setBranchAlias("isotracks_numberOfValidPixelHits"  );
    produces<vector<int> >           ("isotracksnumberOfLostPixelHitsInner" ).setBranchAlias("isotracks_numberOfLostPixelHitsInner"  );
    produces<vector<int> >           ("isotracksnumberOfLostHitsInner" ).setBranchAlias("isotracks_numberOfLostHitsInner"  );
    produces<vector<int> >           ("isotracksnumberOfLostHitsOuter" ).setBranchAlias("isotracks_numberOfLostHitsOuter"  );


}

IsoTrackMaker::~IsoTrackMaker(){}
void IsoTrackMaker::beginRun(const edm::Run&, const edm::EventSetup& es){}
void IsoTrackMaker::beginJob() {}
void IsoTrackMaker::endJob()   {}

void IsoTrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    unique_ptr<vector<LorentzVector> > isotracks_p4                  (new vector<LorentzVector> );
    unique_ptr<vector<float> >         isotracks_dz                  (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_dxy                 (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_dzError             (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_dxyError            (new vector<float>         );
    unique_ptr<vector<int> >           isotracks_charge              (new vector<int>           );
    unique_ptr<vector<int> >           isotracks_particleId          (new vector<int>           );
    unique_ptr<vector<int> >           isotracks_fromPV              (new vector<int>           );
    unique_ptr<vector<bool> >          isotracks_isPFCand            (new vector<bool>          );
    unique_ptr<vector<float> >         isotracks_pfIso_ch            (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_pfIso_nh            (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_pfIso_em            (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_pfIso_db            (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_miniIso_ch          (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_miniIso_nh          (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_miniIso_em          (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_miniIso_db          (new vector<float>         );
    unique_ptr<vector<bool> >          isotracks_isHighPurityTrack   (new vector<bool>          );
    unique_ptr<vector<bool> >          isotracks_isTightTrack        (new vector<bool>          );
    unique_ptr<vector<float> >         isotracks_matchedCaloJetEmEnergy    (new vector<float>   );
    unique_ptr<vector<float> >         isotracks_matchedCaloJetHadEnergy   (new vector<float>   );
    unique_ptr<vector<float> >         isotracks_dEdxStrip           (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_dEdxPixel           (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_deltaEta            (new vector<float>         );
    unique_ptr<vector<float> >         isotracks_deltaPhi            (new vector<float>         );
    unique_ptr<vector<vector<uint16_t> > >  isotracks_crossedEcalStatus  (new vector<vector<uint16_t> >   );
    unique_ptr<vector<vector<uint32_t> > >  isotracks_crossedHcalStatus  (new vector<vector<uint32_t> >   );
    unique_ptr<vector<int> >  isotracks_trackerLayersWithMeasurement (new vector<int>   );
    unique_ptr<vector<int> >  isotracks_pixelLayersWithMeasurement   (new vector<int>   );
    unique_ptr<vector<int> >  isotracks_numberOfValidPixelHits       (new vector<int>   );
    unique_ptr<vector<int> >  isotracks_numberOfLostPixelHitsInner   (new vector<int>   );
    unique_ptr<vector<int> >  isotracks_numberOfLostHitsInner        (new vector<int>   );
    unique_ptr<vector<int> >  isotracks_numberOfLostHitsOuter        (new vector<int>   );

    //get pfcandidates
    Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
    iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);

    Handle<pat::IsolatedTrackCollection> isoTracksHandle;
    iEvent.getByToken(isoTracksToken, isoTracksHandle);
    const pat::IsolatedTrackCollection *isoTracks = isoTracksHandle.product();

    // sort by isotrack pt
    vector<std::pair<int, float> > pts;
    for( pat::IsolatedTrackCollection::const_iterator it = isoTracks->begin(); it != isoTracks->end(); it++ ) {
        pts.push_back(std::pair<int, float>((int)(it-isoTracks->begin()), it->pt()));
    }
    sort(pts.begin(), pts.end(), [](const std::pair<int, float>& v1, const std::pair<int, float>& v2) { return v1.second > v2.second; } );

    for(unsigned int iit=0; iit<pts.size(); iit++){
        const pat::IsolatedTrack& it = isoTracks->at(pts[iit].first);

        float pt = it.pt();

        // use same definition of isolated as miniAOD
        bool isIsolated = it.pfIsolationDR03().chargedHadronIso() < 5.0 ||
            it.pfIsolationDR03().chargedHadronIso() / pt < 0.2 ||
            it.miniPFIsolation().chargedHadronIso() / pt < 0.2;
        if(pt < pt_cut_ || (!isIsolated && pt < pt_cut_noIso_))
            continue;

        bool isPFCand = it.packedCandRef().isNonnull() && it.packedCandRef().id()==pfCandidatesHandle.id();

        isotracks_p4          ->push_back( LorentzVector( it.p4()) );
        isotracks_charge      ->push_back( it.charge() );
        isotracks_particleId  ->push_back( it.pdgId() );
        isotracks_fromPV      ->push_back( it.fromPV() );
        isotracks_isPFCand    ->push_back( isPFCand );
        isotracks_dz          ->push_back( it.dz() );
        isotracks_dxy         ->push_back( it.dxy() );
        isotracks_dzError     ->push_back( it.dzError() );
        isotracks_dxyError    ->push_back( it.dxyError() );
        isotracks_pfIso_ch    ->push_back( it.pfIsolationDR03().chargedHadronIso() );
        isotracks_pfIso_nh    ->push_back( it.pfIsolationDR03().neutralHadronIso() );
        isotracks_pfIso_em    ->push_back( it.pfIsolationDR03().photonIso() );
        isotracks_pfIso_db    ->push_back( it.pfIsolationDR03().puChargedHadronIso() );
        isotracks_miniIso_ch  ->push_back( it.miniPFIsolation().chargedHadronIso() );
        isotracks_miniIso_nh  ->push_back( it.miniPFIsolation().neutralHadronIso() );
        isotracks_miniIso_em  ->push_back( it.miniPFIsolation().photonIso() );
        isotracks_miniIso_db  ->push_back( it.miniPFIsolation().puChargedHadronIso() );
        isotracks_isHighPurityTrack ->push_back(it.isHighPurityTrack() );
        isotracks_isTightTrack      ->push_back(it.isTightTrack() );
        isotracks_matchedCaloJetEmEnergy  ->push_back(it.matchedCaloJetEmEnergy() );
        isotracks_matchedCaloJetHadEnergy ->push_back(it.matchedCaloJetHadEnergy() );
        isotracks_dEdxStrip   ->push_back(it.dEdxStrip() );
        isotracks_dEdxPixel   ->push_back(it.dEdxPixel() );
        isotracks_deltaEta    ->push_back(it.deltaEta() );
        isotracks_deltaPhi    ->push_back(it.deltaPhi() );
        isotracks_crossedEcalStatus ->push_back(it.crossedEcalStatus() );
        isotracks_crossedHcalStatus ->push_back(it.crossedHcalStatus() );
        isotracks_trackerLayersWithMeasurement    ->push_back(it.hitPattern().trackerLayersWithMeasurement() );
        isotracks_pixelLayersWithMeasurement    ->push_back(it.hitPattern().pixelLayersWithMeasurement() );
        isotracks_numberOfValidPixelHits    ->push_back(it.hitPattern().numberOfValidPixelHits() );
        isotracks_numberOfLostPixelHitsInner    ->push_back(it.hitPattern().numberOfLostPixelHits(reco::HitPattern::MISSING_INNER_HITS) );
        isotracks_numberOfLostHitsInner    ->push_back(it.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) );
        isotracks_numberOfLostHitsOuter    ->push_back(it.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) );

    }//loop over candidate collection

    iEvent.put(std::move(isotracks_p4),                   "isotracksp4"       );
    iEvent.put(std::move(isotracks_dz),                   "isotracksdz"       );
    iEvent.put(std::move(isotracks_dxy),                  "isotracksdxy"       );
    iEvent.put(std::move(isotracks_dzError),              "isotracksdzError"   );
    iEvent.put(std::move(isotracks_dxyError),             "isotracksdxyError"  );
    iEvent.put(std::move(isotracks_charge),               "isotrackscharge"    );
    iEvent.put(std::move(isotracks_particleId),           "isotracksparticleId");
    iEvent.put(std::move(isotracks_fromPV),               "isotracksfromPV"    );
    iEvent.put(std::move(isotracks_isPFCand),             "isotracksisPFCand"    );
    iEvent.put(std::move(isotracks_pfIso_ch),             "isotrackspfisoch"    );
    iEvent.put(std::move(isotracks_pfIso_nh),             "isotrackspfisonh"    );
    iEvent.put(std::move(isotracks_pfIso_em),             "isotrackspfisoem"    );
    iEvent.put(std::move(isotracks_pfIso_db),             "isotrackspfisodb"    );
    iEvent.put(std::move(isotracks_miniIso_ch),           "isotracksminiisoch"    );
    iEvent.put(std::move(isotracks_miniIso_nh),           "isotracksminiisonh"    );
    iEvent.put(std::move(isotracks_miniIso_em),           "isotracksminiisoem"    );
    iEvent.put(std::move(isotracks_miniIso_db),           "isotracksminiisodb"    );
    iEvent.put(std::move(isotracks_isHighPurityTrack),    "isotracksisHighPurityTrack" );
    iEvent.put(std::move(isotracks_isTightTrack),         "isotracksisTightTrack"      );
    iEvent.put(std::move(isotracks_matchedCaloJetEmEnergy),  "isotracksmatchedcalojetemenergy"  );
    iEvent.put(std::move(isotracks_matchedCaloJetHadEnergy), "isotracksmatchedcalojethadenergy" );
    iEvent.put(std::move(isotracks_dEdxStrip),            "isotracksdedxstrip"    );
    iEvent.put(std::move(isotracks_dEdxPixel),            "isotracksdedxpixel"    );
    iEvent.put(std::move(isotracks_deltaEta),            "isotracksdeltaeta"    );
    iEvent.put(std::move(isotracks_deltaPhi),            "isotracksdeltaphi"    );
    iEvent.put(std::move(isotracks_crossedEcalStatus),   "isotrackscrossedecalstatus"    );
    iEvent.put(std::move(isotracks_crossedHcalStatus),   "isotrackscrossedhcalstatus"    );
    iEvent.put(std::move(isotracks_trackerLayersWithMeasurement), "isotrackstrackerLayersWithMeasurement"    );
    iEvent.put(std::move(isotracks_pixelLayersWithMeasurement),   "isotrackspixelLayersWithMeasurement"    );
    iEvent.put(std::move(isotracks_numberOfValidPixelHits),       "isotracksnumberOfValidPixelHits"    );
    iEvent.put(std::move(isotracks_numberOfLostPixelHitsInner),   "isotracksnumberOfLostPixelHitsInner"    );
    iEvent.put(std::move(isotracks_numberOfLostHitsInner),        "isotracksnumberOfLostHitsInner"    );
    iEvent.put(std::move(isotracks_numberOfLostHitsOuter),        "isotracksnumberOfLostHitsOuter"    );
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoTrackMaker);
