//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFJetMaker
//
//*\class PFJetMaker PFJetMaker.cc CMS3/NtupleMakerMaker/src/PFJetMaker.cc
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/PFJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig){
    using namespace std;
    using namespace edm;

    pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));
    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    // product of this EDProducer
    produces<vector<LorentzVector> > ( branchprefix+"p4"                               ).setBranchAlias( aliasprefix_+"_p4"                               );
    produces<vector<float> >         ( branchprefix+"undoJEC"                          ).setBranchAlias( aliasprefix_+"_undoJEC"                          );
    produces<vector<float> >         ( branchprefix+"chargedHadronE"                   ).setBranchAlias( aliasprefix_+"_chargedHadronE"                   );
    produces<vector<float> >         ( branchprefix+"neutralHadronE"                   ).setBranchAlias( aliasprefix_+"_neutralHadronE"                   );
    produces<vector<float> >         ( branchprefix+"chargedEmE"                       ).setBranchAlias( aliasprefix_+"_chargedEmE"                       );
    produces<vector<float> >         ( branchprefix+"neutralEmE"                       ).setBranchAlias( aliasprefix_+"_neutralEmE"                       );
    produces<vector<float> >         ( branchprefix+"photonE"                          ).setBranchAlias( aliasprefix_+"_photonE"                          );
    produces<vector<float> >         ( branchprefix+"electronE"                        ).setBranchAlias( aliasprefix_+"_electronE"                        );
    produces<vector<float> >         ( branchprefix+"muonE"                            ).setBranchAlias( aliasprefix_+"_muonE"                            );
    produces<vector<float> >         ( branchprefix+"hfHadronE"                        ).setBranchAlias( aliasprefix_+"_hfHadronE"                        );
    produces<vector<float> >         ( branchprefix+"hfEmE"                            ).setBranchAlias( aliasprefix_+"_hfEmE"                            );
    produces<vector<int> >           ( branchprefix+"chargedHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_chargedHadronMultiplicity"        );
    produces<vector<int> >           ( branchprefix+"neutralHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_neutralHadronMultiplicity"        );
    produces<vector<int> >           ( branchprefix+"photonMultiplicity"               ).setBranchAlias( aliasprefix_+"_photonMultiplicity"               );
    produces<vector<int> >           ( branchprefix+"electronMultiplicity"             ).setBranchAlias( aliasprefix_+"_electronMultiplicity"             );
    produces<vector<int> >           ( branchprefix+"muonMultiplicity"                 ).setBranchAlias( aliasprefix_+"_muonMultiplicity"                 );
    produces<vector<int>   >         ( branchprefix+"chargedMultiplicity"              ).setBranchAlias( aliasprefix_+"_chargedMultiplicity"              );
    produces<vector<int>   >         ( branchprefix+"neutralMultiplicity"              ).setBranchAlias( aliasprefix_+"_neutralMultiplicity"              );
    produces<vector<float> >         ( branchprefix+"area"                             ).setBranchAlias( aliasprefix_+"_area"                             );
    produces<vector<float> >         ( branchprefix+"pileupJetId"                      ).setBranchAlias( aliasprefix_+"_pileupJetId"                      );
    produces<vector<int> >           ( branchprefix+"partonFlavour"                    ).setBranchAlias( aliasprefix_+"_partonFlavour"                    );
    produces<vector<int> >           ( branchprefix+"hadronFlavour"                    ).setBranchAlias( aliasprefix_+"_hadronFlavour"                    );
    produces<vector<vector<LorentzVector> > > ( branchprefix+"pfcandmup4"              ).setBranchAlias( aliasprefix_+"_pfcandmup4"                       );
    produces<vector<int>  >          ( branchprefix+"npfcands"                         ).setBranchAlias( aliasprefix_+"_npfcands"                         );
    produces<vector<int>   >         ( branchprefix+"totalMultiplicity"                ).setBranchAlias( aliasprefix_+"_totalMultiplicity"                );
    produces<vector<float> >         ( branchprefix+"ptDistribution"                   ).setBranchAlias( aliasprefix_+"_ptDistribution"                   );
    produces<vector<float> >         ( branchprefix+"axis1"                            ).setBranchAlias( aliasprefix_+"_axis1"                            );
    produces<vector<float> >         ( branchprefix+"axis2"                            ).setBranchAlias( aliasprefix_+"_axis2"                            );

    // Embedded b-tagging information (miniAOD only)
    produces<vector<float> >         (branchprefix+"pfCombinedInclusiveSecondaryVertexV2BJetTag" ).setBranchAlias(aliasprefix_+"_pfCombinedInclusiveSecondaryVertexV2BJetTag");
    produces<vector<float> >         (branchprefix+"pfDeepCSVJetTagsprobbPlusprobbb"             ).setBranchAlias(aliasprefix_+"_pfDeepCSVJetTagsprobbPlusprobbb"            );
    produces<vector<TString> >       (branchprefix+"bDiscriminatorNames"                         ).setBranchAlias(aliasprefix_+"_bDiscriminatorNames"                     );
    produces<vector<vector<float>> > (branchprefix+"bDiscriminators"                             ).setBranchAlias(aliasprefix_+"_bDiscriminators"                         );

    pfJetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
}

// Destructor
PFJetMaker::~PFJetMaker(){
}

// ------------ method called once each job just before starting event loop  ------------
void PFJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PFJetMaker::endJob() {}

void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

    using namespace std;
    using namespace edm;
    using namespace reco;

    // create containers
    unique_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
    unique_ptr<vector<float> >         pfjets_undoJEC                   (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_chargedHadronE            (new vector<float>          );  
    unique_ptr<vector<float> >         pfjets_neutralHadronE            (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_chargedEmE                (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_neutralEmE                (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_photonE                   (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_electronE                 (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_muonE                     (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_hfHadronE                 (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_hfEmE                     (new vector<float>          );
    unique_ptr<vector<int>   >         pfjets_chargedHadronMultiplicity (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_neutralHadronMultiplicity (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_chargedMultiplicity       (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_neutralMultiplicity       (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_photonMultiplicity        (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_electronMultiplicity      (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_muonMultiplicity          (new vector<int>            );
    unique_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
    unique_ptr<vector<float> >         pfjets_pileupJetId               (new vector<float>          );  
    unique_ptr<vector<int>   >         pfjets_partonFlavour             (new vector<int>            );  
    unique_ptr<vector<int>   >         pfjets_hadronFlavour             (new vector<int>            );  
    unique_ptr<vector<int>   >         pfjets_npfcands                  (new vector<int>            );
    unique_ptr<vector<int>   >         pfjets_totalMultiplicity         (new vector<int>            );
    unique_ptr<vector<float> >         pfjets_ptDistribution            (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_axis1                     (new vector<float>          );
    unique_ptr<vector<float> >         pfjets_axis2                     (new vector<float>          );
    unique_ptr<vector<vector<LorentzVector> > > pfjets_pfcandmup4       (new vector<vector<LorentzVector> > );

    unique_ptr<vector<float> >           pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag (new vector<float>  );
    unique_ptr<vector<float> >           pfjets_pfDeepCSVJetTagsprobbPlusprobbb             (new vector<float>  );
    unique_ptr<vector<TString> >         pfjets_bDiscriminatorNames                         (new vector<TString>        );
    unique_ptr<vector<vector<float> > >  pfjets_bDiscriminators                             (new vector<vector<float> > );
    
    //get pfcandidates
    Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
    iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
    pfCandidates  = pfCandidatesHandle.product();

    //PfJets
    Handle<View<pat::Jet> > pfJetsHandle;
    iEvent.getByToken(pfJetsToken, pfJetsHandle);

    for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

        pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )       );
        // pfjets_mass                      ->push_back( pfjet_it->mass()                      );
        pfjets_undoJEC                   ->push_back( pfjet_it->jecFactor("Uncorrected")    );
        pfjets_chargedHadronE            ->push_back( pfjet_it->chargedHadronEnergy()       );
        pfjets_neutralHadronE            ->push_back( pfjet_it->neutralHadronEnergy()       );
        pfjets_chargedEmE                ->push_back( pfjet_it->chargedEmEnergy()           );
        pfjets_neutralEmE                ->push_back( pfjet_it->neutralEmEnergy()           );
        pfjets_photonE                   ->push_back( pfjet_it->photonEnergy()              );
        pfjets_electronE                 ->push_back( pfjet_it->electronEnergy()            );
        pfjets_muonE                     ->push_back( pfjet_it->muonEnergy()                );
        pfjets_hfHadronE                 ->push_back( pfjet_it->HFHadronEnergy()            );
        pfjets_hfEmE                     ->push_back( pfjet_it->HFEMEnergy()                );
        pfjets_chargedMultiplicity       ->push_back( pfjet_it->chargedMultiplicity()       );
        pfjets_neutralMultiplicity       ->push_back( pfjet_it->neutralMultiplicity()       );
        pfjets_chargedHadronMultiplicity ->push_back( pfjet_it->chargedHadronMultiplicity() );
        pfjets_neutralHadronMultiplicity ->push_back( pfjet_it->neutralHadronMultiplicity() );
        pfjets_photonMultiplicity        ->push_back( pfjet_it->photonMultiplicity()        );
        pfjets_electronMultiplicity      ->push_back( pfjet_it->electronMultiplicity()      );
        pfjets_muonMultiplicity          ->push_back( pfjet_it->muonMultiplicity()          );
        pfjets_area                      ->push_back(pfjet_it->jetArea()                   );
        float pileupJetId = -999; // hedging our beg because this variable isn't yet in the miniAOD we are testing
        if ( pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
        if ( pfjet_it->hasUserFloat("fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
        pfjets_pileupJetId               ->push_back( pileupJetId                          );
        pfjets_partonFlavour             ->push_back(pfjet_it->partonFlavour()             );
        pfjets_hadronFlavour             ->push_back(pfjet_it->hadronFlavour()             );

        //
        // int idx = pfjet_it - pfJetsHandle->begin();

        // // Embedded b-tag info
        // // Default is set automatically to -1000. if no value is found
        const vector<pair<string, float>> bDiscriminatorPairs = pfjet_it->getPairDiscri();
        vector <float> bDiscriminatorPerjet;
        bDiscriminatorPerjet.clear();
        for (size_t bDiscriminator_ind = 0; bDiscriminator_ind < bDiscriminatorPairs.size(); bDiscriminator_ind++ ){
            if (pfjet_it == pfJetsHandle->begin()) pfjets_bDiscriminatorNames->push_back( bDiscriminatorPairs.at(bDiscriminator_ind).first );
            bDiscriminatorPerjet.push_back(pfjet_it->bDiscriminator(string(pfjets_bDiscriminatorNames->at(bDiscriminator_ind))));
        }
        pfjets_bDiscriminators->push_back(bDiscriminatorPerjet);
        pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag->push_back( pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
        pfjets_pfDeepCSVJetTagsprobbPlusprobbb->push_back( pfjet_it->bDiscriminator("pfDeepCSVJetTags:probb") + pfjet_it->bDiscriminator("pfDeepCSVJetTags:probbb") );

        std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 
        vector<LorentzVector> pfcandmup4;
        for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){
            unsigned int ipf = cand_it->key();
            pat::PackedCandidate pfc = pfCandidates->at(ipf);
            if (!pfc.isGlobalMuon() && !pfc.isStandAloneMuon()) continue;
            // LorentzVector pfmup4 = LorentzVector(pfc.p4());
            pfcandmup4.push_back(LorentzVector(pfc.p4()));
            // idx aliasprefix_ ipf pfmup4.Pt() pfc.pdgId() pfc.isGlobalMuon() pfc.isStandAloneMuon()
        } 
        pfjets_pfcandmup4->push_back( pfcandmup4 );
        pfjets_npfcands->push_back(pfjet_cands.size());

        // Do calculation of top-tagger variables
        const bool computeTopTaggerVariables = true;
        if (computeTopTaggerVariables) {
            int totalMult = 0;
            float ptD     = 0;
            float axis1   = 0;
            float axis2   = 0;
            if (pfjet_it->numberOfDaughters() != 0) {
                float sum_weight(0.0), sum_dEta(0.0), sum_dPhi(0.0), sum_dEta2(0.0), sum_dPhi2(0.0), sum_dEta_dPhi(0.0), sum_pt(0.0);

                // loop over the jet constituents (packed candidate situation)
                for (auto part : pfjet_it->getJetConstituentsQuick()) {
                    if (part->charge()) { // charged particles
                        auto p = dynamic_cast<const pat::PackedCandidate*>(part);
                        if (!p) std::cout << "ERROR: QGTagging variables cannot be computed for these jets!" << std::endl;
                        if (!(p->fromPV() > 1 && p->trackHighPurity())) continue;
                        ++totalMult;
                    } else { // neutral particles
                        if (part->pt() < 1.0) continue;
                        ++totalMult;
                    } // charged, neutral particles

                    float dEta   = part->eta() - pfjet_it->eta();
                    float dPhi   = reco::deltaPhi(part->phi(), pfjet_it->phi());
                    float partPt = part->pt();
                    float weight = partPt*partPt;

                    sum_weight    += weight;
                    sum_pt        += partPt;
                    sum_dEta      += dEta      * weight;
                    sum_dPhi      += dPhi      * weight;
                    sum_dEta2     += dEta*dEta * weight;
                    sum_dEta_dPhi += dEta*dPhi * weight;
                    sum_dPhi2     += dPhi*dPhi * weight;
                } // pfjet_it->getJetConstituentsQuick()

                // calculate axis2 and ptD
                if (sum_weight > 0) {
                    ptD = sqrt(sum_weight)/sum_pt;
                    float ave_dEta  = sum_dEta  / sum_weight;
                    float ave_dPhi  = sum_dPhi  / sum_weight;
                    float ave_dEta2 = sum_dEta2 / sum_weight;
                    float ave_dPhi2 = sum_dPhi2 / sum_weight;
                    float a = ave_dEta2 - ave_dEta*ave_dEta;
                    float b = ave_dPhi2 - ave_dPhi*ave_dPhi;
                    float c = -(sum_dEta_dPhi/sum_weight - ave_dEta*ave_dPhi);
                    float delta = sqrt(fabs( (a-b)*(a-b) + 4*c*c ));
                    if(a+b-delta > 0) axis2 = sqrt(0.5*(a+b-delta));
                    else              axis2 = 0.0;
                    if(a+b+delta > 0) axis1 = sqrt(0.5*(a+b+delta));
                    else              axis1 = 0.0;
                }
            }

            pfjets_totalMultiplicity->push_back(totalMult);
            pfjets_ptDistribution   ->push_back(ptD);
            pfjets_axis1            ->push_back(axis1);
            pfjets_axis2            ->push_back(axis2);
        } else {
            pfjets_ptDistribution   ->push_back(pfjet_it->constituentPtDistribution());
            pfjets_totalMultiplicity->push_back(pfjet_it->numberOfDaughters());
            pfjets_axis1            ->push_back(-1);
            pfjets_axis2            ->push_back(-1);
        }

    }

    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(std::move(pfjets_p4                        ), branchprefix+"p4"                        );
    iEvent.put(std::move(pfjets_undoJEC                   ), branchprefix+"undoJEC"                   );
    iEvent.put(std::move(pfjets_chargedHadronE            ), branchprefix+"chargedHadronE"            );
    iEvent.put(std::move(pfjets_neutralHadronE            ), branchprefix+"neutralHadronE"            );
    iEvent.put(std::move(pfjets_chargedEmE                ), branchprefix+"chargedEmE"                );
    iEvent.put(std::move(pfjets_neutralEmE                ), branchprefix+"neutralEmE"                );
    iEvent.put(std::move(pfjets_photonE                   ), branchprefix+"photonE"                   );
    iEvent.put(std::move(pfjets_electronE                 ), branchprefix+"electronE"                 );
    iEvent.put(std::move(pfjets_muonE                     ), branchprefix+"muonE"                     );
    iEvent.put(std::move(pfjets_hfHadronE                 ), branchprefix+"hfHadronE"                 );
    iEvent.put(std::move(pfjets_hfEmE                     ), branchprefix+"hfEmE"                     );  
    iEvent.put(std::move(pfjets_chargedMultiplicity       ), branchprefix+"chargedMultiplicity"       );
    iEvent.put(std::move(pfjets_neutralMultiplicity       ), branchprefix+"neutralMultiplicity"       );
    iEvent.put(std::move(pfjets_chargedHadronMultiplicity ), branchprefix+"chargedHadronMultiplicity" );
    iEvent.put(std::move(pfjets_neutralHadronMultiplicity ), branchprefix+"neutralHadronMultiplicity" );
    iEvent.put(std::move(pfjets_photonMultiplicity        ), branchprefix+"photonMultiplicity"        );
    iEvent.put(std::move(pfjets_electronMultiplicity      ), branchprefix+"electronMultiplicity"      );
    iEvent.put(std::move(pfjets_muonMultiplicity          ), branchprefix+"muonMultiplicity"          );
    iEvent.put(std::move(pfjets_area                      ), branchprefix+"area"                      );
    iEvent.put(std::move(pfjets_pileupJetId               ), branchprefix+"pileupJetId"               );
    iEvent.put(std::move(pfjets_partonFlavour             ), branchprefix+"partonFlavour"             );
    iEvent.put(std::move(pfjets_hadronFlavour             ), branchprefix+"hadronFlavour"             );
    iEvent.put(std::move(pfjets_pfcandmup4                ), branchprefix+"pfcandmup4"                );
    iEvent.put(std::move(pfjets_npfcands                  ), branchprefix+"npfcands"                  );
    iEvent.put(std::move(pfjets_totalMultiplicity         ), branchprefix+"totalMultiplicity"         );
    iEvent.put(std::move(pfjets_ptDistribution            ), branchprefix+"ptDistribution"            );
    iEvent.put(std::move(pfjets_axis1                     ), branchprefix+"axis1"                     );
    iEvent.put(std::move(pfjets_axis2                     ), branchprefix+"axis2"                     );

    iEvent.put(std::move(pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag), branchprefix+"pfCombinedInclusiveSecondaryVertexV2BJetTag");  
    iEvent.put(std::move(pfjets_pfDeepCSVJetTagsprobbPlusprobbb            ), branchprefix+"pfDeepCSVJetTagsprobbPlusprobbb"            );
    iEvent.put(std::move(pfjets_bDiscriminatorNames                        ), branchprefix+"bDiscriminatorNames"     );
    iEvent.put(std::move(pfjets_bDiscriminators                            ), branchprefix+"bDiscriminators"         );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
