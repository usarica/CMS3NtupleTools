#include "common_includes.h"
#include "MELAEvent.h"
#include "PDGHelpers.h"
#include "HiggsComparators.h"
#include "TopComparators.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>
#include "TDirectory.h"


using namespace std;
using namespace IvyStreamHelpers;
using namespace PDGHelpers;


void addLHECandMass(TString inpath, TString outpath, int nevents, std::string strVVMode, int VVDecayMode){
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  MELAEvent::CandidateVVMode VVMode = MELAEvent::getCandidateVVModeFromString(strVVMode);

  BaseTree sample_tree(inpath+"/allevents.root", "cms3ntuple/Events", "", "");

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireCoreGenInfo(false);
  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(true);

  genInfoHandler.bookBranches(&sample_tree);
  genInfoHandler.wrapTree(&sample_tree);

  // Unmute all branches that were silenced before
  sample_tree.unmuteAllBranches();

  TFile* foutput = TFile::Open(outpath+"/allevents_modified.root", "recreate");
  foutput->cd();
  TDirectory* curdir = gDirectory;
  TDirectory* subdir = foutput->mkdir("cms3ntuple");
  subdir->cd();

  float LHECandMass;

  TTree* outputTree = sample_tree.getSelectedTree()->CloneTree(0);
  outputTree->Branch("LHECandMass", &LHECandMass);
  outputTree->SetAutoSave(0);

  if (nevents<0) nevents = sample_tree.getSelectedNEvents();
  for (int ev=0; ev<nevents; ev++){
    HelperFunctions::progressbar(ev, nevents);
    sample_tree.getSelectedEvent(ev);

    genInfoHandler.constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler.getGenInfo();

    unsigned int nup = genInfo->extras.lheparticles_id.size();
    std::vector<MELAParticle*> particleList; particleList.reserve(nup);
    for (unsigned int a=0; a<nup; a++){
      TLorentzVector part_p4(
        genInfo->extras.lheparticles_px.at(a),
        genInfo->extras.lheparticles_py.at(a),
        genInfo->extras.lheparticles_pz.at(a),
        genInfo->extras.lheparticles_E.at(a)
      );
      particleList.push_back(new MELAParticle(genInfo->extras.lheparticles_id.at(a), part_p4));
      particleList.back()->setGenStatus(genInfo->extras.lheparticles_status.at(a));
    }
    for (unsigned int a=0; a<nup; a++){
      if (genInfo->extras.lheparticles_mother0_index.at(a)>=0) particleList.at(a)->addMother(particleList.at(genInfo->extras.lheparticles_mother0_index.at(a)));
      if (genInfo->extras.lheparticles_mother1_index.at(a)>=0) particleList.at(a)->addMother(particleList.at(genInfo->extras.lheparticles_mother1_index.at(a)));
    }

    MELAEvent genEvent;
    std::vector<MELAParticle*> writtenGenCands;
    std::vector<MELAParticle*> writtenGenTopCands;
    {
      for (MELAParticle* genPart:particleList){
        if (isATopQuark(genPart->id)){
          writtenGenTopCands.push_back(genPart);
          if (genPart->genStatus==1) genEvent.addIntermediate(genPart);
        }
        if (isAHiggs(genPart->id)){
          writtenGenCands.push_back(genPart);
          if (VVMode==MELAEvent::UndecayedMode && (genPart->genStatus==1 || genPart->genStatus==2)) genEvent.addIntermediate(genPart);
        }
        if (genPart->genStatus==1){
          if (isALepton(genPart->id)) genEvent.addLepton(genPart);
          else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
          else if (isAPhoton(genPart->id)) genEvent.addPhoton(genPart);
          else if (isAKnownJet(genPart->id) && !isATopQuark(genPart->id)) genEvent.addJet(genPart);
        }
        else if (genPart->genStatus==-1) genEvent.addMother(genPart);
      }
    }
    genEvent.constructTopCandidates();
    {
      std::vector<MELATopCandidate_t*> matchedTops;
      for (auto* writtenGenTopCand:writtenGenTopCands){
        MELATopCandidate_t* tmpCand = TopComparators::matchATopToParticle(genEvent, writtenGenTopCand);
        if (tmpCand) matchedTops.push_back(tmpCand);
      }
      for (MELATopCandidate_t* tmpCand:genEvent.getTopCandidates()){
        if (std::find(matchedTops.begin(), matchedTops.end(), tmpCand)==matchedTops.end()) tmpCand->setSelected(false);
      }
    }
    genEvent.constructVVCandidates(VVMode, VVDecayMode);
    genEvent.addVVCandidateAppendages();

    MELACandidate* genCand=nullptr;
    for (auto* writtenGenCand:writtenGenCands){
      MELACandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(genEvent, writtenGenCand);
      if (tmpCand){
        if (!genCand) genCand = tmpCand;
        else genCand = HiggsComparators::candComparator(genCand, tmpCand, HiggsComparators::BestZ1ThenZ2ScSumPt, VVMode);
      }
    }
    if (!genCand) genCand = HiggsComparators::candidateSelector(genEvent, HiggsComparators::BestZ1ThenZ2ScSumPt, VVMode);

    if (genCand) LHECandMass = genCand->m();
    else LHECandMass = -1;

    outputTree->Fill();
    for (MELAParticle*& tmpPart:particleList){ if (tmpPart) delete tmpPart; }
  }

  subdir->WriteTObject(outputTree);
  subdir->Close();
  foutput->Close();
}
