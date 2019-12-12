#include <MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h>
#include <CMS3/NtupleMaker/interface/CMS3MELAHelpers.h>


using namespace std;


namespace CMS3MELAHelpers{
  std::shared_ptr<Mela> melaHandle(nullptr);

  int getSqrts(int year){
    switch (year){
    case 2011:
      return 7;
    case 2012:
      return 8;
    case 2015:
    case 2016:
    case 2017:
    case 2018:
      return 13;
    default:
      return -1;
    }
  }
  void setupMela(int year, float mh, TVar::VerbosityLevel verbosity){
    int sqrts = getSqrts(year);
    if (!melaHandle && sqrts>0 && mh>=0.f) melaHandle = make_shared<Mela>(sqrts, mh, verbosity);
  }
  void clearMela(){ melaHandle.reset((Mela*) nullptr); }


  GMECBlock::GMECBlock(){}
  GMECBlock::~GMECBlock(){ clearMELABranches(); }

  void GMECBlock::setRefTrees(std::vector<TTree*> const& reftrees_){ reftrees = reftrees_; }

  void GMECBlock::buildMELABranches(std::vector<std::string> const& MElist, bool isGen){
    /***************************/
    /***** LHE ME BRANCHES *****/
    /***************************/
    if (isGen){
      // Loop over the ME list to book the non-copy branches
      for (auto const& strLHEME:MElist){
        MELAOptionParser* me_opt;
        // First find out if the option has a copy specification
        // These copy options will be evaulated in a separate loop
        if (strLHEME.find("Copy")!=string::npos){
          me_opt = new MELAOptionParser(strLHEME);
          lheme_copyopts.push_back(me_opt);
          continue;
        }

        // Create a hypothesis for each option
        MELAHypothesis* me_hypo = new MELAHypothesis(melaHandle.get(), strLHEME);
        lheme_units.push_back(me_hypo);

        // Check if the option is aliased
        me_opt = me_hypo->getOption();
        if (me_opt->isAliased()) lheme_aliased_units.push_back(me_hypo);

        // Create a computation for each hypothesis
        MELAComputation* me_computer = new MELAComputation(me_hypo);
        lheme_computers.push_back(me_computer);

        // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
        GMECHelperFunctions::addToMELACluster(me_computer, lheme_clusters);

        this->bookMELABranches(me_opt, me_computer, false);
      }
      // Resolve copy options
      for (MELAOptionParser* me_opt:lheme_copyopts){
        MELAHypothesis* original_hypo=nullptr;
        MELAOptionParser* original_opt=nullptr;

        // Find the original options
        for (auto* me_aliased_unit:lheme_aliased_units){
          if (me_opt->testCopyAlias(me_aliased_unit->getOption()->getAlias())){
            original_hypo = me_aliased_unit;
            original_opt = original_hypo->getOption();
            break;
          }
        }
        if (!original_opt) continue;
        else me_opt->pickOriginalOptions(original_opt);

        // Create a new computation for the copy options
        MELAComputation* me_computer = new MELAComputation(original_hypo);
        me_computer->setOption(me_opt);
        lheme_computers.push_back(me_computer);

        // The rest is the same story...
        // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
        GMECHelperFunctions::addToMELACluster(me_computer, lheme_clusters);

        // Create the necessary branches for each computation
        // Notice that no tree is passed, so no TBranches are created.
        this->bookMELABranches(me_opt, me_computer, true);
      }
      // Loop over the computations to add any contingencies to aliased hypotheses
      for (auto& me_computer:lheme_computers) me_computer->addContingencies(lheme_aliased_units);

      if (!lheme_clusters.empty()){
        cout << "GMECBlock::buildMELABranches: LHE ME clusters:" << endl;
        for (auto const* me_cluster:lheme_clusters){
          cout << "\t- Cluster " << me_cluster->getName() << " has " << me_cluster->getComputations()->size() << " computations registered." << endl;
        }
      }
    }
    /****************************/
    /***** RECO ME BRANCHES *****/
    /****************************/
    else{
      // Loop over the ME list to book the non-copy branches
      for (auto const& strRecoME:MElist){
        MELAOptionParser* me_opt;
        // First find out if the option has a copy specification
        // These copy options will be evaulated in a separate loop
        if (strRecoME.find("Copy")!=string::npos){
          me_opt = new MELAOptionParser(strRecoME);
          recome_copyopts.push_back(me_opt);
          continue;
        }

        // Create a hypothesis for each option
        MELAHypothesis* me_hypo = new MELAHypothesis(melaHandle.get(), strRecoME);
        recome_units.push_back(me_hypo);

        // Check if the option is aliased
        me_opt = me_hypo->getOption();
        if (me_opt->isAliased()) recome_aliased_units.push_back(me_hypo);

        // Create a computation for each hypothesis
        MELAComputation* me_computer = new MELAComputation(me_hypo);
        recome_computers.push_back(me_computer);

        // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
        GMECHelperFunctions::addToMELACluster(me_computer, recome_clusters);

        this->bookMELABranches(me_opt, me_computer, false);
      }
      // Resolve copy options
      for (MELAOptionParser* me_opt:recome_copyopts){
        MELAHypothesis* original_hypo=nullptr;
        MELAOptionParser* original_opt=nullptr;

        // Find the original options
        for (auto* me_aliased_unit:recome_aliased_units){
          if (me_opt->testCopyAlias(me_aliased_unit->getOption()->getAlias())){
            original_hypo = me_aliased_unit;
            original_opt = original_hypo->getOption();
            break;
          }
        }
        if (!original_opt) continue;
        else me_opt->pickOriginalOptions(original_opt);

        // Create a new computation for the copy options
        MELAComputation* me_computer = new MELAComputation(original_hypo);
        me_computer->setOption(me_opt);
        recome_computers.push_back(me_computer);

        // The rest is the same story...
        // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
        GMECHelperFunctions::addToMELACluster(me_computer, recome_clusters);

        // Create the necessary branches for each computation
        // Notice that no tree is passed, so no TBranches are created.
        this->bookMELABranches(me_opt, me_computer, true);
      }
      // Loop over the computations to add any contingencies to aliased hypotheses
      for (auto& me_computer:recome_computers) me_computer->addContingencies(recome_aliased_units);

      if (!recome_clusters.empty()){
        cout << "GMECBlock::buildMELABranches: Reco ME clusters:" << endl;
        for (auto const* me_cluster:recome_clusters){
          cout << "\t- Cluster " << me_cluster->getName() << " has " << me_cluster->getComputations()->size() << " computations registered." << endl;
        }
      }
    }

  }
  void GMECBlock::bookMELABranches(MELAOptionParser* me_opt, MELAComputation* computer, bool doCopy){
    if (!me_opt){
      cerr << "GMECBlock::bookMELABranches: Did not receive a valid me_opt. Something went wrong." << endl;
      throw std::exception();
    }

    std::vector<MELABranch*>* me_branches = (me_opt->isGen() ? &(this->lheme_branches) : &(this->recome_branches));

    std::vector<TTree*> trees;
    if (reftrees.empty() || doCopy) trees.push_back(nullptr);
    else trees = reftrees;

    if (me_opt->doBranch()){
      for (auto* tree:trees){
        string basename = me_opt->getName();
        if (me_opt->isGen()) basename = string("Gen_") + basename;
        MELABranch* tmpbranch;
        Float_t defVal=1.;
        if (me_opt->hasPAux()){
          tmpbranch = new MELABranch(
            tree, TString((string("pAux_") + basename).c_str()),
            defVal, computer
          );
          me_branches->push_back(tmpbranch);
        }
        if (me_opt->hasPConst()){
          tmpbranch = new MELABranch(
            tree, TString((string("pConst_") + basename).c_str()),
            defVal, computer
          );
          me_branches->push_back(tmpbranch);
        }
        defVal = me_opt->getDefaultME();
        tmpbranch = new MELABranch(
          tree, TString((string("p_") + basename).c_str()),
          defVal, computer
        );
        me_branches->push_back(tmpbranch);
        cout << "GMECBlock::bookMELABranches: Constructed branch with base name " << basename << endl;
      }
    }
  }

  void GMECBlock::clearMELABranches(bool isGen){
#define CLEAR_MELA_BRANCHES_CMD(thelist) \
for (auto*& v:thelist) delete v; \
thelist.clear();

    if (isGen){
      CLEAR_MELA_BRANCHES_CMD(lheme_branches);
      CLEAR_MELA_BRANCHES_CMD(lheme_clusters);
      CLEAR_MELA_BRANCHES_CMD(lheme_computers);
      CLEAR_MELA_BRANCHES_CMD(lheme_copyopts);
      CLEAR_MELA_BRANCHES_CMD(lheme_originalopts);
      // Do not delete me_aliased_units. They are deleted together with me_units.
      CLEAR_MELA_BRANCHES_CMD(lheme_units);
    }
    else{
      CLEAR_MELA_BRANCHES_CMD(recome_branches);
      CLEAR_MELA_BRANCHES_CMD(recome_clusters);
      CLEAR_MELA_BRANCHES_CMD(recome_computers);
      CLEAR_MELA_BRANCHES_CMD(recome_copyopts);
      CLEAR_MELA_BRANCHES_CMD(recome_originalopts);
      // Do not delete me_aliased_units. They are deleted together with me_units.
      CLEAR_MELA_BRANCHES_CMD(recome_units);
    }

#undef CLEAR_MELA_BRANCHES_CMD
  }
  void GMECBlock::clearMELABranches(){
    clearMELABranches(true);
    clearMELABranches(false);
  }

  void GMECBlock::computeMELABranches(bool isGen){
    if (!melaHandle) return;
    // Sequantial computation
    updateMELAClusters_Common("Common", isGen);
    updateMELAClusters_J1JECJER("J1JECNominal", isGen); updateMELAClusters_J1JECJER("J1JECUp", isGen); updateMELAClusters_J1JECJER("J1JECDn", isGen); updateMELAClusters_J1JECJER("J1JERUp", isGen); updateMELAClusters_J1JECJER("J1JERDn", isGen);
    updateMELAClusters_J2JECJER("J2JECNominal", isGen); updateMELAClusters_J2JECJER("J2JECUp", isGen); updateMELAClusters_J2JECJER("J2JECDn", isGen); updateMELAClusters_J2JECJER("J2JERUp", isGen); updateMELAClusters_J2JECJER("J2JERDn", isGen);
    updateMELAClusters_LepWH("LepWH", isGen);
    updateMELAClusters_LepZH("LepZH", isGen);
    updateMELAClusters_NoInitialQ("NoInitialQ", isGen);
    updateMELAClusters_NoInitialG("NoInitialG", isGen);
    updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZ", isGen);
    updateMELAClusters_BestLOAssociatedW("BestLOAssociatedW", isGen);
    updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBF", isGen);
    updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximation", isGen);
    updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximation", isGen);
    updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximation", isGen);
    updateMELAClusters_NoAssociatedG("NoAssociatedG", isGen);
    updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedG", isGen);
    // Reverse sequence
    updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedGLast", isGen);
    updateMELAClusters_NoAssociatedG("NoAssociatedGLast", isGen);
    updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximationLast", isGen);
    updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximationLast", isGen);
    updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximationLast", isGen);
    updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBFLast", isGen);
    updateMELAClusters_BestLOAssociatedW("BestLOAssociatedWLast", isGen);
    updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZLast", isGen);
    updateMELAClusters_NoInitialG("NoInitialGLast", isGen);
    updateMELAClusters_NoInitialQ("NoInitialQLast", isGen);
    updateMELAClusters_LepZH("LepZHLast", isGen);
    updateMELAClusters_LepWH("LepWHLast", isGen);
    updateMELAClusters_J2JECJER("J2JECNominalLast", isGen); updateMELAClusters_J2JECJER("J2JECUpLast", isGen); updateMELAClusters_J2JECJER("J2JECDnLast", isGen);
    updateMELAClusters_J1JECJER("J1JECNominalLast", isGen); updateMELAClusters_J1JECJER("J1JECUpLast", isGen); updateMELAClusters_J1JECJER("J1JECDnLast", isGen);
    updateMELAClusters_Common("CommonLast", isGen);
  }
  void GMECBlock::computeMELABranches(){
    computeMELABranches(true);
    computeMELABranches(false);
  }

  void GMECBlock::pushMELABranches(bool isGen){
    std::vector<MELABranch*>& me_branches = (isGen ? lheme_branches : recome_branches);
    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    // Pull + push...
    for (MELABranch* v:me_branches) v->setVal();
    // ...then reset
    for (MELACluster* v:me_clusters) v->reset();
  }
  void GMECBlock::pushMELABranches(){
    pushMELABranches(true);
    pushMELABranches(false);
  }

  void GMECBlock::getBranchValues(std::unordered_map<std::string, float>& io_rcd, bool isGen){
    std::vector<MELABranch*> const& me_branches = (isGen ? lheme_branches : recome_branches);
    for (MELABranch* const& br:me_branches) io_rcd[br->bname.Data()] = br->getVal();
  }
  void GMECBlock::getBranchValues(std::unordered_map<std::string, float>& io_rcd){
    getBranchValues(io_rcd, true);
    getBranchValues(io_rcd, false);
  }

  // Common ME computations that do not manipulate the LHE candidate
  void GMECBlock::updateMELAClusters_Common(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }
  }
  // Common ME computations for JECNominal, Up and Down variations, case where ME requires 1 jet
  void GMECBlock::updateMELAClusters_J1JECJER(const string clustertype, bool isGen){
    constexpr int njecnum = 3;
    const int jecnumsel = -1
      + int(clustertype=="J1JECNominal")*1
      + int(clustertype=="J1JECUp")*2
      + int(clustertype=="J1JECDn")*3;
    if (jecnumsel<0) return;

    //cout << "Begin GMECBlock::updateMELAClusters_J1JECJER(" << clustertype << "," << isGen << ")" << endl;

    // First determine if any of the candidates has only one jet
    bool doSkip=true;
    int nMelaStored = melaHandle->getNCandidates();
    bool setJECcand = (nMelaStored == njecnum);
    if (!(setJECcand || nMelaStored==0)) return;

    for (int jecnum=0; jecnum<njecnum; jecnum++){
      if (setJECcand) melaHandle->setCurrentCandidateFromIndex(jecnum);
      else if (jecnum>0) break;

      MELACandidate* melaCand = melaHandle->getCurrentCandidate();
      if (!melaCand) continue;

      unsigned int nGoodJets=melaCand->getNAssociatedJets();
      doSkip = doSkip && (nGoodJets!=1);
    }
    if (doSkip) return; // If none of the candidates have exactly 1 jet, skip the computations

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

    for (int jecnum=0; jecnum<njecnum; jecnum++){
      if (jecnum!=jecnumsel) continue;

      if (setJECcand) melaHandle->setCurrentCandidateFromIndex(jecnum);
      else if (jecnum>0) break;

      MELACandidate* melaCand = melaHandle->getCurrentCandidate();
      if (!melaCand) continue;

      const unsigned int nGoodJets=std::min(1, melaCand->getNAssociatedJets());
      for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet

        for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(disableJet==firstjet); // Disable the other jets

        for (MELACluster* theCluster:me_clusters){
          if (theCluster->getName()==clustertype){
            // Re-compute all related hypotheses first...
            theCluster->computeAll();
            // ...then force an update the cluster
            theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
          }
        } // End loop over clusters

      } // End loop over first jet
        // Turn associated jets back on
      for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
    } // End jecnum loop

      //cout << "End GMECBlock::updateMELAClusters_J1JECJER(" << clustertype << "," << isGen << ")" << endl;
  }
  // Common ME computations for JECNominal, Up and Down variations, case where ME requires 2 jets
  void GMECBlock::updateMELAClusters_J2JECJER(const string clustertype, bool isGen){
    constexpr int njecnum = 3;
    const int jecnumsel = -1
      + int(clustertype=="J2JECNominal")*1
      + int(clustertype=="J2JECUp")*2
      + int(clustertype=="J2JECDn")*3;
    if (jecnumsel<0) return;

    //cout << "Begin GMECBlock::updateMELAClusters_J2JECJER(" << clustertype << "," << isGen << ")" << endl;

    int nMelaStored = melaHandle->getNCandidates();
    bool setJECcand = (nMelaStored == njecnum);
    if (!(setJECcand || nMelaStored==0)) return;

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

    for (int jecnum=0; jecnum<njecnum; jecnum++){
      if (jecnum!=jecnumsel) continue;

      if (setJECcand) melaHandle->setCurrentCandidateFromIndex(jecnum);
      else if (jecnum>0) break;

      MELACandidate* melaCand = melaHandle->getCurrentCandidate();
      if (!melaCand) continue;

      const unsigned int nGoodJets=melaCand->getNAssociatedJets();
      for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet
        for (unsigned int secondjet = firstjet+1; secondjet < nGoodJets; secondjet++){ // Loop over second jet
                                                                                       // Disable jets and tops
          for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected((disableJet==firstjet || disableJet==secondjet)); // Disable the other jets
          unsigned int nDisabledStableTops=0;
          for (MELATopCandidate_t* einTop:melaCand->getAssociatedTops()){
            if (einTop->getNDaughters()==3) einTop->setSelected(false); // All unstable tops are disabled in the loop for jets (where "jet"=="stable top") since we are looping over jecnum
            else{
              einTop->setSelected((nDisabledStableTops==firstjet || nDisabledStableTops==secondjet)); // Disable the other stable tops
              nDisabledStableTops++;
            }
          }

          for (MELACluster* theCluster:me_clusters){
            if (theCluster->getName()==clustertype){
              // Re-compute all related hypotheses first...
              theCluster->computeAll();
              // ...then force an update the cluster
              theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
            }
          } // End loop over clusters

        } // End loop over second jet
      } // End loop over first jet
        // Turn associated jets/tops back on
      for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
      for (MELATopCandidate_t* einTop:melaCand->getAssociatedTops()) einTop->setSelected(true); // Turn all tops back on
    } // End jecnum loop

      //cout << "End GMECBlock::updateMELAClusters_J2JECJER(" << clustertype << "," << isGen << ")" << endl;
  }
  // Common ME computations for leptonic WH: Loops over possible fake neutrinos
  void GMECBlock::updateMELAClusters_LepWH(const string clustertype, bool isGen){
    int nMelaStored = melaHandle->getNCandidates();
    bool setLepHypoCand = (nMelaStored == 1);
    if (!(setLepHypoCand || nMelaStored==0)) return;

    if (setLepHypoCand) melaHandle->setCurrentCandidateFromIndex(0);
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

    for (MELAParticle* nu:melaCand->getAssociatedNeutrinos()){
      // Notice: Looping over Ws does not make much sense unless you have more than one lepton since the fake neutrino is already calculated from the available lepton with W mass constraint.
      // Such a loop over Ws only makes sense if there are more than one lepton in the event, but in that case, it still does not make sense to cross-match neutrinos and leptons.
      for (MELAParticle* dnu:melaCand->getAssociatedNeutrinos()) dnu->setSelected((dnu==nu)); // Disable all neutrinos other than index==inu
      for (MELACluster* theCluster:me_clusters){
        if (theCluster->getName()==clustertype){
          // Re-compute all related hypotheses first...
          theCluster->computeAll();
          // ...then force an update the cluster
          theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
        }
      } // End loop over clusters
    } // End loop over possible neutrinos
      // Re-enable all neutrinos
    for (MELAParticle* dnu:melaCand->getAssociatedNeutrinos()) dnu->setSelected(true);
  }
  // Common ME computations for leptonic ZH: Picks best Z3
  void GMECBlock::updateMELAClusters_LepZH(const string clustertype, bool isGen){
    int nMelaStored = melaHandle->getNCandidates();
    bool setLepHypoCand = (nMelaStored == 1);
    if (!(setLepHypoCand || nMelaStored==0)) return;

    if (setLepHypoCand) melaHandle->setCurrentCandidateFromIndex(0);
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

    std::vector<MELAParticle*> associatedVs;
    for (MELAParticle* associatedV:melaCand->getAssociatedSortedVs()){
      if (!PDGHelpers::isAZBoson(associatedV->id)) continue;
      if (!PDGHelpers::isALepton(associatedV->getDaughter(0)->id)) continue;
      bool passSelection=true;
      for (MELAParticle* dauV:associatedV->getDaughters()) passSelection &= dauV->passSelection;
      if (!passSelection) continue;
      associatedVs.push_back(associatedV);
    }

    double dZmass=-1;
    MELAParticle* chosenZ=nullptr;
    // Choose the Z by mass closest to mZ (~equivalent to ordering by best SM ME but would be equally valid for BSM MEs as well)
    for (MELAParticle* associatedV:associatedVs){
      if (!chosenZ || fabs(associatedV->m()-PDGHelpers::Zmass)<dZmass){ dZmass=fabs(associatedV->m()-PDGHelpers::Zmass); chosenZ=associatedV; }
    }
    if (chosenZ){
      // Disable every associated Z boson (and not its daughters!) unless it is the chosen one
      for (MELAParticle* associatedV:associatedVs) associatedV->setSelected((associatedV==chosenZ));

      for (MELACluster* theCluster:me_clusters){
        if (theCluster->getName()==clustertype){
          // Re-compute all related hypotheses first...
          theCluster->computeAll();
          // ...then force an update the cluster
          theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
        }
      } // End loop over clusters

    } // End if chosenZ>=0
      // Re-enable every associated Z boson and its daughters unless it is the chosen one
    for (MELAParticle* associatedV:associatedVs) associatedV->setSelected(true);
  }
  // ME computations that require no quark initial state
  void GMECBlock::updateMELAClusters_NoInitialQ(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of quark mothers
    std::vector<int> motherIds;
    for (int imot=0; imot<melaCand->getNMothers(); imot++){
      motherIds.push_back(melaCand->getMother(imot)->id);
      if (PDGHelpers::isAQuark(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  }
  // ME computations that require no gluon initial state
  void GMECBlock::updateMELAClusters_NoInitialG(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    std::vector<int> motherIds;
    for (int imot=0; imot<melaCand->getNMothers(); imot++){
      motherIds.push_back(melaCand->getMother(imot)->id);
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  }
  // ME computations that require no gluons as associated particles
  void GMECBlock::updateMELAClusters_NoAssociatedG(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    std::vector<int> ajetIds;
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
      ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  }
  // ME computations that require no gluon initial state and no gluons as associated particles
  void GMECBlock::updateMELAClusters_NoInitialGNoAssociatedG(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    std::vector<int> motherIds;
    std::vector<int> ajetIds;
    for (int imot=0; imot<melaCand->getNMothers(); imot++){
      motherIds.push_back(melaCand->getMother(imot)->id);
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
      ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  }
  // ME computations that require best Z, W or VBF topology at LO (no gluons)
  void GMECBlock::updateMELAClusters_BestLOAssociatedZ(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    std::vector<int> motherIds;
    std::vector<int> ajetIds;
    for (int imot=0; imot<melaCand->getNMothers(); imot++){
      motherIds.push_back(melaCand->getMother(imot)->id);
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
      ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    std::vector<MELAParticle*> associatedVs; // Vector of Zs to loop over
    for (MELAParticle* Vtmp:melaCand->getAssociatedSortedVs()){
      if (Vtmp && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
        bool passSelection=true;
        for (MELAParticle* dauVtmp:Vtmp->getDaughters()) passSelection &= dauVtmp->passSelection;
        if (!passSelection) continue;

        associatedVs.push_back(Vtmp);
      }
    }

    // Give precedence to leptonic V decays
    bool hasALepV=false;
    for (MELAParticle* Vtmp:associatedVs){
      const int& Vtmp_dauid = Vtmp->getDaughter(0)->id;
      if (
        PDGHelpers::isALepton(Vtmp_dauid)
        ||
        PDGHelpers::isANeutrino(Vtmp_dauid)
        ){
        hasALepV=true;
        break;
      }
    }

    std::vector<MELAParticle*> modifiedVs;
    MELAParticle* bestVbyMass=nullptr;
    float bestVMassDiff=-1;
    for (MELAParticle* Vtmp:associatedVs){
      if (
        hasALepV &&
        PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
        ){
        for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(false);
        modifiedVs.push_back(Vtmp);
      }
      else if (!bestVbyMass || fabs(Vtmp->m()-PDGHelpers::Zmass)<bestVMassDiff){
        bestVMassDiff=fabs(Vtmp->m()-PDGHelpers::Zmass);
        bestVbyMass = Vtmp;
      }
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
    for (MELAParticle* Vtmp:modifiedVs){ for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(true); }
  }
  void GMECBlock::updateMELAClusters_BestLOAssociatedW(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    std::vector<int> motherIds;
    std::vector<int> ajetIds;
    for (int imot=0; imot<melaCand->getNMothers(); imot++){
      motherIds.push_back(melaCand->getMother(imot)->id);
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
      ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    std::vector<MELAParticle*> associatedVs; // Vector of Ws to loop over
    for (MELAParticle* Vtmp:melaCand->getAssociatedSortedVs()){
      if (Vtmp && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
        bool passSelection=true;
        for (MELAParticle* dauVtmp:Vtmp->getDaughters()) passSelection &= dauVtmp->passSelection;
        if (!passSelection) continue;

        associatedVs.push_back(Vtmp);
      }
    }
    // Give precedence to leptonic V decays
    bool hasALepV=false;
    for (MELAParticle* Vtmp:associatedVs){
      const int& Vtmp_dauid = Vtmp->getDaughter(0)->id;
      if (
        PDGHelpers::isALepton(Vtmp_dauid)
        ||
        PDGHelpers::isANeutrino(Vtmp_dauid)
        ){
        hasALepV=true;
        break;
      }
    }

    std::vector<MELAParticle*> modifiedVs;
    MELAParticle* bestVbyMass=nullptr;
    float bestVMassDiff=-1;
    for (MELAParticle* Vtmp:associatedVs){
      if (
        hasALepV &&
        PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
        ){
        for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(false);
        modifiedVs.push_back(Vtmp);
      }
      else if (!bestVbyMass || fabs(Vtmp->m()-PDGHelpers::Wmass)<bestVMassDiff){
        bestVMassDiff = fabs(Vtmp->m()-PDGHelpers::Wmass);
        bestVbyMass = Vtmp;
      }
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
    for (MELAParticle* Vtmp:modifiedVs){ for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(true); }
  }
  void GMECBlock::updateMELAClusters_BestLOAssociatedVBF(const string clustertype, bool isGen){
    // Same as updateMELAClusters_NoInitialGNoAssociatedG, but keep a separate function for future studies
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Manipulate the candidate
    // Assign 0 to the id of gluon mothers
    std::vector<int> motherIds;
    std::vector<int> ajetIds;
    for (int imot=0; imot<melaCand->getNMothers(); imot++){
      motherIds.push_back(melaCand->getMother(imot)->id);
      if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
    }
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
      ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
      if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
    }

    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    // Restore the candidate properties
    for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
    for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  }
  // ME computations that can approximate the NLO QCD (-/+ MiNLO extra jet) phase space to LO QCD in signal VBF or VH
  // Use these for POWHEG samples
  // MELACandidateRecaster has very specific use cases, so do not use these functions for other cases.
  void GMECBlock::updateMELAClusters_BestNLOVHApproximation(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Check if any clusters request this computation
    bool clustersRequest=false;
    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        clustersRequest=true;
        break;
      }
    }
    if (!clustersRequest) return;

    // Need one recaster for each of ZH and WH, so distinguish by the cluster name
    TVar::Production candScheme;
    if (clustertype.find("BestNLOZHApproximation")!=string::npos) candScheme = TVar::Had_ZH;
    else if (clustertype.find("BestNLOWHApproximation")!=string::npos) candScheme = TVar::Had_WH;
    else return;

    MELACandidateRecaster recaster(candScheme);
    MELACandidate* candModified=nullptr;
    MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(melaCand, candScheme);
    if (bestAV){
      recaster.copyCandidate(melaCand, candModified);
      recaster.deduceLOVHTopology(candModified);
      melaHandle->setCurrentCandidate(candModified);
    }
    else return; // No associated Vs found. The algorithm won't work.

    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    delete candModified;
    melaHandle->setCurrentCandidate(melaCand); // Go back to the original candidate
  }
  void GMECBlock::updateMELAClusters_BestNLOVBFApproximation(const string clustertype, bool isGen){
    MELACandidate* melaCand = melaHandle->getCurrentCandidate();
    if (!melaCand) return;

    // Check if any clusters request this computation
    bool clustersRequest=false;
    std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        clustersRequest=true;
        break;
      }
    }
    if (!clustersRequest) return;

    // Need one recaster for VBF
    TVar::Production candScheme;
    if (clustertype.find("BestNLOVBFApproximation")!=string::npos) candScheme = TVar::JJVBF;
    else return;

    MELACandidateRecaster recaster(candScheme);
    MELACandidate* candModified=nullptr;
    recaster.copyCandidate(melaCand, candModified);
    recaster.reduceJJtoQuarks(candModified);
    melaHandle->setCurrentCandidate(candModified);

    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then update the cluster
        theCluster->update();
      }
    }

    delete candModified;
    melaHandle->setCurrentCandidate(melaCand); // Go back to the original candidate
  }


}
