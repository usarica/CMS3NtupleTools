#ifndef CMS3MELAHELPERS_H
#define CMS3MELAHELPERS_H

#include <MelaAnalytics/GenericMEComputer/interface/GMECHelperFunctions.h>
#include <MelaAnalytics/EventContainer/interface/MELAEvent.h>
#include <string>
#include <vector>
#include <memory>
#include "TTree.h"


// Global helpers
namespace CMS3MELAHelpers{
  extern std::shared_ptr<Mela> melaHandle;

  void setupMela(int sqrts, float mh, TVar::VerbosityLevel verbosity);
  void clearMela();

  using namespace BranchHelpers;
  class GMECBlock{
  protected:
    std::vector<MELAOptionParser*> recome_originalopts; // Be careful: Only for reading
    std::vector<MELAOptionParser*> recome_copyopts;
    std::vector<MELAHypothesis*> recome_units;
    std::vector<MELAHypothesis*> recome_aliased_units;
    std::vector<MELAComputation*> recome_computers;
    std::vector<MELACluster*> recome_clusters;
    std::vector<MELABranch*> recome_branches;

    std::vector<MELAOptionParser*> lheme_originalopts; // Be careful: Only for reading
    std::vector<MELAOptionParser*> lheme_copyopts;
    std::vector<MELAHypothesis*> lheme_units;
    std::vector<MELAHypothesis*> lheme_aliased_units;
    std::vector<MELAComputation*> lheme_computers;
    std::vector<MELACluster*> lheme_clusters;
    std::vector<MELABranch*> lheme_branches;

    std::vector<TTree*> reftrees;

  public:
    GMECBlock(){}
    ~GMECBlock(){}

    void setRefTrees(std::vector<TTree*> const& reftrees_);

    void buildMELABranches(std::vector<std::string> const& lheMElist, std::vector<std::string> const& recoMElist);
    void bookMELABranches(MELAOptionParser* me_opt, MELAComputation* computer, bool doCopy);
    void clearMELABranches();

    void computeMELABranches(bool isGen);
    void pushMELABranches(bool isGen);

    void updateMELAClusters_Common(const std::string clustertype, bool isGen);
    void updateMELAClusters_J1JECJER(const std::string clustertype, bool isGen);
    void updateMELAClusters_J2JECJER(const std::string clustertype, bool isGen);
    void updateMELAClusters_LepWH(const std::string clustertype, bool isGen);
    void updateMELAClusters_LepZH(const std::string clustertype, bool isGen);
    void updateMELAClusters_NoInitialQ(const std::string clustertype, bool isGen);
    void updateMELAClusters_NoInitialG(const std::string clustertype, bool isGen);
    void updateMELAClusters_NoAssociatedG(const std::string clustertype, bool isGen);
    void updateMELAClusters_NoInitialGNoAssociatedG(const std::string clustertype, bool isGen);
    void updateMELAClusters_BestLOAssociatedZ(const std::string clustertype, bool isGen);
    void updateMELAClusters_BestLOAssociatedW(const std::string clustertype, bool isGen);
    void updateMELAClusters_BestLOAssociatedVBF(const std::string clustertype, bool isGen);
    void updateMELAClusters_BestNLOVHApproximation(const std::string clustertype, bool isGen);
    void updateMELAClusters_BestNLOVBFApproximation(const std::string clustertype, bool isGen);
  };

}


#endif
