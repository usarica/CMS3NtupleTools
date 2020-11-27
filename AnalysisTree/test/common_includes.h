#include <cassert>

#include "SamplesCore.h"
#include "RunLumiEventBlock.h"
#include "OffshellSampleHelpers.h"
#include "OffshellTriggerHelpers.h"

#include "SimEventHandler.h"
#include "GenInfoHandler.h"
#include "VertexHandler.h"
#include "MuonHandler.h"
#include "ElectronHandler.h"
#include "PhotonHandler.h"
#include "SuperclusterHandler.h"
#include "JetMETHandler.h"
#include "DileptonHandler.h"
#include "IsotrackHandler.h"
#include "FSRHandler.h"
#include "OverlapMapHandler.h"
#include "PFCandidateHandler.h"
#include "EventFilterHandler.h"

#include "ParticleObjectHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "IsotrackSelectionHelpers.h"
#include "BtagHelpers.h"

#include "MuonScaleFactorHandler.h"
#include "ElectronScaleFactorHandler.h"
#include "PhotonScaleFactorHandler.h"
#include "PUJetIdScaleFactorHandler.h"
#include "BtagScaleFactorHandler.h"
#include "TriggerScaleFactorHandler.h"

#include "METCorrectionHandler.h"

#include "DiscriminantClasses.h"

#include "BaseTreeLooper.h"

#include "HostHelpersCore.h"
#include "HelperFunctions.h"
#include "StatisticsHelpers.h"
#include "MELAStreamHelpers.hh"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"


using namespace std;
using namespace MELAStreamHelpers;
