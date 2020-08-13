#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"


bool VertexSelectionHelpers::testGoodVertex(reco::Vertex const& vtx){
  return (vtx.isValid() && !vtx.isFake() && vtx.ndof()>vtx_ndof_thr && vtx.position().Rho()<=vtx_rho_thr && std::abs(vtx.position().Z())<=vtx_z_thr);
}
