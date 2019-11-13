#ifndef CMS3_VERTEXSELECTIONHELPERS_H
#define CMS3_VERTEXSELECTIONHELPERS_H

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


namespace VertexSelectionHelpers{
  constexpr double vtx_rho_thr = 2.;
  constexpr double vtx_z_thr = 24.;
  constexpr double vtx_ndof_thr = 4.;

  bool testGoodVertex(reco::Vertex const& vtx);
}


#endif
