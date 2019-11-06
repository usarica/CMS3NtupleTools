#include <vector>
#include "CMS3/NtupleMaker/interface/GenInfo.h"
#include "CMS3/NtupleMaker/interface/TriggerInfo.h"
#include "DataFormats/Common/interface/Wrapper.h"


namespace{
  struct dummy_CMS3_dict{
    edm::Wrapper<GenInfo> dummy_geninfo_wrapper;
    edm::Wrapper<TriggerInfo> dummy_triggerinfo_wrapper;
    edm::Wrapper< std::vector<TriggerInfo> > dummy_vtriggerinfo_wrapper;
  };
}
