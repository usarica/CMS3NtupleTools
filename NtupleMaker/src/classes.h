#include <vector>
#include "DataFormats/Common/interface/Wrapper.h"
#include "CMS3/NtupleMaker/interface/GenInfo.h"
#include "CMS3/NtupleMaker/interface/IsotrackInfo.h"
#include "CMS3/NtupleMaker/interface/TriggerInfo.h"
#include "CMS3/NtupleMaker/interface/TriggerObjectInfo.h"
#include "CMS3/NtupleMaker/interface/METFilterInfo.h"
#include "CMS3/NtupleMaker/interface/METInfo.h"
#include "CMS3/NtupleMaker/interface/METShiftInfo.h"


namespace{
  struct dummy_CMS3_dict{
    edm::Wrapper<GenInfo> dummy_geninfo_wrapper;

    edm::Wrapper<IsotrackInfo> dummy_isotrackinfo_wrapper;
    edm::Wrapper< std::vector<IsotrackInfo> > dummy_visotrackinfo_wrapper;

    edm::Wrapper<TriggerInfo> dummy_triggerinfo_wrapper;
    edm::Wrapper< std::vector<TriggerInfo> > dummy_vtriggerinfo_wrapper;

    edm::Wrapper<TriggerObjectInfo> dummy_triggerobjectinfo_wrapper;
    edm::Wrapper< std::vector<TriggerObjectInfo> > dummy_vtriggerobjectinfo_wrapper;

    edm::Wrapper<METFilterInfo> dummy_metfilterinfo_wrapper;

    edm::Wrapper<METInfo> dummy_metinfo_wrapper;

    edm::Wrapper<METShiftInfo> dummy_metshiftinfo_wrapper;
  };
}
