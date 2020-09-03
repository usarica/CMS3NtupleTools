#include <CMS3/Dictionaries/interface/JetMETEnums.h>
#include <CMS3/NtupleMaker/interface/METShiftInfo.h>


METShiftInfo::METShiftInfo() :
  metshifts(JetMETEnums::nMETShiftTypes)
{}

METShiftInfo::METShiftInfo(const METShiftInfo& other) :
  metshifts(other.metshifts)
{}
