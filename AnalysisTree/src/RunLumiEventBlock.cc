#include "RunLumiEventBlock.h"


RunLumiEventBlock::RunLumiEventBlock(){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
}
RunLumiEventBlock::RunLumiEventBlock(unsigned int RunNumber_, unsigned int LuminosityBlock_, unsigned long long EventNumber_){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) NAME = NAME##_;
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
}
RunLumiEventBlock::RunLumiEventBlock(RunLumiEventBlock const& other){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
}
bool RunLumiEventBlock::operator==(RunLumiEventBlock const& other) const{
  bool res = true;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) res &= (NAME == other.NAME);
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
  return res;
}
