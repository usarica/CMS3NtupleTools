#ifndef RUNEVENTLUMIBLOCK_H
#define RUNEVENTLUMIBLOCK_H

#include "CMS3/Dictionaries/interface/CommonTypedefs.h"


#define RUNLUMIEVENT_VARIABLES \
RUNLUMIEVENT_VARIABLE(cms3_runnumber_t, RunNumber, 0) \
RUNLUMIEVENT_VARIABLE(cms3_lumiblock_t, LuminosityBlock, 0) \
RUNLUMIEVENT_VARIABLE(cms3_eventnumber_t, EventNumber, 0)


class RunLumiEventBlock{
protected:
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE

public:
  RunLumiEventBlock();
  RunLumiEventBlock(unsigned int RunNumber_, unsigned int LuminosityBlock_, unsigned long long EventNumber_);
  RunLumiEventBlock(RunLumiEventBlock const& other);

  bool operator==(RunLumiEventBlock const& other) const;

};


#endif
