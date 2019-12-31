#ifndef SIMPLEINTERFERENCETRIGPHASE_H
#define SIMPLEINTERFERENCETRIGPHASE_H

#include "Discriminant.h"


class SimpleInterferenceTrigPhase : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  SimpleInterferenceTrigPhase();
};

typedef SimpleInterferenceTrigPhase Cintkin_t;
typedef SimpleInterferenceTrigPhase Cggint_t;
typedef SimpleInterferenceTrigPhase CjjVBFint_t;
typedef SimpleInterferenceTrigPhase CjjVHL1ZGsint_t;

#endif
