#ifndef SIMPLEAVERAGEINTERFERENCETRIGPHASE_H
#define SIMPLEAVERAGEINTERFERENCETRIGPHASE_H

#include "Discriminant.h"


class SimpleAverageInterferenceTrigPhase : public Discriminant{
protected:
  void eval(const std::vector<float>& vars, const float& valReco);

public:
  SimpleAverageInterferenceTrigPhase();

};

typedef SimpleAverageInterferenceTrigPhase CjjVHint_t;

#endif
