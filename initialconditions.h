#ifndef INIT
#define INIT

#include "state.h"

#include <memory>
#include <vector>

double intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad);
double weibullOfAge(double weibullLambda, double weibullKappa, double age);
void setStartingDistribution(const State& state);
void setInitialInfecteds(const State& state);

#endif
