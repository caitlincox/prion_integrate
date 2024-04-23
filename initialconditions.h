#ifndef INIT
#define INIT

#include "state.h"

#include <memory>
#include <vector>

double intrinsicGrowthRate(const ModelParams& modParms);
double weibullOfAge(double age, double weibullLambda, double weibullKappa);
void setStartingDistribution(const State& state);
void setInitialInfecteds(const State& state);
double gammaDist(double infectionLoad, double aveInitInfectionLoad, double shapeParam);


#endif
