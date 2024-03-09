#ifndef INIT
#define INIT

#include "state.h"

#include <memory>
#include <vector>

double intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad);
double weibullOfAge(double weibullLambda, double weibullKappa, double age);
std::unique_ptr<std::vector<double>> startingDistribution(const State& state);
std::unique_ptr<std::vector<double>> initialInfecteds(const State& state);

#endif
