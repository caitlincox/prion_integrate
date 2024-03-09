#include "death.h"

#include <cmath>
#include <cstdio>

//Constructor
Death::Death(const State& state) {
    weibullSurvivorshipKappa_ = state.modParms.kappa;
    weibullSurvivorshipLambda_ = state.compParms.lambda;
}

//Given age, calculate survival coefficient
double Death::weibullDeathRate(double age) {
    return(weibullSurvivorshipKappa_ * weibullSurvivorshipLambda_ * pow(weibullSurvivorshipLambda_ * age, weibullSurvivorshipKappa_ - 1));
}
