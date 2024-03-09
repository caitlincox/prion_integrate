#include "death.h"
#include <cmath>
#include <cstdio>

//Constructor
Death::Death(double kappa, double aveLifespan) {
    weibullSurvivorshipKappa_ = kappa;
    weibullSurvivorshipLambda_ = findLambda(aveLifespan, kappa);
}

//Given age, calculate survival coefficient
double Death::weibullDeathRate(double age) {
    return(weibullSurvivorshipKappa_ * weibullSurvivorshipLambda_ * pow(weibullSurvivorshipLambda_ * age, weibullSurvivorshipKappa_ - 1));
}

double Death::findLambda(double aveLifespan, double kappa) {
  printf("Not implemented");
  return 0;
}
