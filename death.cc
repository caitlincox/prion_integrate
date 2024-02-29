#include "death.h"
#include <cmath>

//Constructor
Death::Death(double kappa, double aveLifespan){
    weibullSurvivorshipKappa = kappa;
    weibullSurvivorshipLambda = findLambda(aveLifespan, kappa);
}

//Given age, calculate survival coefficient
double Death::weibullDeathRate(double age){
    return(weibullSurvivorshipKappa * weibullSurvivorshipLambda * pow(weibullSurvivorshipLambda * age, weibullSurvivorshipKappa - 1));
}

double Death::findLambda(double aveLifespan, double kappa){return 0;} //placeholder for now