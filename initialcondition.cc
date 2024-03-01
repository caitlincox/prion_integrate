#include "initialconditions.h"
#include <cmath>
double InitialConditions::intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad) {
    return log(1/aveInitInfectionLoad)/aveLifespan;
}

//Set up the initial distribution of susceptibles. For now we set up with Weibull distribution, in accorance with Weibull death.
std::unique_ptr<std::vector<double>> InitialConditions::startingDistribution(double lambda, double kappa, double maxAge, double deltaT,
                                                                             double popSize) {
    int numBuckets = maxAge/deltaT; //Type cast for rounding.
    auto myvec = std::make_unique<std::vector<double>>();
    myvec->resize(numBuckets);
    auto myptr = *myvec;
    for (int i = 0; i < numBuckets, i++;){
        myptr[i] = weibullOfAge(lambda, kappa, deltaT * i) * popSize;
    }
    return myvec;
}

//The way Stringer et al. did initial distribution was using a steady state Weibull. 
//This is just an implementation of the Weibull distribution density function.
double InitialConditions::weibullOfAge(double lambda, double kappa, double age){
    double ratio1 = kappa/lambda;
    double ratio2 = age/lambda;
    double firstComponent = ratio1 * pow(ratio2, kappa - 1); //breaking the calculation into parts
    double secondComponent = -pow(ratio2,kappa);
    return firstComponent * exp(secondComponent);
}