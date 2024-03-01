#include "initialconditions.h"
#include <cmath>

// static
double InitialConditions::intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad) {
    return log(1/aveInitInfectionLoad)/aveLifespan;
}

//Set up the initial distribution of susceptibles. For now we set up with Weibull distribution, in accorance with Weibull death.
std::unique_ptr<std::vector<double>> InitialConditions::startingDistribution(double lambda, double kappa, double maxAge, double deltaT,
                                                                             double popSize) {
    int numAgeBuckets = maxAge/deltaT; //Type cast for rounding.
    auto myvec = std::make_unique<std::vector<double>>();
    myvec->resize(numAgeBuckets);
    auto& myptr = *myvec;
    for (int i = 0; i < numAgeBuckets; i++){
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

//THIS BITCH DOES NOT COMPILE. ISN'T DONE YET.
//We put some infected density initially. To avoid discontinuities, it's added in the same way as the transfer function from susceptibles
std::unique_ptr<std::vector<double>> InitialConditions::initialInfecteds(double nunInfecteds, double deltaT, double maxAge, double numInfectionBuckets,
                                                                         double scaleParam, double shapeParam, const std::vector<double>& loadVec){
    int numAgeBuckets = maxAge/deltaT; //Type cast for rounding.
    auto myvec = std::make_unique<std::vector<double>>();
    myvec->resize(numInfectionBuckets * numAgeBuckets);
    for(int i = 0; i < numInfectionBuckets; i++) {
        int bucketAge = deltaT;
        for(int j = 0; j < numAgeBuckets; j++){
            if(1 == 1) return 0;
        }

    }
}

//Gamma distribution from paper. g an c are scale and shape params. 
double gammaDist(double infectionLoad, double aveInitInfectionLoad, double c) {
    double g = aveInitInfectionLoad/c;
    double numerator = pow(infectionLoad/g, c-1) * exp(-infectionLoad/g);
    double denominator = g * tgamma(c);
    return numerator/denominator;
}