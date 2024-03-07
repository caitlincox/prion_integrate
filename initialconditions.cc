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

//We put some infected density initially. To avoid discontinuities, it's added in the same way as the transfer function from susceptibles
std::unique_ptr<std::vector<double>> InitialConditions::initialInfecteds(double numInfecteds, double deltaT, double maxAge, double numInfectionBuckets,
                                                                         double scaleParam, double shapeParam, const std::vector<double>& loadVec,
                                                                          std::vector<double>& susceptibles, double numSusceptibles, double aveLoad,
                                                                          double c){
    int numAgeBuckets = maxAge/deltaT; //Type cast for rounding.
    auto myvec = std::make_unique<std::vector<double>>();
    myvec->resize(numInfectionBuckets * numAgeBuckets);
    auto& myptr = *myvec;
    for(int i = 0; i < numInfectionBuckets; i++) { //loop over infection levels
        double gammaVal = gammaDist(loadVec[i], aveLoad, c); //Get perportion that lands in this infection level
        for(int j = 0; j < numAgeBuckets; j++){ //loop over ages
            if(j == 0 | i == 0) myptr[i * numInfectionBuckets] = 0; //if age = 0 or infection load = 0, no infecteds
            else {
                double weight = susceptibles[i] / numSusceptibles; //get perportion of susceptibles at this age
                double infectedVal = weight * gammaVal * numInfecteds; //get perportion of new infecteds at this bucket & age
                myptr[i * numInfectionBuckets + j] = infectedVal; //put that value in the vector
            }
        }
    }
    return myvec;
}

//Gamma distribution from paper. g an c are scale and shape params. 
double gammaDist(double infectionLoad, double aveInitInfectionLoad, double c) {
    double g = aveInitInfectionLoad/c;
    double numerator = pow(infectionLoad/g, c-1) * exp(-infectionLoad/g);
    double denominator = g * tgamma(c);
    return numerator/denominator;
}