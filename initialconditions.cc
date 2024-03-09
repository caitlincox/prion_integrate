#include "initialconditions.h"
#include <cmath>

double intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad) {
    return log(1/aveInitInfectionLoad)/aveLifespan;
}


//Gamma distribution from paper. g an c are scale and shape params. 
double gammaDist(double infectionLoad, double aveInitInfectionLoad, double c) {
    double g = aveInitInfectionLoad/c;
    double numerator = pow(infectionLoad/g, c-1) * exp(-infectionLoad/g);
    double denominator = g * tgamma(c);
    return numerator/denominator;
}

//Set up the initial distribution of susceptibles. For now we set up with Weibull distribution, in accorance with Weibull death.
std::unique_ptr<std::vector<double>> startingDistribution(const State& state) {
    auto dist = std::make_unique<std::vector<double>>();
    dist->resize(state.compParms.ageSize);
    auto& distRef = *dist;
    for (int i = 0; i < state.compParms.ageSize; i++) {
        distRef[i] = weibullOfAge(state.compParms.lambda, state.modParms.kappa, state.intParms.deltaTime * i) * state.compParms.popSize;
    }
    return dist;
}

//The way Stringer et al. did initial distribution was using a steady state Weibull. 
//This is just an implementation of the Weibull distribution density function.
double weibullOfAge(double lambda, double kappa, double age) {
    double ratio1 = kappa/lambda;
    double ratio2 = age/lambda;
    double firstComponent = ratio1 * pow(ratio2, kappa - 1); //breaking the calculation into parts
    double secondComponent = -pow(ratio2,kappa);
    return firstComponent * exp(secondComponent);
}

//We put some infected density initially. To avoid discontinuities, it's added
//in the same way as the transfer function from susceptibles
std::unique_ptr<std::vector<double>> initialInfecteds(const State& state) {
    auto dist = std::make_unique<std::vector<double>>();
    dist->resize(state.intParms.numInfectionLoadBuckets * state.compParms.ageSize);
    auto& distRef = *dist;
    auto& loadVec = *state.compParms.columnLoads;
    double c = state.compParms.intrinsicGrowthRate;
    std::vector<double>& susceptibles = *state.susceptibles->getCurrentState();
    size_t numSusceptibles = susceptibles.size();
    size_t numInfecteds = state.compParms.infectedPop;
    for(int i = 0; i < state.intParms.numInfectionLoadBuckets; i++) { //loop over infection levels
        double gammaVal = gammaDist(loadVec[i], state.compParms.aveLoad, c); //Get perportion that lands in this infection level
        for(int j = 0; j < state.compParms.ageSize; j++) { //loop over ages
            if(j == 0 | i == 0) distRef[i * state.intParms.numInfectionLoadBuckets] = 0; //if age = 0 or infection load = 0, no infecteds
            else {
                double weight = susceptibles[i] / numSusceptibles; //get perportion of susceptibles at this age
                double infectedVal = weight * gammaVal * numInfecteds; //get perportion of new infecteds at this bucket & age
                distRef[i * state.intParms.numInfectionLoadBuckets + j] = infectedVal; //put that value in the vector
            }
        }
    }
    return dist;
}
