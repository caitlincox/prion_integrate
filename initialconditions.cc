#include "initialconditions.h"

#include <cassert>
#include <cmath>

double intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad) {
    return log(1/aveInitInfectionLoad)/aveLifespan;
}

// TODO: which paper?
// Gamma distribution from paper. g an c are scale and shape params. 
double gammaDist(double infectionLoad, double aveInitInfectionLoad, double c) {
    double g = aveInitInfectionLoad/c;
    double numerator = pow(infectionLoad/g, c-1) * exp(-infectionLoad/g);
    double denominator = g * tgamma(c);
    return numerator/denominator;
}

//Set up the initial distribution of susceptibles. For now we set up with Weibull distribution, in accorance with Weibull death.
void setStartingDistribution(const State& state) {
    std::vector<double>& dist = *state.susceptibles->getCurrentState();
    assert(dist.size() == state.compParms.ageSize);
    double lambda = state.compParms.lambda;
    assert(state.compParms.lambda != 0.0);
    for (int i = 0; i < state.compParms.ageSize; i++) {
        dist[i] = weibullOfAge(state.compParms.lambda, state.modParms.kappa,
                state.intParms.deltaTime * i) * state.modParms.initialSusceptiblePop;
    }
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
void setInitialInfecteds(const State& state) {
    auto& infecteds = *state.infecteds;
    std::vector<double>& loadVec = *state.compParms.columnLoads;
    double c = state.compParms.intrinsicGrowthRate;
    std::vector<double>& susceptibles = *state.susceptibles->getCurrentState();
    size_t numSusceptibles = state.modParms.initialSusceptiblePop;
    size_t numInfecteds = state.intParms.numInfectionLoadBuckets;
    // loop over ages
    for(int xLoad = 0; xLoad < state.intParms.numInfectionLoadBuckets; xLoad++) {
        // Get proportion that lands in this infection level
        // TODO: The gamma distribution function needs work: It needs more
        // parameters to cause the total infecteds to add to numInfecteds.
        // double gammaVal = gammaDist(loadVec[xLoad], state.compParms.aveLoad, c);
        //
        // For now, this should cause the initial infecteds to add to numInfecteds:
        // This is just a constant amount in each bucket.
        double gammaVal = 1.0 / state.intParms.numInfectionLoadBuckets;
        // Loop over infection levels.  if age = 0, no infecteds.
        for(int xAge = 1; xAge < state.compParms.ageSize; xAge++) {
            // Get proportion of susceptibles at this age.
            double weight = susceptibles[xAge] / (numSusceptibles - susceptibles[0]);
            // Get proportion of new infecteds at this bucket & age.
            double infectedVal = weight * gammaVal * numInfecteds;
            infecteds.setIndex(xAge, xLoad,  infectedVal);
        }
    }
}
