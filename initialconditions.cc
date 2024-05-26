#include "initialconditions.h"

#include <cassert>
#include <cmath>

double intrinsicGrowthRate(const ModelParams& modParms) {
    return log(1/modParms.aveInitInfectionLoad)/modParms.aveInfectiousPeriod;
}

// TODO: which paper?
// Gamma distribution from paper. g an c are scale and shape params. 
double gammaDist(double infectionLoad, double aveInitInfectionLoad, double shapeParam) {
    double g = aveInitInfectionLoad/shapeParam;
    double numerator = pow(infectionLoad/g, shapeParam-1) * exp(-infectionLoad/g);
    double denominator = g * tgamma(shapeParam);
    return numerator/denominator;
}

//Set up the initial distribution of susceptibles. For now we set up with Weibull distribution, in accorance with Weibull death.
void setStartingDistribution(const State& state) {
    std::vector<double>& dist = *state.susceptibles->getCurrentState();
    assert(dist.size() == state.compParms.ageSize);
    assert(state.compParms.lambda != 0.0);
    double totalSusPop = 0.0;
    for (uint32_t i = 0; i < state.compParms.ageSize; i++) {
        dist[i] = weibullOfAge(state.intParms.deltaTime * i, state.compParms.lambda, state.modParms.kappa)
         * state.modParms.initialSusceptiblePop;
        totalSusPop += dist[i] * state.intParms.deltaTime;
    }
    // Normalize it. Weibull of age is a survival fcn, so it won't integrate to 1.
    double normalizer = state.modParms.initialSusceptiblePop/totalSusPop;
    for (uint32_t i = 0; i < state.compParms.ageSize; i++) {
        dist[i] *= normalizer;
    }
}

//This is a Weibull survival function (1 - CDF of Weibull). It's the steady-state distribution
//if we assume death rate as a fcn of age is Weibull distributed.
double weibullOfAge(double age, double lambda, double kappa){
    return exp(-pow(age * lambda, kappa));
}

//TO DO: 
//We put some infected density initially. To avoid discontinuities, it's added
//in the same way as the transfer function from susceptibles
void setInitialInfecteds(const State& state) {
    double deltaArea;
    //Shorthand variables to make code less verbose
    auto& infecteds = *state.infecteds;
    std::vector<double>& loadVec = *state.compParms.columnLoads, susceptibles = *state.susceptibles->getCurrentState();
    double shapeParam = state.modParms.gammaShapeParam, numSusceptibles = state.modParms.initialSusceptiblePop,
        numInfecteds = state.modParms.initialInfectedPop; 
    //Counter for debugging
    double totalInfected = 0.0;
    //double deltaArea = state.intParms.deltaTime * deltaInfection; //wrong bc transformation not measure preserving
    for(uint32_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        // Get proportion that lands in this infection level
        double gammaVal = gammaDist(loadVec[xLoad], state.modParms.aveInitInfectionLoad, shapeParam);
        deltaArea = state.intParms.deltaTime * state.compParms.deltaInfectionForLoad->at(xLoad);
        // Loop over infection levels. Skip age 0 -- no infections
        for(size_t xAge = 1; xAge < state.compParms.ageSize; xAge++) {
            // Get normalized density of susceptibles at this age.
            double weight = susceptibles[xAge] / numSusceptibles;
            // Get density of new infecteds at this bucket & age.
            double infectedVal = weight * gammaVal * numInfecteds;
            totalInfected += infectedVal * deltaArea;
            infecteds.setIndex(xAge, xLoad,  infectedVal);
        }
    }
}

