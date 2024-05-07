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
    double lambda = state.compParms.lambda;
    assert(state.compParms.lambda != 0.0);
    double totalSusPop = 0.0;
    for (int i = 0; i < state.compParms.ageSize; i++) {
        dist[i] = weibullOfAge(state.intParms.deltaTime * i, state.compParms.lambda, state.modParms.kappa)
         * state.modParms.initialSusceptiblePop;
        totalSusPop += dist[i] * state.intParms.deltaTime;
    }
    // Normalize it.
    double normalizer = state.modParms.initialSusceptiblePop/totalSusPop;
    for (int i = 0; i < state.compParms.ageSize; i++) {
        dist[i] *= normalizer;
    }
}

//The way Stringer et al. did initial distribution was using a steady state Weibull. 
//This is just an implementation of the Weibull distribution density function.
//double weibullOfAge(double lambda, double kappa, double age) {
 //   double ratio1 = kappa/lambda;
 //   double ratio2 = age/lambda;
 //   double firstComponent = ratio1 * pow(ratio2, kappa - 1); //breaking the calculation into parts
 //   double secondComponent = -pow(ratio2,kappa);
 //   return firstComponent * exp(secondComponent);
//}

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
    for(int xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        // Get proportion that lands in this infection level
        // TODO: The gamma distribution function needs work: It needs more
        // parameters to cause the total infecteds to add to numInfecteds.
        double gammaVal = gammaDist(loadVec[xLoad], state.modParms.aveInitInfectionLoad, shapeParam);
        if(xLoad == state.compParms.infectionSize - 1) {
            deltaArea = state.intParms.deltaTime * (1 - loadVec[xLoad]);
        }else {
            deltaArea = state.intParms.deltaTime * (loadVec[xLoad + 1] - loadVec[xLoad]);
            }
        // Loop over infection levels. Skip age 0 -- no infections
        for(int xAge = 1; xAge < state.compParms.ageSize; xAge++) {
            // Get normalized density of susceptibles at this age.
            double weight = susceptibles[xAge] / numSusceptibles;
            // Get density of new infecteds at this bucket & age.
            double infectedVal = weight * gammaVal * numInfecteds;
            totalInfected += infectedVal * deltaArea;
            infecteds.setIndex(xAge, xLoad,  infectedVal);
        }
    }
}

