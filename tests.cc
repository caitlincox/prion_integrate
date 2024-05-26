#include "tests.h"
#include "initialconditions.h"

#include <cassert>

void assertAproxEqual(double a, double b) {
    assert(std::abs(a - b) <= 0.01 * (a + b)/2.0);
}

void verifySusceptibleTotal(const State& state, double expectedTotal) {
    double total = 0.0;
    std::vector<double>& dist = *state.susceptibles->getCurrentState();
    for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        total += dist[xAge];
    }
    assertAproxEqual(total * state.intParms.deltaTime, expectedTotal);
}

void verifyInfectedTotal(const State& state, double expectedTotal) {
    double total = 0.0;
    double deltaInfection;
    for (size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        double currentLoad = state.compParms.columnLoads->at(xLoad);
        if(xLoad == state.compParms.infectionSize - 1) {
            deltaInfection = 1 - currentLoad;
        }else {
            deltaInfection = state.compParms.columnLoads->at(xLoad + 1) - currentLoad;
        }
        for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            double popDensityAtAge = state.infecteds->getIndex(xAge, xLoad);
            if (xAge == 0) {
                assert(popDensityAtAge == 0.0);
            }
            total += popDensityAtAge * state.intParms.deltaTime * deltaInfection;
        }
    }
    assertAproxEqual(total, expectedTotal);
}

void runTests(const State& state, bool expectConstantPop) {
    if (expectConstantPop) {
        double expectedPop = state.modParms.initialSusceptiblePop +
                state.modParms.initialInfectedPop;
        double actualPop = state.compParms.susceptiblePop + state.compParms.infectedPop;
        assertAproxEqual(expectedPop, actualPop);
    }
}

//Tests whether a distribution integrates to 1
void testDistributionMeasure(double (*dist) (double, double, double), double param1, 
                      double param2, double deltaX, double xMin, double xMax) {
    int numSlices = (xMax - xMin) / deltaX;
    double total = 0;
    double x = xMin;
    for(int i = 0; i < numSlices; i++){
        total += dist(x, param1, param2) * deltaX;
        x += deltaX;
    }
    assertAproxEqual(total, 1.0);
}


void runInitTests() {
    //testDistributionMeasure(&weibullOfAge, 0.2216, 2, 0.02, 0, 10);
    testDistributionMeasure(&gammaDist, 0.10, 5, 0.001, 0, 1);
};

void testInfection(State&state, NewInfections& infections) {
    testInfectionNumbers(state, infections);
    for (int i = 0; i < state.compParms.ageSize; i++) {
        double susNum = infections.getSusceptibles()->getCurrentState()->at(i) * state.intParms.deltaTime;
        testRowSum(*infections.getInfecteds(), i, susNum, state);
    }
}

// Tests whether the moved susceptibles = moved infecteds
void testInfectionNumbers(State& state, NewInfections& infections) {
    double susTotal = 0;
    for(int xAge = 0; xAge < state.compParms.ageSize; xAge++){
        susTotal += infections.getSusceptibles()->getCurrentState()->at(xAge) * state.intParms.deltaTime;
    }
    double infTotal = 0.0;
    double deltaInfection;
    for (size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        double currentLoad = state.compParms.columnLoads->at(xLoad);
        if(xLoad == state.compParms.infectionSize - 1) {
            deltaInfection = 1 - currentLoad;
        }else {
            deltaInfection = state.compParms.columnLoads->at(xLoad + 1) - currentLoad;
        }
        for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            double popDensityAtAge = infections.getInfecteds()->getIndex(xAge, xLoad);
            if (xAge == 0) {
                assert(popDensityAtAge == 0.0);
            }
            infTotal += popDensityAtAge * state.intParms.deltaTime * deltaInfection;
        }
    }
    assertAproxEqual(susTotal, infTotal);
}

// Tests whether an age row of infecteds (or infection) adds up to intended density
void testRowSum(Infecteds& infecteds, size_t ageIndex, double expectedSum, State& state) {
    double sum = 0;
    for (int i = 0; i < state.compParms.infectionSize; i++) {
        double deltaLoad = state.compParms.deltaInfectionForLoad->at(i);
        sum += infecteds.getIndex(ageIndex, i) * state.intParms.deltaTime * deltaLoad;
    }
    assertAproxEqual(sum, expectedSum);
}

// Tests if density of infecteds sums within a reasonable tolerance given our step sizes
// by equally distributing individuals over the range of infection loads we look at
void testInfectedsDensity(State& state) {
        Infecteds ass = Infecteds(state.compParms.ageSize, state.compParms.infectionSize); //shit's fucked
        for (int i = 0; i < state.compParms.infectionSize; i++) { 
            for (int j = 0; j < state.compParms.ageSize; j++) {
                ass.setIndex(j, i, 1);
            }
        }
}