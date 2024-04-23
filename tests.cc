#include "tests.h"
#include "initialconditions.h"

#include <cassert>

void assertAproxEqual(double a, double b) {
    assert(std::abs(a - b) < 0.01 * (a + b)/2.0);
}

void verifySusceptibleTotal(const State& state, double expectedTotal) {
    double total = 0.0;
    std::vector<double>& dist = *state.susceptibles->getCurrentState();
    for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        total += dist[xAge] * state.intParms.deltaTime;
    }
    assertAproxEqual(total, expectedTotal);
}

void verifyInfectedTotal(const State& state, double expectedTotal) {
    double total = 0.0;
    for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        for (size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
            double popDensityAtLoad = state.infecteds->getIndex(xAge, xLoad);
            if (xAge == 0) {
                assert(popDensityAtLoad == 0.0);
            }
            total += popDensityAtLoad * state.intParms.deltaTime * state.compParms.deltaLogInfection;
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


void runInitTests(){
    testDistributionMeasure(&weibullOfAge, 0.2216, 2, 0.02, 0, 10);
    testDistributionMeasure(&gammaDist, 0.10, 5, 0.001, 0, 1);
};