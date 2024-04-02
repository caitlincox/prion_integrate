#include "tests.h"

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
        for (size_t xLoad = 0; xLoad < state.intParms.numInfectionLoadBuckets; xLoad++) {
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
