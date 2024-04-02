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
        for (size_t xLoad = 1; xLoad < state.intParms.numInfectionLoadBuckets; xLoad++) {
            double popAtLoad = state.infecteds->getIndex(xAge, xLoad);
            if (xAge == 0) {
                assert(popAtLoad == 0.0);
            }
            total += popAtLoad * state.intParms.deltaTime * state.intParms.deltaLogInfection;
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
