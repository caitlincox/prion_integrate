#include "tests.h"

#include <cassert>

void assertAproxEqual(double a, double b) {
    assert(std::abs(a - b) < 0.0001 * (a + b)/2.0);
}

void verifySusceptibleTotal(const State& state, double expectedTotal) {
    double total = 0.0;
    std::vector<double>& dist = *state.susceptibles->getCurrentState();
    for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        total += dist[xAge];
    }
    assertAproxEqual(total, expectedTotal);
}

void verifyInfectedTotal(const State& state, double expectedTotal) {
    // TODO: write me.
}

void runTests(const State& state, bool expectConstantPop) {
    if (expectConstantPop) {
        double expectedPop = state.modParms.initialSusceptiblePop +
                state.modParms.initialInfectedPop;
        double actualPop = state.compParms.susceptiblePop + state.compParms.infectedPop;
        assertAproxEqual(expectedPop, actualPop);
    }
}
