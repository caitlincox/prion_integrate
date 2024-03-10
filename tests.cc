#include "tests.h"

#include <cassert>

void verifyScuseptibleTotal(const State& state) {
    // TODO: write me.
}

void verifyInfectedTotal(const State& state) {
    // TODO: write me.
}

void assertAproxEqual(double a, double b) {
    assert(std::abs(a - b) < 0.0001 * (a + b)/2.0);
}

void runTests(const State& state, bool expectConstantPop) {
    verifyScuseptibleTotal(state);
    verifyInfectedTotal(state);
    if (expectConstantPop) {
        double expectedPop = state.modParms.initialSusceptiblePopSize +
                state.modParms.initialInfectedPopSize;
        double actualPop = state.compParms.susceptiblePop + state.compParms.infectedPop;
        assertAproxEqual(expectedPop, actualPop);
    }
}
