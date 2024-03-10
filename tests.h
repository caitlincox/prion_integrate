#ifndef TESTS
#define TESTS

#include "state.h"

void runTests(const State& state, bool expectConstantPop);
void verifySusceptibleTotal(const State& state, double expectedTotal);
void verifyInfectedTotal(const State& state, double expectedTotal);

#endif
