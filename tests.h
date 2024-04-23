#ifndef TESTS
#define TESTS

#include "state.h"

void runTests(const State& state, bool expectConstantPop);
void verifySusceptibleTotal(const State& state, double expectedTotal);
void verifyInfectedTotal(const State& state, double expectedTotal);
void testDistributionMeasure(double (&dist) (double, double, double), double param1, 
                      double param2, double deltaX, double xMin, double xMax);
void runInitTests();

#endif
