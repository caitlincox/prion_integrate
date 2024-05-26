#ifndef TESTS
#define TESTS

#include "state.h"
#include "new_infections.h"

void runTests(const State& state, bool expectConstantPop);
void verifySusceptibleTotal(const State& state, double expectedTotal);
void verifyInfectedTotal(const State& state, double expectedTotal);
void testDistributionMeasure(double (&dist) (double, double, double), double param1, 
                      double param2, double deltaX, double xMin, double xMax);
void runInitTests();
void testInfectionNumbers(State& state, NewInfections& newInfections);
void testInfection(State& state, NewInfections& infections);
void testRowSum(Infecteds& infecteds, size_t ageIndex, double expectedSum, State& state);
void testDeltaInfectionsAddToOne(const State& state);


#endif
