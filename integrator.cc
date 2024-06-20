#include "integrator.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "tests.h"

Integrator::Integrator(
    std::unique_ptr<BirthScheme> births,
    std::unique_ptr<Death> deaths,
    std::unique_ptr<State> state,
    std::unique_ptr<NewInfections> newInfections,
    bool expectConstantPop)
{
    births_ = std::move(births);
    deaths_ = std::move(deaths);
    state_ = std::move(state);
    newInfections_ = std::move(newInfections);
    expectConstantPop_ = expectConstantPop;
}

// Determine if we will pass a time epoch boundary this time step.  This is
// used to do occasional compute-intensive processing such as writing graphs of
// the current state.
uint32_t computeTimeEpoch(float time, float deltaTime, float totalTime, uint32_t numEpochs) {
    return static_cast<uint32_t>(time * numEpochs / totalTime);
}

void Integrator::run() {
    uint32_t epoch = 0;
    std::vector<double> graphSus;
    std::vector<double> graphInf;
    graphSus.push_back(integrateSusceptibles(*state_->susceptibles));
    graphInf.push_back(integrateInfecteds(*state_->infecteds, *state_));
    for (double time = 0.0; time < state_->intParms.totalTime;
           time += state_->intParms.deltaTime) {
        state_->updateComputedParameters();
        uint32_t nextEpoch = computeTimeEpoch(time, state_->intParms.deltaTime,
                state_->intParms.totalTime, 100);
        graphSus.push_back(state_->compParms.susceptiblePop);
        graphInf.push_back(integrateInfecteds(*state_->infecteds, *state_));
        if (nextEpoch != epoch) {
            printf("Writing graph %u\n", epoch);
            std::string epochStr = std::to_string(epoch);
            //printArray(*state_->susceptibles->getCurrentState(), "sus" + std::to_string(time));
            //state_->writeSusceptiblesPBM("data/suseptibles_" + epochStr + ".pbm", 2);
            //state_->writeInfectedsPGM("data/infecteds_" + epochStr + ".pgm");
            if(epoch == 98){
             printArray(*state_->susceptibles->getCurrentState(), "sus_snap");
             printArray(graphSus, "sus.csv");
             printArray(graphInf, "inf.csv");
            }
            epoch = nextEpoch;
        }
        runTests(*state_, expectConstantPop_);
        timeStep();
// temp
//        testInfection(*state_, *newInfections_);
    }
}

// We kill off animals at a max age.  This has contributions from both susceptibles and infecteds.
void Integrator::computeAgeDeaths() {
    ComputedParams& compParms = state_->compParms;
    Infecteds* infecteds = state_->infecteds.get();
    std::vector<double>& susVec = *state_->susceptibles->getCurrentState();
    // Initialize age deaths to the oldest suseptibles.
    compParms.ageDeaths = susVec[compParms.ageSize - 1];
    susVec[compParms.ageSize - 1] = 0.0;
    size_t xMaxAge = compParms.ageSize - 1;
    for (size_t xLoad = 0; xLoad < compParms.infectionSize; xLoad++) {
        compParms.ageDeaths += infecteds->getIndex(xMaxAge, xLoad) * compParms.deltaInfectionForLoad->at(xLoad);
    }
    compParms.ageDeaths *= state_->intParms.deltaTime; //Scale by delta time
}

// Compute deaths from infection.
void Integrator::computeInfectionDeaths() {
    ComputedParams& compParms = state_->compParms;
    Infecteds* infecteds = state_->infecteds.get();
    size_t xMaxAge = compParms.ageSize - 1;
    compParms.infectionDeaths = 0.0;
    size_t xMaxLoad = compParms.infectionSize - 1;
    // The special case of the bucket that would die both from age and
    // infection is added to the age deaths.
    for (size_t xAge = 0; xAge < xMaxAge; xAge++) {
        compParms.infectionDeaths += infecteds->getIndex(xAge, xMaxLoad) * compParms.deltaInfectionForLoad->at(xMaxLoad);
    }
    // Note that infections in a column are from 0 to 1, so total area is just deltaTime.
    compParms.infectionDeaths *= state_->intParms.deltaTime; //scale by delta time
}

// Redundant code. Will refactor, but for now, easier to debug so we don't
// have to pass class members to a class function.
void Integrator::computeNewInfectionDeaths() {
    ComputedParams& compParms = state_->compParms;
    Infecteds* newInfecteds = newInfections_->getInfecteds();
    size_t xMaxAge = compParms.ageSize - 1;
    double newInfectionDeaths = 0.0;
    size_t xMaxLoad = compParms.infectionSize - 1;
    double deltaArea = compParms.deltaInfectionForLoad->at(xMaxLoad) * state_->intParms.deltaTime;
    // The special case of the bucket that would die both from age and
    // infection is added to the age deaths.
    for (size_t xAge = 0; xAge < xMaxAge; xAge++) {
        newInfectionDeaths += newInfecteds->getIndex(xAge, xMaxLoad) * deltaArea;
    }
    // Note that infections in a column are from 0 to 1, so total area is just deltaTime.
    compParms.infectionDeaths += newInfectionDeaths;
}

double Integrator::integrateInfecteds(Infecteds& infecteds, State& state) {
    double total = 0.0;
    for (int xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        double deltaArea = state.intParms.deltaTime * state.compParms.deltaInfectionForLoad->at(xLoad);
        for (int xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            total += infecteds.getIndex(xAge, xLoad) * deltaArea;
        }
    }
    return total;
}

double Integrator::integrateSusceptibles(Susceptibles& susceptibles) {
    double total = 0.0;
    for (int xAge = 0; xAge < state_->compParms.ageSize; xAge++){
        total += susceptibles.getCurrentState()->at(xAge) * state_->intParms.deltaTime;
    }
    return total;
}
//Calculates amount to be moved prior to modifing population
void Integrator::calculateDeltas() {
    // Calculate number to be infected. We do this because we're using forwards
    // Euler. Add the delta back in after the timestep below.
    newInfections_->prepInfecteds(*state_);
    double wot = integrateSusceptibles(*newInfections_->getSusceptibles());
    //Calculate the ones that will die due to exceeding max infection or max age
    computeAgeDeaths();
    computeInfectionDeaths();
    computeNewInfectionDeaths();

    // Now kill off population due to natural deaths.
    // Later would like to pre-calculate deaths.
    deaths_->kill(*state_);
    adjustSusceptibles();
    newInfections_->calculateDeltaInfecteds(*state_);
    adjustInfecteds();
}

// forces to be non-negative
void Integrator::adjustSusceptibles() {
    for(int xAge = 0; xAge < state_->compParms.ageSize; xAge++) {
        if (state_->susceptibles->getCurrentState()->at(xAge) < newInfections_->getSusceptibles()->getCurrentState()->at(xAge)) {
            newInfections_->getSusceptibles()->getCurrentState()->at(xAge) = state_->susceptibles->getCurrentState()->at(xAge);
        }
    }
}

void Integrator::adjustInfecteds() {
    double expectedVal = integrateSusceptibles(*newInfections_->getSusceptibles());
    double outputVal = integrateInfecteds(*newInfections_->getInfecteds(), *state_);
    if(outputVal != 0.0) { 
        newInfections_->getInfecteds()->rescaleInfecteds(expectedVal/outputVal);
    }
    assertAproxEqual(expectedVal, integrateInfecteds(*newInfections_->getInfecteds(), *state_));
}

void Integrator::timeStep() {
    //Nameing things to make more readable
    ComputedParams& compParms = state_->compParms;
    size_t xMaxAge = compParms.ageSize - 1;
    Infecteds* infecteds = state_->infecteds.get();
    std::vector<double>& susVec = *state_->susceptibles->getCurrentState();

    calculateDeltas();

    // Compute births.  Add it in after the time step.
    double births = births_->calculateBirth(*state_);
    newInfections_->moveInfecteds(*state_);

    // Now advance time by deltaTime.
    // Move all susceptiblees to one higher age index.
    double* susPtr = &susVec[0];
    memmove(susPtr + 1, susPtr, xMaxAge * sizeof(double));
    // Move all the infectes both by a deltaTime and increase infection.
    // This simply moves data by 1 in both age and infection dimensions.
    double* infectedPtr = infecteds->getInfectionRow(0);
    size_t numMoving = xMaxAge * compParms.infectionSize - 1;
    memmove(infectedPtr + compParms.infectionSize + 1, infectedPtr, numMoving * sizeof(double));
    // Zero out the min infection buckets and min age buckets.
    for (size_t xAge = 0; xAge < compParms.ageSize; xAge++) {
        infecteds->setIndex(xAge, 0, 0.0);
    }
    for (size_t xLoad = 0; xLoad < compParms.infectionSize; xLoad++) {
        infecteds->setIndex(0, xLoad, 0.0);
    }
    // Add in births to suseptibles.
    if (births < 0) births = 0;
    susVec[0] = births;
    // Add newly infecteds
   double adjustmentFactor = 1.0/state_->compParms.deltaInfection;
    for (size_t xAge = 1; xAge < compParms.ageSize; xAge++) {
        for (size_t xLoad = 1; xLoad < compParms.infectionSize; xLoad++) {
            infecteds->setIndex(xAge, xLoad, infecteds->getIndex(xAge, xLoad) * adjustmentFactor);
       }
    }
}
