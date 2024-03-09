#ifndef INTEGRATOR
#define INTEGRATOR

#include <memory>

#include "birth.h"
#include "death.h"
#include "infecteds.h"
#include "state.h"
#include "susceptibles.h"

//Integrate with forward Euler scheme
class Integrator {
public:
    Integrator(
        std::unique_ptr<BirthScheme> births,
        std::unique_ptr<Death> deaths,
        std::unique_ptr<State> state);
    void run();

private:
    // Update state_ for the current time step.
    void update();
    void findSums();

    IntegrationParams integrationParams_;
    ModelParams modelParams_;
    std::unique_ptr<BirthScheme> births_;
    std::unique_ptr<Death> deaths_;
    std::unique_ptr<State> state_;
    // Sums computed per step.
    double totalInfection_;
    double infectedPop_;
    double susceptiblePop_;
    double transferRate_;
    // Computed in initialcontinsion.cc.
    double intrinsicGrowthRate_;
    double firstBucketLogLoad_;
    // These are the infections loads for column i in the infecteds table.
    // The first column has all 0's meaning no infecteds have zero load.
    // Also there are no 0-age infectes, so the 0 row is also 0's.
    std::unique_ptr<std::vector<double>> columnLoads_;
};

#endif
