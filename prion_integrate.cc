#include "death.h"
#include "birth.h"
#include "new_infections.h"
#include "integrator.h"
#include "tests.h"

int main() {
    // TEMP TESTS
    runInitTests();
    
    IntegrationParams integrationParams = {
        .deltaTime = 0.01, //0.02
        .totalTime = 10.0,
        // For now, use ageSize == infectionSize.  There is an argument for
        // 2*ageSize or more, because infected animals can spread disease even
        // if they will die of old age.
        .infectionSizeToAgesSizeRatio = 0.5, // cHANGED!!!!!    
    };
    double totalPop = 200;  // 200
    double infectedPop = 1.0; // 1.0
    double susceptiblePop = totalPop - infectedPop;
    ModelParams modelParams = {
        .maxAge = 10.0, //10.0
        .aveLifespan = 5.0, //4.0
        .aveInfectiousPeriod = 2.0, //2.0 //1.8
        .aveInitInfectionLoad = 0.0001, //0.0001
        .beta = 0.05, //0.04
        .kappa = 2.0,
        .initialSusceptiblePop = susceptiblePop,
        .initialInfectedPop = infectedPop,
        .gammaShapeParam = 5.0,
        .infectionDeathLoadMultiplier = 0.0
    };
    auto state = std::make_unique<State>(integrationParams, modelParams);
    auto births = newReplacementBirthScheme();
    auto deaths = std::make_unique<Death>(*state); 
    //auto deaths = std::make_unique<PredatorScheme>(*state, 10.0, 0.5); //0.17, 0.51
    //auto deaths = std::make_unique<sexScheme>(*state, 0.25, 0.25);
    auto infections = std::make_unique<NewInfections>(*state);
    auto integrator = Integrator(std::move(births), std::move(deaths), std::move(state), std::move(infections), true);
    const State& state_ref = integrator.getState();
    verifySusceptibleTotal(state_ref, susceptiblePop);
    integrator.dryRun(20.0); //stabilize distribution
    integrator.run();
    return 0;
}