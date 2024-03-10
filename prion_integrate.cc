#include "death.h"
#include "integrator.h"

int main() {
    IntegrationParams integrationParams = {
        .deltaTime = 0.01,
        .numInfectionLoadBuckets = 100,
        .deltaLogInfection = 0.1,
        .totalTime = 15.0,
    };
    ModelParams modelParams = {
        .maxAge = 10.0,
        .aveLifespan = 5.0,
        .aveInfectiousPeriod = 2.0,
        .aveInitInfectionLoad = 0.1,
        .beta = 0.04,
        .kappa = 2.0,
        .initialSusceptiblePop = 10'000.0 - 1.0,
        .initialInfectedPop = 1.0
    };
    auto state = std::make_unique<State>(integrationParams, modelParams);
    auto births = newReplacementBirthScheme();
    auto deaths = std::make_unique<Death>(*state);
    auto integrator = Integrator(std::move(births), std::move(deaths), std::move(state), true);
    integrator.run();
    return 0;
}
