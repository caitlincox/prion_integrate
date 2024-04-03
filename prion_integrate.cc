#include "death.h"
#include "integrator.h"

int main() {
    IntegrationParams integrationParams = {
        .deltaTime = 0.1,
        .totalTime = 15.0,
        // For now, use ageSize == infectionSize.  There is an argument for
        // 2*ageSize or more, because infected animals can spread disease even
        // if they will die of old age.
        .infectionSizeToAgesSizeRatio = 1.0,
    };
    ModelParams modelParams = {
        .maxAge = 10.0,
        .aveLifespan = 4.0,
        .aveInfectiousPeriod = 2.0,
        .aveInitInfectionLoad = 0.1,
        .beta = 0.04,
        .kappa = 2.0,
        .initialSusceptiblePop = 10'000.0 - 1.0,
        .initialInfectedPop = 1.0,
        .gammaShapeParam = 5.0
    };
    auto state = std::make_unique<State>(integrationParams, modelParams);
    auto births = newReplacementBirthScheme();
    auto deaths = std::make_unique<Death>(*state);
    auto integrator = Integrator(std::move(births), std::move(deaths), std::move(state), true);
    const State& state_ref = integrator.getState();
    state_ref.writeInfectedsPGM("initial_infecteds.pgm");
    state_ref.writeSusceptiblesPBM("initial_suseptibles.pbm", 2);
// temp
return 0;
    integrator.run();
    return 0;
}
