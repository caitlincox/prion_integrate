#include "death.h"
#include "birth.h"
#include "new_infections.h"
#include "integrator.h"
#include "tests.h"

int main() {
     // TEMP TESTS
    runInitTests();
    //
    
    IntegrationParams integrationParams = {
        .deltaTime = 0.005, //0.02
        .totalTime = 15.0,
        // For now, use ageSize == infectionSize.  There is an argument for
        // 2*ageSize or more, because infected animals can spread disease even
        // if they will die of old age.
        .infectionSizeToAgesSizeRatio = 0.5, // cHANGED!!!!!    
    };
    ModelParams modelParams = {
        .maxAge = 10.0, //10.0
        .aveLifespan = 4.0, //4.0
        .aveInfectiousPeriod = 2.0, //2.0
        .aveInitInfectionLoad = 0.1,
        .beta = 0.04, //0.04
        .kappa = 2.0,
        .initialSusceptiblePop = 10'000.0 - 1.0,
        .initialInfectedPop = 1.0,
        .gammaShapeParam = 5.0
    };
    auto state = std::make_unique<State>(integrationParams, modelParams);
    auto births = newReplacementBirthScheme();
    auto deaths = std::make_unique<Death>(*state);
    auto infections = std::make_unique<NewInfections>(*state);
    auto integrator = Integrator(std::move(births), std::move(deaths), std::move(state), std::move(infections), true);
    const State& state_ref = integrator.getState();
    state_ref.writeInfectedsPGM("initial_infecteds.pgm");
    state_ref.writeSusceptiblesPBM("initial_suseptibles.pbm", 2);
    const State& state_ref_3 = integrator.getState();
    verifySusceptibleTotal(state_ref_3, 9999.0);
    double yo = 0;
    integrator.run();
    return 0;
}
