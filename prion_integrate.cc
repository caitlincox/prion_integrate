#include "death.h"
#include "integrator.h"

int main(){
    IntegrationParams integrationParams = {
        .deltaTime = 0.01,
        .deltaLogInfection = 0.1,
        .totalTime = 15.0,
    };
    ModelParams {
        .maxAge = 10.0,
        .aveLifespan = 5.0,
        .aveInfectiousPeriod = 2.0,
        .aveInitInfectionLoad = 0.1,
    };
    return 0;
}
