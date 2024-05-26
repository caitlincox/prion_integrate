#include "state.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio> 

#include "initialconditions.h"
#include "tests.h"

namespace {

// Weibull survivorship lambda
double findLambda(double aveLifespan, double kappa) {
  double lamNumerator = tgamma((kappa + 1)/kappa);
  return lamNumerator / aveLifespan;
}

}

// Initialize computed parameters, including those computed up-front.
void State::initializeComputedParameters() {
    compParms.intrinsicGrowthRate = intrinsicGrowthRate(modParms);
    compParms.deltaLogInfection = compParms.intrinsicGrowthRate * intParms.deltaTime;
    compParms.ageSize = modParms.maxAge / intParms.deltaTime;
    compParms.infectionSize = compParms.ageSize * intParms.infectionSizeToAgesSizeRatio;
    compParms.lambda = findLambda(modParms.aveLifespan, modParms.kappa);
    setStartingDistribution(*this);
    // We represent load buckets on a natural log scale.  This model throws
    // away loads that are very small, which otherwise would show up as animals
    // with longer incubation periods.  To model longer incubation periods,
    // increase the number of infection load buckets.
    //
    // The first bucket at load index 0 has load e^w, the minimum infection
    // load modeled.  We have to compute this bucket's load.  Call 'b' the
    // number of infection load buckets.  When load reaches 1, the animal dies,
    // so the max bucket is the one just before reaching 1.  We can compute 1 =
    // e^w * e^(C*dt*b) => 0 = w + C*dt*b => w = -C*dt*b.  The first bucketk
    // has load e^w.  After that the load is (e^w)*e^(C*dt*bucketIndex),
    // where bucketIndex is on the infection load axis.  At bucketIndex == 0,
    // we have the minimum represented load of e^w.  The bucket index with
    // infection load == 1 is not represented, as the animals would be dead.
    compParms.firstBucketLogLoad = -compParms.intrinsicGrowthRate *
        intParms.deltaTime * compParms.infectionSize;
    compParms.columnLoads = std::make_unique<std::vector<double>>();
    compParms.columnLoads->resize(compParms.infectionSize);
    // Compute the infection load for each bucket.
    auto& columnLoads = *compParms.columnLoads;
    // First bucket is e^w.
    columnLoads[0] = exp(compParms.firstBucketLogLoad);
    for (size_t i = 0; i < compParms.infectionSize; i++) {
        columnLoads[i] = exp(compParms.firstBucketLogLoad + i * intParms.deltaTime *
            compParms.intrinsicGrowthRate);
    }
    //Now create a delta load vector as well
    compParms.deltaInfectionForLoad = std::make_unique<std::vector<double>>();
    compParms.deltaInfectionForLoad->resize(compParms.infectionSize);
    auto& deltaInfection = *compParms.deltaInfectionForLoad;
    deltaInfection[compParms.infectionSize - 1] = 1 - columnLoads[compParms.infectionSize - 1];
    for (size_t i = 0; i < compParms.infectionSize - 1; i++) {
        deltaInfection[i] = columnLoads[i + 1] - columnLoads[i];
    }
    testDeltaInfectionsAddToOne(*this);
}

State::State(IntegrationParams integrationParams, ModelParams modelParams) {
    size_t ageSize = modelParams.maxAge / integrationParams.deltaTime;
    size_t infectionSize = ageSize * integrationParams.infectionSizeToAgesSizeRatio;
    susceptibles = std::make_unique<Susceptibles>(ageSize);
    infecteds = std::make_unique<Infecteds>(ageSize, infectionSize);
    intParms = integrationParams;
    modParms = modelParams;
    initializeComputedParameters();
    setInitialInfecteds(*this);
    updateComputedParameters();
    verifySusceptibleTotal(*this, modParms.initialSusceptiblePop);
    verifyInfectedTotal(*this, modParms.initialInfectedPop);
}

void State::updateComputedParameters() {
    size_t maxAgeStep = modParms.maxAge/intParms.deltaTime;
    compParms.totalInfection = 0.0;
    compParms.infectedPop = 0.0;
    compParms.susceptiblePop = 0.0;
    compParms.transferRate = 0.0;
    compParms.aveLoad = 0.0;
    compParms.maxInfectedsPopDensity = 0.0;
    compParms.maxSusceptiblesPopDensity = 0.0;
    double dt = intParms.deltaTime;
    double it; //delta infection
    double totalLoad = 0.0;
    auto& columnLoads = *compParms.columnLoads;
    auto& susceptiblesVec = *susceptibles->getCurrentState();
    for (size_t xAge = 0; xAge < maxAgeStep; xAge++) {
        double popDensityAtAge = susceptiblesVec[xAge];
        compParms.susceptiblePop += popDensityAtAge * dt;
        if (popDensityAtAge > compParms.maxSusceptiblesPopDensity) {
            compParms.maxSusceptiblesPopDensity = popDensityAtAge;
        }
        for (size_t xLoad = 0; xLoad < compParms.infectionSize; xLoad++) {
            double popDensityAtLoad = infecteds->getIndex(xAge, xLoad);
            if (popDensityAtLoad > compParms.maxInfectedsPopDensity) {
                compParms.maxInfectedsPopDensity = popDensityAtLoad;
            }
            it = compParms.deltaInfectionForLoad->at(xLoad);
            double infectedPopAtLoad = popDensityAtLoad * dt * it;
            compParms.infectedPop += infectedPopAtLoad;
            compParms.transferRate += modParms.beta * columnLoads[xLoad] * popDensityAtLoad;
            totalLoad += columnLoads[xLoad] * infectedPopAtLoad;
        }
    }
    compParms.popSize = compParms.susceptiblePop + compParms.infectedPop;
    compParms.aveLoad = totalLoad / compParms.infectedPop;
}

// See https://en.wikipedia.org/wiki/Netpbm for a description of this simple format.
// To view the greymap on linux, use 'gimp <filename>', where mfilename must end
// in .pgm.
void State::writeInfectedsPGM(const std::string& filename) const {
    size_t rows = compParms.infectionSize;
    size_t columns = compParms.ageSize;
    FILE* file = fopen(filename.c_str(), "w");
    fprintf(file, "P2\n%lu %lu\n255\n", columns, rows);
    for (size_t xLoad = 0; xLoad < rows; xLoad++) {
        bool firstTime = true;
        for (size_t xAge = 0; xAge < columns; xAge++) {
            if (!firstTime) {
                fputc(' ', file);
            }
            firstTime = false;
            double popDensityAtLoad = infecteds->getIndex(xAge, rows - 1 - xLoad);
            fprintf(file, "%u", (uint32_t)(popDensityAtLoad * 255 / compParms.maxInfectedsPopDensity));
        }
        fputc('\n', file);
    }
    fclose(file);
}

namespace {

void setBit(bool* bitmap, size_t rows, size_t columns, size_t x, size_t y, uint32_t width) {
    assert(x < columns && y < rows);
    size_t left = x > width? x - width : 0;
    size_t bottom = y > width? y - width : 0;
    size_t right = x + width < columns? x + width : columns - 1;
    size_t top = y + width < rows? y + width : rows - 1;
    for (size_t xp = left; xp <= right; xp++) {
        for (size_t yp = bottom; yp <= top; yp++) {
            bitmap[yp * columns + xp] = true;
        }
    }
}

}  // namespace

// See https://en.wikipedia.org/wiki/Netpbm for a description of this simple format.
// To view the bitmap on linux, use 'gimp <filename>', where mfilename must end
// in .pbm.
void State::writeSusceptiblesPBM(const std::string& filename, uint32_t width) const {
    size_t columns = compParms.ageSize;
    size_t rows = columns / 3;
    bool* bitmap = new bool[rows * columns];
    memset(bitmap, '\0', rows * columns * sizeof(bool));
    std::vector<double>& dist = *susceptibles->getCurrentState();
    for (size_t xAge = 0; xAge < columns; xAge++) {
        double popAtAge = dist[xAge];
        size_t rowIndex = (rows - 1) * (1.0 - popAtAge / compParms.maxSusceptiblesPopDensity); 
        setBit(bitmap, rows, columns, xAge, rowIndex, width);
    }
    FILE* file = fopen(filename.c_str(), "w");
    fprintf(file, "P1\n%lu %lu\n", columns, rows);
    for (size_t xRow = 0; xRow < rows; xRow++) {
        for (size_t xCol = 0; xCol < columns; xCol++) {
            bool filled = bitmap[xRow * columns + xCol];
            char c = filled? '1' : '0';
            fputc(c, file);
        }
        fputc('\n', file);
    }
    fclose(file);
    delete[] bitmap;
}
