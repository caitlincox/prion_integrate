#ifndef INIT
#define INIT

#include <memory>
#include <vector>

class InitialConditions {
public:
    static double intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad);
    double weibullOfAge(double weibullLambda, double weibullKappa, double age);
    std::unique_ptr<std::vector<double>> startingDistribution(double weibullLambda, double weibullKappa, double maxAge, double deltaT,
                                                              double popSize);
    std::unique_ptr<std::vector<double>> initialInfecteds(double numInfecteds, double deltaT, double maxAge, double numInfectionBuckets,
                                                                         double scaleParam, double shapeParam, const std::vector<double>& loadVec,
                                                                          std::vector<double>& susceptibles, double numSusceptibles, double aveLoad,
                                                                          double c);
};

#endif
