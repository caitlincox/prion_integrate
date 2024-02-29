#ifndef DEATH
#define DEATH
class Death{
public:
    double weibullDeathRate(double age);
    double findLambda(double aveLifespan, double kappa);
    Death(double kappa, double aveLifespan); //constructor

    //getters and setters
    double getSurvivalKappa(){return weibullSurvivorshipKappa;}
    double getSurvivalLambda(){return weibullSurvivorshipLambda;}
    void setSurvivalKappa(double kappa){weibullSurvivorshipKappa = kappa;}
    void setSurvivalLambda(double lambda){weibullSurvivorshipLambda = lambda;}
private:
    double weibullSurvivorshipLambda;
    double weibullSurvivorshipKappa;
};
#endif