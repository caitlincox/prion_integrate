#include "birth.h"
#include "integrator.h"

//Stupidest possible thing. Useful for debugging/checking intuition. Constant birth not perportional to population.
class ConstantBirthScheme : public BirthScheme {
 public:
    ConstantBirthScheme(double birthConstant) {
        birthConstant_ = birthConstant;
    }

    double calculateBirth(const State& state) override {
        return birthConstant_;
    };

 private:
    double birthConstant_;
};

std::unique_ptr<BirthScheme> newConstantBirthScheme(double birthConstant) {
    return std::make_unique<ConstantBirthScheme>(birthConstant);
}

//SIR-inspired constant population. Sets birth equal to death. Was used in Stringer et al.
class ReplacementBirthScheme : public BirthScheme {

public:
    //We set the birth pt such as to force our left Riemann sum integration to preserve the total.
    double calculateBirth(const State& state) override {
        return (state.compParms.ageDeaths + state.compParms.infectionDeaths + state.compParms.naturalDeaths) / state.intParms.deltaTime;
    }
};

std::unique_ptr<BirthScheme> newReplacementBirthScheme() {
    return std::make_unique<ReplacementBirthScheme>();
}

//Biologically more accurate birth, but we don't change the birthrate by age.
class BirthRateBirthScheme : public BirthScheme {

public:
    double calculateBirth(const State& state) override {
      return state.compParms.ageDeaths + state.compParms.infectionDeaths +
             state.compParms.naturalDeaths;;
    }

    //getters and setters
    double getBirthRate() {return birthRate;}
    void setBirthRate(double value) {birthRate = value;}
    
private:
    double birthRate;
};
