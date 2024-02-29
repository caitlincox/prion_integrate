#include "birth.h"
#include "integrator.h"

double ConstantBirthScheme::calculateBirth() {
    return birthConstant;
}

//Integrate over age and infection load and return perportional birth rate
double BirthRateBirthScheme::calculateBirth() {
    
}

//Have births at time t equal deaths at time t. Integrate over age/infection load, and re-add the difference.
double ReplacementBirthScheme::calculateBirth() {}