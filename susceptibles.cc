#include "susceptibles.h"
#include "initialconditions.h"

#include <cassert>

void Susceptibles::setBirths(double births) {
  assert((*currentState_)[0] == 0.0);
  (*currentState_)[0] = births;
}
