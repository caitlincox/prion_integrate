# prion_integrate

# TODO

* Implement function to infect susceptibles at a given time step.
* Output data in a manner that can be viewed easily.
* Lotss of debugging.
* What is the Deaths class for?  We kill by deleting buckets that are too old
  or too infected.
* Fix gammaVal computation such that the totals of all infected after initial
  state is equal to initialInfecteds.
* Is aveLoad over the entire pop or just infecteds?  Fix in State::updateComputedParameters
* Watch every line execute in debugger.
