# prion\_integrate

# TODO

* Implement Death::kill.
* Implement function to infect susceptibles at a given time step.
  See TODO in integrator.cc.
* Implement Death class.  Currently, animals only die throug hitting max age or infection >= 1.0.
* Watch every line execute in debugger.

# Notes

Curerntly, initial infecteds do not look at all like the steady-state.  The
steady-state should have lower infection loads for younger pops, and more for
older pops.  Also, initial infections is heavier in the young.  Given that a
real initial infection is likely to start with 1 animal, it really doesn't
matter how we distribute initial load.
