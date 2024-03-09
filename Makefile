HDRS= \
birth.h \
death.h \
infecteds.h \
initialconditions.h \
integrator.h \
state.h \
susceptibles.h

SRCS= \
birth.cc \
death.cc \
infecteds.cc \
initialconditions.cc \
integrator.cc \
prion_integrate.cc \
susceptibles.cc

prion_integrate: $(HDRS) $(SRCS)
	g++ -g -o prion_integrate $(SRCS) -lm

clean:
	rm -f prion_integrate
