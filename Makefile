HDRS= \
birth.h \
death.h \
infecteds.h \
initialconditions.h \
integrator.h \
state.h \
susceptibles.h \
tests.h

SRCS= \
birth.cc \
death.cc \
infecteds.cc \
initialconditions.cc \
integrator.cc \
prion_integrate.cc \
state.cc \
susceptibles.cc \
tests.cc

prion_integrate: $(HDRS) $(SRCS)
	g++ -Wall -g -o prion_integrate $(SRCS) -lm
	mkdir -p ./data

clean:
	rm -f prion_integrate
