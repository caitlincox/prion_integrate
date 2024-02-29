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
initialcondition.cc \
integrator.cc \
susceptibles.cc \
testing.cc

prion_integrate: $(HDRS) $(SRCS)
	gcc -g -o prion_integrate $(SRCS) -lm
