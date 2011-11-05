
DEBUG=-g

CXX=/opt/local/bin/g++-mp-4.4
CXXFLAGS=$(DEBUG) -W -Wall -pedantic -O3 -funroll-loops
LDFLAGS=$(DEBUG)
LDLIBS=-L/opt/local/lib -lboost_graph 

SRCS=pathmeasure.cpp mersenne.cpp userintf.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

pathmeasure : pathmeasure.o mersenne.o userintf.o
	$(CXX) $(LDFLAGS) -o pathmeasure $(LDLIBS) $(OBJS)

pathmeasure.o : pathmeasure.cpp randomc.h

mersenne.o : mersenne.cpp randomc.h
userintf.o : userintf.cpp randomc.h

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) pathmeasure
