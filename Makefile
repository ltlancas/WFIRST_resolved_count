CXX	= g++
CXXFLAGS= -O3 -std=gnu++11 #-g -O3 -Wall -Wextra
PROGRAM	= getN
DEPS	= utils.h

$(PROGRAM): main.o utils.o
	$(CXX) $(CXXFLAGS) -o $(PROGRAM) main.o utils.o

main.o: main.cpp utils.h
	$(CXX) $(CXXFLAGS) -c main.cpp

utils.o: 
	$(CXX) -c utils.cpp

clean:
		/bin/rm -f *.o
		/bin/rm -f $(PROGRAM)
		/bin/rm -f core
