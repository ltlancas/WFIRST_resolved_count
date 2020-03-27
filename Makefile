CXX                = /usr/local/gfortran/bin/gcc
CXXFLAGS    = -fopenmp -O2 -std=c++11 #-g -O3 -Wall -Wextra
INC                = -I/usr/local/lib/
PROGRAM        = getN
OBJS 					 = obj.o #can change this
LIBS           =  -lstdc++ -lpthread -lfftw3_threads -lfftw3 -lm -L/usr/local/lib/

$(PROGRAM): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS)  -o $(PROGRAM) $(LIBS) $(INC)

obj.o: main.cpp # put your file name here,if you change OBJS above change obj.o
	$(CXX)  -c main.cpp -o obj.o

.PHONY:
clean:
		/bin/rm -f $(OBJS)
		/bin/rm -f $(PROGRAM)
		/bin/rm -f core
