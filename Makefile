CXX	= g++
CXXFLAGS= -O3 -std=gnu++11 #-g -O3 -Wall -Wextra
#INC                = -I/usr/local/lib/
PROGRAM	= getN
OBJS	= obj.o #can change this
LIBS	= -lm #-L/usr/local/lib/



$(PROGRAM): $(OBJS)
	$(CXX) $(OBJS) -o $(PROGRAM)

obj.o: utils.h
	$(CXX) utils.h -o obj.o

clean:
		/bin/rm -f $(OBJS)
		/bin/rm -f $(PROGRAM)
		/bin/rm -f core
