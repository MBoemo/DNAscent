OBJS = event_detection.o
CXX = gcc
DEBUG = -g
CXXFLAGS = -Wall -c -O2 -std=c++11 -fopenmp $(DEBUG)
MAIN_LIBRARY = libEvents.a

$(MAIN_LIBRARY) : $(OBJS)
	ar rcs libEvents.a $(OBJS)

event_detection.o : event_detection.h
	$(CXX) $(CXXFLAGS) event_detection.c

clean:
	rm $(OBJS) $(MAIN_LIBRARY)
