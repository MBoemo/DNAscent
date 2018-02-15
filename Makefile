OBJS = event_detection.o event_handling.o common.o build_model.o data_IO.o Osiris_train.o Osiris.o
CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS = -L Penthus/ -l Penthus -fopenmp
CXXFLAGS = -Wall -c -O2 -std=c++11 $(DEBUG)
CFLAGS = -Wall -c -std=c99 -O2 $(DEBUG)
LFLAGS = -Wall -O2 $(DEBUG)

MAIN_EXECUTABLE = bin/Osiris

$(MAIN_EXECUTABLE) : $(OBJS)
	$(CXX) $(LFLAGS) $(OBJS) -o bin/Osiris $(LIBFLAGS)

event_detection.o : src/scrappie/event_detection.h src/scrappie/event_detection.c
	$(CC) $(CFLAGS) src/scrappie/event_detection.c

event_handling.o : src/scrappie/event_detection.h src/common.h src/poreModels.h src/data_IO.h src/event_handling.h src/event_handling.cpp
	$(CXX) $(CXXFLAGS) src/event_handling.cpp $(LIBFLAGS)

common.o : src/error_handling.h src/common.h src/common.cpp
	$(CXX) $(CXXFLAGS) src/common.cpp $(LIBFLAGS)

data_IO.o : src/error_handling.h src/data_IO.h src/data_IO.cpp
	$(CXX) $(CXXFLAGS) src/data_IO.cpp $(LIBFLAGS)

build_model.o : src/poreSpecificParameters.h src/build_model.h src/build_model.cpp src/data_IO.h
	$(CXX) $(CXXFLAGS) src/build_model.cpp $(LIBFLAGS)

Osiris_train.o : src/common.h src/build_model.h src/data_IO.h src/error_handling.h src/event_handling.h src/poreModels.h src/Osiris_train.h src/Osiris_train.cpp
	$(CXX) $(CXXFLAGS) src/Osiris_train.cpp $(LIBFLAGS)

#Osiris_detect.o : src/common.h src/build_model.h src/data_IO.h src/error_handling.h src/event_handling.h src/Osiris_detect.h src/Osiris_detect.cpp
#	$(CXX) $(CXXFLAGS) src/Osiris_detect.cpp $(LIBFLAGS)

Osiris.o : src/Osiris.cpp src/Osiris_train.h src/Osiris_fixedPos.h src/data_IO.h src/build_model.h src/event_handling.h
	$(CXX) $(CXXFLAGS) src/Osiris.cpp $(LIBFLAGS)

clean:
	rm $(OBJS) $(MAIN_EXECUTABLE)
