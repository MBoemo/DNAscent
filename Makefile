OBJS = build_model.o data_IO.o Osiris_train.o Osiris.o
CC = g++
DEBUG = -g
PENFLAGS = -L Penthus/ -lPenthus
CFLAGS = -Wall -c -std=c++11 $(DEBUG)
LFLAGS = -Wall $(DEBUG)

bin/Osiris : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o bin/Osiris $(PENFLAGS)

data_IO.o : src/data_IO.h src/data_IO.cpp
	$(CC) $(CFLAGS) src/data_IO.cpp

build_model.o : src/build_model.h src/build_model.cpp src/data_IO.h
	$(CC) $(CFLAGS) src/build_model.cpp $(PENFLAGS)

Osiris_train.o : src/Osiris_train.h src/Osiris_train.cpp src/data_IO.h src/build_model.h src/utility.h
	$(CC) $(CFLAGS) src/Osiris_train.cpp $(PENFLAGS)

Osiris.o : src/Osiris.cpp src/Osiris_train.h src/data_IO.h src/build_model.h
	$(CC) $(CFLAGS) src/Osiris.cpp $(PENFLAGS)

clean:
	rm $(OBJS)
