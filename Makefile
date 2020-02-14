SOURCES = main.cpp

OBJECTS = main.o

filter: main.cpp
	g++ -o filter main.cpp

CXXFLAGS = -std=c++11 -D_GNU_SOURCE -Wall
