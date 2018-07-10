CC=g++
CPPFLAGS=-Wextra -Wall -march=native -O3
SOURCES=SparseAssembler.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=SparseAssembler

all: $(SOURCES) $(EXECUTABLE)

%.o:%.cpp
	$(CC) $(CPPFLAGS) -c $<

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LINKFLAGS) $(OBJECTS) -o $@

clean:
	rm -f *.o core.*
