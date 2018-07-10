CC = g++
CCFLAGS = -Wextra -Wall -march=native -O3
LINKFLAGS = 
SOURCES = SparseAssembler.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = SparseAssembler

all: $(SOURCES) $(EXECUTABLE)

%.o:%.c
	$(CC) $(CCFLAGS) -c $<

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LINKFLAGS) $(OBJECTS) -o $@

clean:
	rm -f *.o core.*
