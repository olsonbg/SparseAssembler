CC=g++
OPTFLAGS=-march=native -O3
#EXTRAFLAGS=-Wextra -Wall
CPPFLAGS=$(OPTFLAGS) $(EXTRAFLAGS) -DUSE_ZLIB -DUSE_BZIP2 -DUSE_LZMA $(BOOST_INCLUDE_OPTS)
LINKFLAGS=-lz -lbz2 -llzma -lboost_iostreams $(BOOST_LINK_OPTS)
SOURCES=ReadFirstLine.cpp OpenFile.cpp MagicNumber.cpp SparseAssembler.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=SparseAssembler

all: $(SOURCES) $(EXECUTABLE)

%.o:%.cpp
	$(CC) $(CPPFLAGS) -c $<

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LINKFLAGS) $(OBJECTS) -o $@

clean:
	rm -f *.o core.*
