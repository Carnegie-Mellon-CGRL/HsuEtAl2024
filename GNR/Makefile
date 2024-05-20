CC = g++ -std=c++11
CFLAGS = 
LDFLAGS=
LDLIBS = -lgsl -lgslcblas -lm -lboost_program_options -target x86_64-apple-darwin20.3.0
SOURCES= vessel.cpp functions.cpp main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=gnr

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -c -o $@ $(LDLIBS)

clean:
	rm -f *.o *.mod *~ $(EXECUTABLE)

