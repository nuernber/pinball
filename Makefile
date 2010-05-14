
LIBS = -lglut -lGLU -lGL 
CFLAGS = -Wall

ALL = pinball Timer
SOURCE = $(addsuffix .cpp,$(ALL))
OBJECTS = $(addsuffix .o,$(ALL))

.PHONY: clean

pinball: $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJECTS): Timer.h $(SOURCE)

$(SOURCE):
	$(CXX) -c $(CFLAGS) $@

clean:
	$(RM) *.o $(ALL)
