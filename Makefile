CC = g++
CFLAGS =

OBJECTS = spectrum2rgb.o

SPECTRUM2RGB = spectrum2rgb.c spectrum2rgb.h

default: hydro

starcolours: $(OBJECTS)
	$(CC) $(CFLAGS) -o spectrum2rgb $(OBJECTS) -lm -O3 -g

#-O3: optimize this thing

spectrum2rgb.o: $(SPECTRUM2RGB)
	$(CC) $(CFLAGS) -c spectrum2rgb.c

clean:
	rm -f $(OBJECTS)	 
