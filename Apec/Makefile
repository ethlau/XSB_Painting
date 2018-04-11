CC=gcc
CFLAGS = -O2 -openmp
#LIBS = -L/Users/ethlau/Research/program/cfitsio -lcfitsio
LIBS = -lcfitsio
INCLUDES = 

APEC_SRCS = Apec.c atomdb_make_spectrum.c calc_continuum.c calc_lines.c messages.c readapec.c read_continuum_data.c read_fits_spectrum.c read_line_data.c read_parameters.c gaussianLine.c 

APEC_OBJS = $(patsubst %.c,%.o,$(APEC_SRCS))

%.o: %.c messages.h readapec_spectrum.h 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
test : main.c $(APEC_OBJS)  messages.h  readapec.h  readapec_spectrum.h 
	$(CC) $(CFLAGS) $(INCLUDES) $(APEC_OBJS) main.c -o test $(LIBS)
clean:
	rm -f *.o test
