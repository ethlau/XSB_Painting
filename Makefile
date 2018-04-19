CC = icc
CXX= icpc

CFLAGS = -O2 
CXXFLAGS = -O2 -std=c++11
#CXXFLAGS = -O0 -g
#CXXFLAGS += -DMPI_PARALLEL
#CXXFLAGS += -DOUTPUT_NFW_DENSITY
#CXXFLAGS += -DLONGIDS
#CXXFLAGS += -DEXCLUDE_SUBHALOS
#CXXFLAGS += -DMASS_SELECTION=1e14
#CXXFLAGS += -DROCKSTAR_CONCENTRATION
FLINE = -D'LINE_FITS_FILE="/home/fas/nagai/etl28/programs/Xrays/atomdb/atomdb_v3.0.9/apec_line.fits"'
FCOCO = -D'COCO_FITS_FILE="/home/fas/nagai/etl28/programs/Xrays/atomdb/atomdb_v3.0.9/apec_coco.fits"'

CFLAGS += -I/home/fas/nagai/etl28/programs/cfitsio
CXXFLAGS += -I/home/fas/nagai/etl28/programs/cfitsio

CLIBS = -lgsl -lgslcblas -lm -L/home/fas/nagai/etl28/programs/cfitsio -lcfitsio 

APEC_SRCS = Apec.c atomdb_make_spectrum.c calc_continuum.c calc_lines.c messages.c readapec.c read_continuum_data.c read_fits_spectrum.c read_line_data.c read_parameters.c gaussianLine.c

APEC_OBJS = $(patsubst %.c,Apec/%.o,$(APEC_SRCS))

IO_SRCS = read_halo.cpp read_rockstar.cpp read_binary.cpp

IO_OBJS = $(patsubst %.cpp,io/%.o,$(IO_SRCS))

Apec/%.o: Apec/%.c Apec/%.h
	$(CC) $(CFLAGS) $(FCOCO) $(FLINE) -I. -I./Apec -c $< -o $@

io/%.o: io/%.c io/%.h
	$(CXX) $(CXXFLAGS) -I. -I./io -c $< -o $@

sbprof: main.cpp gas_model.o xray.o ConfigParser.o $(APEC_OBJS) $(IO_OBJS)
	$(CXX) $(CXXFLAGS) -o sbprof main.cpp gas_model.o xray.o ConfigParser.o $(IO_OBJS) $(APEC_OBJS) $(CLIBS) 

xray.o: xray.c xray.h
	$(CC) $(CFLAGS) -c xray.c

gas_model.o: gas_model.cpp gas_model.h
	$(CXX) $(CXXFLAGS) -c gas_model.cpp


#read_halo: io/read_lightcone.cpp io/read_halo_rs.cpp io/read_halo_simple.cpp io/read_halo.h
#	$(CXX) $(CXXFLAGS) -c io/read_lightcone.cpp io/read_halo_rs.cpp io/read_halo_simple.cpp

ConfigParser.o: ConfigParser/ConfigParser.c ConfigParser/ConfigParser.h
	$(CC) -c ConfigParser/ConfigParser.c

clean:
	/bin/rm -f *.o io/*.o Apec/*.o sbprof
