SHELL = /bin/sh

FC = gfortran
FCF =  -frecord-marker=4 -O4 -m32 -g -fno-second-underscore -funroll-loops
FCL = -static -m32 -O4 -L.
FCLD = -m32 -O4 -L.
FCLMAP = -Wl,--print-map -static -m32 -O4 -L.
CRLIBM = 
CERNLIB = 
gfortran_lib = 
BOINC = 
NAGLIBC = 
BEAMGAS = 
SIXTRACK = 
LIBS = 
MODS = 
MOD_OBJS = 

all: 

track.o: track.f $(MODS)
	$(FC) -c $(FCF) $(FCFF) $<

sixve.o: sixve.f $(MODS)
	$(FC) -c $(FCF) $(FCFF) $<

dafor: dafor.o
	$(FC) $(FCLD) -o dafor dafor.o

astuce: astuce.o
	$(FC) $(FCLD) -o astuce astuce.o $(NAGLIBC)

sixvefox.o: sixvefox.f
	$(FC) -c $(FCF) $^

myboinc.o: myboinc.f
	$(FC) -c $(FCF) myboinc.f

# Default rules

%.o: %.f
	$(FC) -c $(FCF) $<

%.o %.mod: %.f90
	$(FC) -Wall -c $(FCF) $<

# Main program

: $(MYHDF5LIB) track.o sixve.o sixvefox.o dabnews.o lielib.o $(MOD_OBJS) $(BEAMGAS) $(BOINC) $(CRLIBM)
	$(FC) $(FCLMAP) -o  track.o sixve.o sixvefox.o dabnews.o lielib.o $(BEAMGAS) \
	$(BOINC) $(CERNLIB) $(CRLIBM) $(MOD_OBJS) $(gfortran_lib) $(LIBS) > sixtrack.map

: sixda.o sixsc.o sixscfox.o dabnew.o lielib.o $(MOD_OBJS) $(CRLIBM)
	$(FC) $(FCLMAP) -o  sixda.o sixsc.o sixscfox.o dabnew.o lielib.o \
	$(NAGLIB) $(NAGLIBC) $(CRLIBM) $(MOD_OBJS) $(gfortran_lib) $(LIBS) > sixtrack_da.map

clean: 
	rm -f *.ast *.c *.o *.mod astuce dafor  
