SHELL = /bin/sh

FC = gfortran
FCF =  -frecord-marker=4 -O4 -m32 -g -fno-second-underscore -funroll-loops
FCL = -static -m32 -O4 -L.
FCLD = -m32 -O4 -L.
FCLMAP = -Wl,--print-map -static -m32 -O4 -L.
CRLIBM = 
CERNLIB = -L../lib32 -lgraflib -lgrafX11 -lpacklib_noshift -lkernlib_noshift -lg2c -lX11 -lxcb -lXau -lXdmcp -ldl -lpthread -lgcc_eh
gfortran_lib = 
BOINC = 
NAGLIBC = 
BEAMGAS = 
SIXTRACK = SixTrack_4518_cernlib_coll_gfortran_O4
LIBS = 
MODS = 
MOD_OBJS = 

all: SixTrack_4518_cernlib_coll_gfortran_O4

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

SixTrack_4518_cernlib_coll_gfortran_O4: $(MYHDF5LIB) track.o sixve.o sixvefox.o dabnews.o lielib.o $(MOD_OBJS) $(BEAMGAS) $(BOINC) $(CRLIBM)
	$(FC) $(FCLMAP) -o SixTrack_4518_cernlib_coll_gfortran_O4 track.o sixve.o sixvefox.o dabnews.o lielib.o $(BEAMGAS) \
	$(BOINC) $(CERNLIB) $(CRLIBM) $(MOD_OBJS) $(gfortran_lib) $(LIBS) > sixtrack.map

: sixda.o sixsc.o sixscfox.o dabnew.o lielib.o $(MOD_OBJS) $(CRLIBM)
	$(FC) $(FCLMAP) -o  sixda.o sixsc.o sixscfox.o dabnew.o lielib.o \
	$(NAGLIB) $(NAGLIBC) $(CRLIBM) $(MOD_OBJS) $(gfortran_lib) $(LIBS) > sixtrack_da.map

clean: 
	rm -f *.ast *.c *.o *.mod astuce dafor SixTrack_4518_cernlib_coll_gfortran_O4 
