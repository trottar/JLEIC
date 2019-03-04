#                     
#          R A P G A P
#

#CC_OPTS=-arch i386
#CC_OPTS=-arch x86_64
F77=gfortran
F77_OPTS= -fno-automatic -fno-second-underscore -falign-commons -fd-lines-as-comments
# F77_DEBUG= -g3 -ggdb3

OPTS_CC = -fconserve-space

CC_OPTS = $(OPTS_CC) $(PROC_OPTS) $(DEBUG_CC)

#CERNLIB=$(shell cernlib -lpdflib804 -ljetset74 -lariadne -lmathlib -lkernlib -lpacklib)
CERNLIB=$(shell cernlib pdflib804 jetset74 ariadne mathlib kernlib packlib)

clean:
	rm -f *.x *.o rapgap/*.o rapgap/*.a heracles/*.o heracles/*.a bases/*.o bases/*.a

all: rgEIC.x 

rgEIC.x:  rgEIC.o analys.o Makefile rapgap/librapgap.a heracles/libheracles.a bases/libbases.a
	$(F77) $(F77_OPTS) $(F77_DEBUG) \
		rgEIC.o analys.o -o rgEIC.x \
		-L./rapgap -lrapgap \
		-L./heracles -lheracles \
		-L./rapgap -lrapgap \
		-L./bases -lbases \
		$(CERNLIB)

rapgap/librapgap.a:  rapgap20806.f Makefile fsplit
	cd rapgap; \
		rm -f *.f; \
		rm -f *.o; \
		../fsplit ../rapgap20806.f; \
		rm -f arinit.f; \
		cp ../percut.inc .; \
		$(F77) $(F77_OPTS) $(F77_DEBUG) -c *.f
	rm -f rapgap/librapgap.a
	ar -r rapgap/librapgap.a rapgap/*.o

heracles/libheracles.a:  heracles461.f Makefile fsplit
	cd heracles; \
		../fsplit ../heracles461.f; \
		$(F77) $(F77_OPTS) $(F77_DEBUG) -c *.f
	rm -f heracles/libheracles.a
	ar -r heracles/libheracles.a heracles/*.o
	rm heracles/*.f
	rm heracles/*.o

bases/libbases.a:  bases51.f Makefile fsplit
	cd bases; \
		../fsplit ../bases51.f; \
		$(F77) $(F77_OPTS) $(F77_DEBUG) -c *.f
	rm -f bases/libbases.a
	ar -r bases/libbases.a bases/*.o
	rm bases/*.f
	rm bases/*.o

.f.o:
	$(F77) $(F77_OPTS) $(F77_DEBUG) -c $*.f  

fsplit: fsplit.c
