#
# Makefile for nonspu program:
#       libcfits.a
#
# Oct-96 : original version by 
#
#       JDD/WDP
#       NASA GSFC
#       Oct 1996
#
# 25-Jan-01 : removed conditional drvrsmem.c compilation because this
#             is now handled within the source file itself.
# 09-Mar-98 : modified to conditionally compile drvrsmem.c. Also
# changes to target all (deleted clean), added DEFS, LIBS, added
# DEFS to .c.o, added SOURCES_SHMEM and MY_SHMEM, expanded getcol*
# and putcol* in SOURCES, modified OBJECTS, mv changed to /bin/mv
# (to bypass aliasing), cp changed to /bin/cp, add smem and
# testprog targets. See also changes and comments in configure.in
#

prefix		= /Users/wdeanpesnell/_files/stars/evolve
exec_prefix	= ${prefix}
DESTDIR		= 
CFITSIO_PREFIX	= $(prefix)
CFITSIO_LIB	= $(DESTDIR)$(exec_prefix)/lib
CFITSIO_INCLUDE	= $(DESTDIR)$(prefix)/include
INSTALL_DIRS	= $(DESTDIR)${prefix} ${CFITSIO_LIB} ${CFITSIO_LIB}/pkgconfig ${CFITSIO_INCLUDE}


SHELL =		/bin/sh
RANLIB =	ranlib
CC =		cc
CFLAGS =	-g -O2 #-Dg77Fortran -fPIC -fno-common
#FC =		g77
FC =		gfortran
LDFLAGS =	$(CFLAGS)
DEFS =		-DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MATH_H=1 -DHAVE_LIMITS_H=1 -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DHAVE_FTRUNCATE=1 -DHAVE_LONGLONG=1 -DHAVE_SHMEM_SERVICES=1 -DHAVE_UNION_SEMUN=1 -DHAVE_NET_SERVICES=1 
LIBS =		
FLEX =		flex
BISON =		bison

SHLIB_LD =	cc -dynamiclib
SHLIB_SUFFIX =	.dylib

#.SUFFIXES:	.for

.c.o:
		$(CC) -c $(CFLAGS) $(DEFS) $<

#.for.o:
#		$(FC) -c -o $@ $<
%.o : %.for
		$(FC) -c $(CFLAGS) -o $@ $<
#		$(FC) -c $(CFLAGS) $(DEFS) $<

SOURCES_FOR = 	nonspum.for cmplna.for matset.for matrixes.for pekeris.for \
                nonmisc.for lnanon.for pltdmp.for
		
UTIL_SOURCES = ../../utility/comopen.for ../../utility/newfil.for ../../recipes/locate.for

SOURCES = ${SOURCES_FOR} #${UTIL_SOURCES}

%OBJECTS = 	${SOURCES:.c=.o}

#CORE_OBJECTS = 	${CORE_SOURCES:.c=.o, .for=.o}

OBJECTS_FOR = 	${SOURCES_FOR:.for=.o} # ${UTIL_SOURCES:.for=.o} 

FITSIO_SRC =	f77_wrap1.c f77_wrap2.c f77_wrap3.c f77_wrap4.c

# ============ description of all targets =============
#       -  <<-- ignore error code

all:
		@if [ "x${FC}" = x ]; then \
			${MAKE} all-nofitsio; \
		else \
			${MAKE} stand_alone; \
		fi

nonspu.exe: ${OBJECTS_FOR}
		${FC} -o nonspu.exe $(OBJECTS_FOR) ${LIBS}

analyze_s_model:	analyze_s_model.for summary.for
		${FC} -o analyze_s_model.exe analyze_s_model.for summary.for pltdmp.for
all-nofitsio:
		${MAKE} stand_alone "FITSIO_SRC="

stand_alone:	libcfitsio.a

libcfitsio.a:	${OBJECTS}
		ar rv libcfitsio.a ${OBJECTS}; \
		${RANLIB} libcfitsio.a;

shared: libcfitsio${SHLIB_SUFFIX}

libcfitsio${SHLIB_SUFFIX}: ${OBJECTS}
		${SHLIB_LD} -o $@ ${OBJECTS}

install:	libcfitsio.a $(INSTALL_DIRS)
		@if [ -f libcfitsio.a ]; then \
			/bin/mv libcfitsio.a ${CFITSIO_LIB}; \
		fi; \
		if [ -f libcfitsio${SHLIB_SUFFIX} ]; then \
			/bin/mv libcfitsio${SHLIB_SUFFIX} ${CFITSIO_LIB}; \
		fi; \
		/bin/cp fitsio.h fitsio2.h longnam.h drvrsmem.h ${CFITSIO_INCLUDE}/; \
		/bin/cp cfitsio.pc ${CFITSIO_LIB}/pkgconfig

smem:		smem.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o smem smem.o -L. -lcfitsio -lm

testprog:	testprog.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o testprog testprog.o -L. -lcfitsio -lm ${LIBS}

fpack:		fpack.o fpackutil.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o fpack fpack.o fpackutil.o -L. -lcfitsio -lm ${LIBS}

funpack:	funpack.o fpackutil.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o funpack funpack.o fpackutil.o -L. -lcfitsio -lm ${LIBS}

fitscopy:	fitscopy.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o fitscopy fitscopy.o -L. -lcfitsio -lm ${LIBS}

speed:		speed.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o speed speed.o -L. -lcfitsio -lm ${LIBS}

imcopy:		imcopy.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o imcopy imcopy.o -L. -lcfitsio -lm ${LIBS}

listhead:	listhead.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o listhead listhead.o -L. -lcfitsio -lm ${LIBS}

cookbook:	cookbook.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o cookbook cookbook.o -L. -lcfitsio -lm ${LIBS}


clean:
	-	/bin/rm -f *.o libcfitsio.a libcfitsio${SHLIB_SUFFIX} \
			smem testprog y.output

distclean:	clean
	-	/bin/rm -f Makefile cfitsio.pc config.* configure.lineno

# Make target which outputs the list of the .o contained in the cfitsio lib
# usefull to build a single big shared library containing Tcl/Tk and other
# extensions.  used for the Tcl Plugin. 

listobjs:
	@echo ${OBJECTS_FOR}

listsource:
	@echo ${SOURCES_FOR}

# This target actually builds the objects needed for the lib in the above
# case
objs: ${OBJECTS_FOR}

$(INSTALL_DIRS):
	@if [ ! -d $@ ]; then mkdir -p $@; fi
