
##
#  Check Inline::C {inc,libs} paths in ../bin/*.pl scripts, lines {101,102}.
#  No change is necessary for Strawberry Perl or for {gmp,mpir} installed
#  in /usr/{include,lib[64]}.
##

INCDIR = /usr/local/include
LIBDIR = /usr/local/lib

CC = gcc -std=gnu99 -Wno-attributes

ifeq "${CFLAGS}" ""
  CFLAGS = -O2 -DBIG_SIEVE=0
endif

ifeq "${OS}" ""
  OS = $(shell uname -s)
endif

ifeq "$(strip ${OS})" "Darwin"
  OPENMP =
else
  OPENMP = -fopenmp
endif

ifeq "$(strip ${OS})" "Linux"
  RPATH = -Wl,-rpath=${LIBDIR}
else
  RPATH =
endif

all: pi-perl

pi-perl:
	perl -MExtUtils::MakeMaker -MStorable -MTime::HiRes -e 1
	perl ../bin/pi-hobo.pl CFLAGS="${CFLAGS}"

pi-gmp:
	${CC} ${OPENMP} \
	  ${CFLAGS} -DUSE_GMP pgmp-chudnovsky.c \
	  -I${INCDIR} -L${LIBDIR} ${RPATH} \
	  -o ../bin/pi-gmp.exe -lgmp -lm

pi-mpir:
	${CC} ${OPENMP} \
	  ${CFLAGS} -DUSE_MPIR pgmp-chudnovsky.c \
	  -I${INCDIR} -L${LIBDIR} ${RPATH} \
	  -o ../bin/pi-mpir.exe -lmpir -lm

