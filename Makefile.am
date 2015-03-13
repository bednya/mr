##
## Makefile.am for SM3r2m
##
## written by Andrey Pikelner
##
## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

SUBDIRS =         \
	tsil-1.21 \
	mr        \
	tests


if HAVE_MPREP
SUBDIRS+=math
endif

SUBDIRS+=examples  

EXTRA_DIST = tsil-1.21


