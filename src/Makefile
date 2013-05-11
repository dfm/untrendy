CC ?= gcc
CLIB =
CFLAGS =

FF ?= gfortran
FFLAGS =

FPDIR = dierckx
FITPACK = $(FPDIR)/fpback.o $(FPDIR)/fpbspl.o $(FPDIR)/fpcurf.o \
		  $(FPDIR)/fpdisc.o $(FPDIR)/fpgivs.o $(FPDIR)/fpknot.o \
		  $(FPDIR)/fprati.o $(FPDIR)/fprota.o

.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c

.f.o:
	$(FF) $(FFLAGS) -o $*.o -c $*.f

spline_test: $(FITPACK) spline_test.o
	$(CC) $(CFLAGS) $(CLIB) spline_test.o $(FITPACK) -o spline_test

clean:
	rm -rf $(FITPACK) spline_test.o
