#
#	Makefile for program NAB example programs rfi_nab and rfi_prior_nab
#
#	These implement the receiver function problem. User specific
#	parts of this makefile occur between the dashed lines. Only
#	edit between the dashed lines.
#
#	type `make all' for a complete installation
#	type `make clean' to remove all object files and executables 
#
# setup NA macros for compilers and archiver

include makenab.macros   # Edit this file to change location or
                         # type of compilers (e.g. to gcc, g77)

#---- start user specific make --------------------------------------------

#USERPROGS = rfi_nab rfi_prior_nab
#USERPROGS = rfi_subs rfi_nab gauss_nab 
USERPROGS = freedrift_nab pecube_nab gauss_nab

DIR = .

#RFISUBS = $(DIR)/theo.o $(DIR)/qlayer.o \

pecube_nab: nab pecube_nab.f
	$(F77) $(FFLAGS) lib/nab.o \
        pecube_nab.f -o ../bin/pecube_nab

gauss_nab: nab gauss_nab.f
	$(F77) $(FFLAGS) lib/nab.o \
        gauss_nab.f -o ../bin/gauss_nab

freedrift_nab: nab freedrift_nab.f
	$(F77) $(FFLAGS) lib/nab.o \
        freedrift_nab.f -o ../bin/freedrift_nab

#rfi_prior_nab: nab rfi_prior_nab.f
#	$(F77) $(FFLAGS) lib/nab.o rfi_prior_nab.f \
#        -o ../bin/rfi_prior_nab

cleanuser :
	\rm -f utl/fit2ppd/rfi_subs/*.o
	\rm -f utl/fit2ppd/rfi_subs/rfi_syn/*.o
	\rm -f utl/fit2ppd/rfi_interface.o
	\rm -f ../data/rfi_files/NA_MDL/rfi_models
	\rm -f ../data/rfi_files/NA_SUM/*

lib: $(MACH)

sparc:
	cp  ../data/rfi_files/NA_SUM/rec_func_demo_sun_big_endian rec_func
dec:
	cp  ../data/rfi_files/NA_SUM/rec_func_demo_dec_lit_endian rec_func
sgi:
	cp  ../data/rfi_files/NA_SUM/rec_func_demo_dec_lit_endian rec_func
x86:
	cp  ../data/rfi_files/NA_SUM/rec_func_demo_dec_lit_endian rec_func
mac:
	cp  ../data/rfi_files/NA_SUM/rec_func_demo_dec_lit_endian rec_func


#---- end user specific make ----------------------------------------------
#
# do not change below
#
all:: utl nab $(USERPROGS) mes

nab	: NAB_src/nab.F
	$(F77) $(FFLAGS) -c NAB_src/nab.F 
	mv nab.o lib

utl::	utl/naplot/naplot.f
	cd utl/naplot; make all -f Makefile
	cd utl/fit2ppd; make all -f Makefile
	cd utl/combine; make all -f Makefile
	cd utl/nad; make all -f Makefile


clean:: cleanuser
	\rm -f ../bin/*
	\rm -f ../data/nab.out
	\rm -f ../data/nab.sum
	\rm -f ../data/nab.samples
	\rm -f ../data/sobol.coeff
	\rm -f ../data/model.opt
	\rm -f ../data/ensemble.nad
	\rm -f ../data/demo_results/model.opt
	\rm -f ../data/demo_results/junk
	\rm -f ../data/demo_results/chi_demo.nad
	\rm -f ../data/fit2ppd/junk
	\rm -f ../data/fit2ppd/misfit_demo.nad
	\rm -f ../data/fit2ppd/model.opt
	\rm -f ../data/naplot/*.matrix
	\rm -f ../data/naplot/*.ps
	\rm -f lib/*.o
	\rm -f lib/libnab.a
	\rm -f utl/naplot/Lib/*
	\rm -f utl/naplot/subs_src/*.o
	\rm -f utl/naplot/subs_src/plotlib_src/*.o
	\rm -f utl/fit2ppd/*.o
	\rm -f utl/fit2ppd/model.opt
	\rm -f utl/combine/*.o
	\rm -f utl/nad/*.o
#	cp ../data/nab.in_test ../data/nab.in 
#	cp ../data/pspace.in_rfi ../data/pspace.in 
#	cp ../data/gauss.in_24d ../data/gauss.in 

mes:
	echo 'Make all completed'
