#
#	Make file for program NA.f
#

# Ways to use this make file 
#    > `make all' 
#    > `make clean'
#    
# For a more detail explanation see HTML manual
#
# Setup NA macros for compilers and archiver

include makena.macros	# Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)

#MPI = -DNA_MPI=1           # MPI option (original parallel NA only used if TYPE=mpi)
MPI = -DNA_MPIP=1           # MPI option (continuous parallel NA only used if TYPE=mpi)

#---- start user specific make  (Edit below) ----------------------------------

# edit the three variables below to compile user routines and main program
# The present entries refer to the example program for fitting 
# Seismic receiver functions. (See HTML manual) 

                         # Enter here the directory name containing all 
			 # user code other than main program
			 # Note the users Makefile must exist in this 
			 # directory to compile users routines (see example). 
USER_DIR = freedrift_subs
                         # Enter here the list of object files or libraries 
			 # to be linked in with user code
USER_OBJ = \
         $(USER_DIR)/triple_misfit.o #\
#         $(USER_DIR)/calcmisfit.o \
#         $(USER_DIR)/display_model.o \
#         $(USER_DIR)/display_final.o \
#         $(USER_DIR)/lofw.o \
#         $(USER_DIR)/output_summary.o \
#         $(USER_DIR)/read_vmodelNA.o

                         # Enter here name of user's main routine
			 # Note that the main program name is
			 # assumed to be `USER_MAIN.F' 
USER_MAIN = freedrift_na

user :
	 cd $(USER_DIR); make -f Makefile

all:: $(TYPE)

serial:: user na user_main mes

mpi:: user na_mpi user_main_mpi mes_mpi

user_main:
	$(F77) $(FFLAGS) -o ../bin/$(USER_MAIN) $(USER_MAIN).F \
	lib/na.o lib/utils.o $(USER_OBJ)

user_main_mpi:
	$(MF77) $(MPI) $(FFLAGS) -o ../bin/$(USER_MAIN) $(USER_MAIN).F \
	lib/na.o lib/utils.o #$(USER_OBJ)

# remove user object files and libraries

cleanuser:
	cd $(USER_DIR); \rm -fr *.o
	\rm -f ../bin/$(USER_MAIN)

#---- end user specific make ----Do not edit below-------------------------
#
#
na:     
	$(CC) -c NA_src/utils.c
	mv utils.o lib
	$(F77) $(FFLAGS) -c $(F77_COMPAT) NA_src/na.F
	mv na.o lib

na_mpi:     
	$(CC) -c NA_src/utils.c
	mv utils.o lib
	$(MF77) $(MPI) $(FFLAGS) $(F77_COMPAT) -c NA_src/na.F
	mv na.o lib

# remove all object files and libraries 

distclean: clean

clean: cleanuser
	-rm -f ../bin/*
	-rm -f lib/*

#remove_demo:
#	\rm -fr ../data/demo_results
#	\rm -fr ./freedrift_subs
#	\rm -fr ../data/freedrift_plot
#	\rm -fr ../data/freedrift_files
#	\rm -fr ../data/freedrift_na

mes:
	cat < utl/message 

mes_mpi:
	cat < utl/message_mpi 
