#
#	Make discription file for Receiver function example subroutines
#

#FFLAGS = -O4
#F77 = /opt/SUNWspro/bin/f77

DEMO = ../..

include $(DEMO)/src/makena.macros   # Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)

DIR = .
FDRT = \
	$(DIR)/geo_misfit.o  \
        $(DIR)/fb_misfit.o


# make freedrift_subs 

freedrift_subs: $(FDRT) 

.f.o:
	$(F77) $(FFLAGS) -c $*.f
.c.o:
	$(CC) $(CFLAGS) -c $*.c
