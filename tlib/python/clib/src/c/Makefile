#---------------------------#
DO_MPI= no
#DO_MPI= yes


ifeq ($(DO_MPI), yes)
  CC= mpicc
else
  CC= gcc
endif
MPICC= mpicc
FC= ifort
#---------------------------#



# C main() call fortran subroutine, for ifort.
COMFFLAGS = -nofor_main    
FLIB = -L/usr/local/lib  -lm   

#-----------------------------------------------------------------#
#------------------------- My NOTEBOOKS---------------------------#
#-----------------------------------------------------------------#

ifeq ($(MACHINE),mylaptop)
    MYLIBPATH= /Users/harrygun/workspace/lib/mylib
    INIPATH = /Users/harrygun/workspace/lib/iniparser
    MPIPATH = /Users/harrygun/software/work/mpich2/mpich2-1.2.1p1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/Users/harrygun/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -arch x86_64 -D_MPI_ -I /Users/harrygun/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -arch i386 -arch x86_64 -I /Users/harrygun/workspace/lib/gsl/include \
		 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)

	#CFLAGS = -g -I /Users/harrygun/workspace/lib/gsl/include \
	#	 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)
    endif

endif

#---- gwln ----
ifeq ($(MACHINE),gwln)
    MYLIBPATH= /home/wangxin/workspace/lib/mylib
    INIPATH = /home/wangxin/workspace/lib/iniparser/$(HOSTNAME)
    MPIPATH = /home/wangxin/software/mpich2/gwln1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)/$(HOSTNAME)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -D_MPI_ -I /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -fPIC -I /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    endif

endif

#---- HHPC ----
ifeq ($(MACHINE),HHPC)
    MYLIBPATH= /home/wangxin/lib/mylib
    INIPATH = /home/wangxin/lib/iniparser
    MPIPATH = /usr/mpi/intel/openmpi-1.4.3

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 


    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -fPIC -D_MPI_ -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -fPIC -I/home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    endif

endif







#################################################################
#MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 


LIB  = $(INILIB) $(GSLIB) $(MYLIB)
#lib  = $(inilib) $(gslib) 




.SUFFIXES : .o .c
.c.o :
	$(CC) $(CFLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<





# workspace path
WORKSPACE = ../../workspace
AR=ar
ARFLAGS=-r
RANLIB=ranlib


default: libmysf.a



libmysf.a:  utilis.o wigner.o alm.o
ifeq ($(MACHINE),mylaptop)
	@($(AR) $(ARFLAGS) $@ $^)
	@($(RANLIB) $@)
endif	

ifeq ($(MACHINE),gwln)
	ar -r  $@ $?
endif

ifeq ($(MACHINE),HHPC)
	ar -r  $@ $?
endif

	cp $@ ../../lib/

#mpi_test: mpi_test.o
#	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)


#--------------------------------------------------#
clean:
	rm -f *.o *.c~ *.a


