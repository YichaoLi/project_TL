ifeq ($(MACHINE),mylaptop)
    cflag= "-arch i386 -arch x86_64 -I /Users/harrygun/workspace/software/works/cuba/include -I/Users/harrygun/workspace/lib/gsl/include -I../../clib/src/c -I /Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include"
    libflag= "-L /Users/harrygun/workspace/software/works/cuba/lib -L /Users/harrygun/workspace/lib/mylib/script/clib -L ./  -L/Users/harrygun/workspace/lib/gsl/lib -lgsl -lgslcblas -lm -L../../clib/lib -lmysf"

endif

ifeq ($(MACHINE),gwln)
    cflag="-I /home/wangxin/workspace/software/cuba/include -I /home/wangxin/software/mpich2/gwln1/include -I/home/wangxin/workspace/lib/include  -I ../../clib/src/c"
    libflag="-L /home/wangxin/workspace/software/cuba/lib -L ./ -L/home/wangxin/workspace/lib/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"

endif


ifeq ($(MACHINE),HHPC)
	cflag="-I/home/wangxin/software/cuba/include -I/usr/mpi/intel/openmpi-1.4.3/include -I/home/wangxin/workspace/lib/include -I/home/wangxin/lib/gsl/include -I ../../clib/src/c"
	libflag="-L/home/wangxin/software/cuba/lib -L ./ -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"
endif


all:
	CFLAGS=$(cflag) LDFLAGS=$(libflag) \
	python setup.py build_ext --inplace

clean:
	@echo Cleaning
	@rm -f *.c *.o *.so *~
	@rm -rf build
