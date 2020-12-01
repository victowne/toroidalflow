SRCS =	gem_com.f90 gem_equil.f90 gem_main.f90 gem_outd.f90 gem_fcnt.f90 gem_fft_wrapper.f90 gem_gkps_adi.f90

OBJS =	gem_com.o gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_fft_wrapper.o gem_gkps_adi.o

LIBS = $(DFFTPACK) -mkl
PLIB = gem_pputil.o

#F90 = ftn
F90 = mpif90
#OPT = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
OPT = -FR -r8 -O3
LDFLAGS = 
#INCLUDES = -I/share/apps/soft/intel2019/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include

#all : gem

gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o
	$(F90) -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) 

gem_pputil.o: gem_pputil.f90
	$(F90) -c $(OPT) $(INCLUDES) gem_pputil.f90

gem_com.o: gem_com.f90 gem_pputil.o
	$(F90) -c $(OPT) $(INCLUDES) gem_com.f90

gem_equil.o: gem_equil.f90 gem_pputil.o
	$(F90) -c $(OPT) $(INCLUDES) gem_equil.f90

gem_gkps_adi.o: gem_gkps_adi.f90 gem_com.f90 gem_equil.f90 gem_pputil.f90
	$(F90) -c $(OPT) $(INCLUDES) gem_gkps_adi.f90

gem_main.o: gem_main.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o
	$(F90) -c $(OPT) $(INCLUDES) gem_main.f90

gem_outd.o: gem_outd.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o
	$(F90) -c $(OPT) $(INCLUDES) gem_outd.f90

gem_fcnt.o: gem_fcnt.f90
	$(F90) -c $(OPT) $(INCLUDES) gem_fcnt.f90

gem_fft_wrapper.o: gem_fft_wrapper.f90
	$(F90) -c $(OPT) $(INCLUDES) gem_fft_wrapper.f90

clean:
	rm -f *.o *.lst *.mod gem_main
