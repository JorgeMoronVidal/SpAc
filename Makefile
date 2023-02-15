# CC : C compiler

CC	= mpigxx

FLAGS = -g -march=native -O -DNDEBUG -O3

MUMPS = -I/home/jorge/Dropbox/DOC/MUMPS_TEST/MUMPS_5.5.1/include -Xlinker -rpath=/home/jorge/Dropbox/DOC/MUMPS_TEST/MUMPS_5.5.1/lib -L/home/jorge/Dropbox/DOC/MUMPS_TEST/MUMPS_5.5.1/lib -ldmumps -lmumps_common -lmetis \
    	-Xlinker -rpath=/home/jorge/Dropbox/DOC/MUMPS_TEST/MUMPS_5.5.1/lib -L/home/jorge/Dropbox/DOC/MUMPS_TEST/MUMPS_5.5.1/lib -lpord

MKL = -L/opt/intel/oneapi/mkl/2023.0.0/lib/intel64  -lmkl_intel_lp64 -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_core


OPTC = -fPIC

# $(wildcard *.cpp /xxx/xxx/*.cpp): get all .cpp files from the current directory and dir "/xxx/xxx/"
SRCS := $(wildcard */*.cpp)
# $(patsubst %.cpp,%.o,$(SRCS)): substitute all ".cpp" file name strings to ".o" file name strings
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

all: main class

class:$(OBJS)
	@echo $(OBJS)

%.o : %.cpp
	$(CC) $(FLAGS) -c $< -o $@

main_parallel:
	$(CC) $(OPTC) -o test_mpi $(FLAGS) main_test_mpi.cpp $(MUMPS) $(MKL) $(OBJS) -lgfortran -liomp5

main: class
	$(CC) $(OPTC)  main_test_mpi.cpp -o main $(FLAGS) $(MUMPS) $(OBJS)

main_serial:
	g++ $(OPTC) main_test_serial.cpp -o main $(FLAGS) $(MUMPS) $(OBJS)

clean_class:
	rm $(OBJS)

clean:
	rm main
