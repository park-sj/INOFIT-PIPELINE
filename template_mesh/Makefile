OBJ=main.o mmio.o gpu_solve.o
TEST=test

CFLAGS=-g -std=c++14 -I/home/shkim/Libraries/template/CGAL-5.0.3/include -I/home/shkim/Libraries/template/eigen-3.3.7/ -I/home/shkim/Libraries/template/eigen-3.3.7/unsupported/ -I/home/shkim/Libraries/template/boost_1_73_0/ -I/home/shkim/Libraries/cuda_cg/include

CXXFLAGS=$(CFLAGS)
NVCCFLAGS=-arch sm_20
LDFLAGS=-lgmp -lmpfr -lcuda -lcudart -lCGAL -L/usr/local/cuda/lib64 -L/home/shkim/Libraries/cuda_cg/lib -lyaml-cpp -g 
CXX=g++
CC=gcc
NVCC=nvcc

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) -c $^

%.o:	%.c
	$(CC) $(CFLAGS) -c $^

%.o:	%.cu
	$(NVCC) $(NVCCFLAGS) -c $^

all: $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $(TEST)

clean:
	rm -f *.o 
