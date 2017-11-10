CXX=mpicxx
CXXFLAGS=-std=c++14 -O3 -Wall

all: a1

demo:
	$(CXX) -DDEBUG $(CXXFLAGS) aAlt.cpp -o a1
	mpiexec -np 8 ./a1 8 10 out

clean:
	rm -rf a1
