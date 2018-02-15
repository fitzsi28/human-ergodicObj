CXX=g++
CXXFLAGS= -std=c++11 -Wall -g

main:main.o 
	$(CXX) $(CXXFLAGS) -o main main.o 
main.o: rk4_int.hpp cartpend.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp 

.PHONY: clean
# Erase all hex, map, object, and elf files.
clean :
	$(RM) *.hex *.map *.o *.elf *.dep *.dis