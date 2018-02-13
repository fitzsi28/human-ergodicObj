CXX=g++
CXXFLAGS= -std=c++11 -Wall -g

main:main.o rk4_int.o
	$(CXX) $(CXXFLAGS) -o main main.o rk4_int.o
main.o: rk4_int.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp rk4_int.cpp
rk4_int.o: rk4_int.hpp
	$(CXX) $(CXXFLAGS) -c rk4_int.cpp
.PHONY: clean
# Erase all hex, map, object, and elf files.
clean :
	$(RM) *.hex *.map *.o *.elf *.dep *.dis