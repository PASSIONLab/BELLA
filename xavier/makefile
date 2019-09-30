CC = gcc
CXX = g++
CFLAGS  = 	-I. -O3 -Wall -Wextra -pedantic -ansi -c -Wno-write-strings

demo: demo.cpp
	$(CXX) -march=native -std=c++14 -O3 -fpermissive -o demo demo.cpp -DDEBUG

clean:
	rm -f *.o
	rm -f demo

