
# Copyright (c) 2012, Jyri J. Virkki
# All rights reserved.
#
# This file is under BSD license. See LICENSE file.
#
# By default, builds optimized 64bit libbloom (under ./build)
# Requires GNU Make, so invoke appropriately (make or gmake)
#
# Other build options:
#
#   DEBUG=1 make        to build debug instead of optimized
#
# Other build targets:
#
#   make test           to build and run test code
#   make gcov           to build with code coverage and run gcov
#   make lint           to run lint
#   make clean          the usual
#

TOP := $(shell /bin/pwd)
BUILD_OS := $(shell uname)

BUILD=$(TOP)/build
INC=-I$(TOP) -I$(TOP)/murmur2
CC=g++ -Wall ${OPT}


ifeq ($(DEBUG),1)
OPT=-g $(DEBUGOPT)
else
OPT=-O3
endif


all: $(BUILD)/libbloom.a 

$(BUILD)/libbloom.a: $(BUILD)/murmurhash2.o $(BUILD)/bloom64.o
	(cd $(BUILD) && ar rcs libbloom.a bloom64.o murmurhash2.o)
	ranlib $(BUILD)/libbloom.a

$(BUILD)/test-libbloom: $(BUILD)/libbloom.a $(BUILD)/test.o
	(cd $(BUILD) && $(CC) test.o -L$(BUILD) -lbloom -o test-libbloom)

$(BUILD)/%.o: %.cpp
	mkdir -p $(BUILD)
	$(CC) $(INC) -c $< -o $@

$(BUILD)/murmurhash2.o: murmur2/MurmurHash2.c murmur2/murmurhash2.h
	mkdir -p $(BUILD)
	$(CC) $(INC) -c murmur2/MurmurHash2.c -o $(BUILD)/murmurhash2.o

clean:
	rm -rf $(BUILD)

lint:
	lint -x -errfmt=simple $(INC) $(LIB) *.c murmur2/*.c

test: $(BUILD)/test-libbloom
	$(BUILD)/test-libbloom

gcov:
	$(MAKE) clean
	DEBUG=1 DEBUGOPT="-fprofile-arcs -ftest-coverage" $(MAKE) all
	(cd $(BUILD) && \
		cp ../*.c . && \
		./test-libbloom && \
		gcov -bf bloom64.cpp)
	@echo Remember to make clean to remove instrumented objects
