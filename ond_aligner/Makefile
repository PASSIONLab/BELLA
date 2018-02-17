DEST_DIR = ./bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL= ondaligner.o

all: $(ALL)

ondaligner.o: align.c align.h
	gcc $(CFLAGS) -lm -o $@ $<

clean:
	rm $(ALL)