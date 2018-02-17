DEST_DIR = ./bin

CC=gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL= ondaligner

all: $(ALL)
	mkdir $(DEST_DIR)
	mv $(ALL) $(DEST_DIR)

ondaligner: driver.c align.c align.h
	$(CC) $(CFLAGS) -lm -o $@ driver.c align.c

clean:
	rm -r $(DEST_DIR)