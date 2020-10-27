CFLAGS = -Wall
LIBS = -lm -lSDL2 -lconfig
CC = gcc

all: visualize simulate

visualize: visual.c
	$(CC) $(CFLAGS) $(LIBS) $? -o $@

simulate: trial.c
	$(CC) $(CFLAGS) $(LIBS) $? -o $@

