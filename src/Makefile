# A standard compiler on most unix systems.

CC = gcc

# your general purpose optimization flag

CFLAGS = -O2

# Use this next section as your CFLAGS for more expensive
# optimizations on the 486 and above family of processors.

# CFLAGS = -s -static -O3 -fomit-frame-pointer -m486  \
#	-malign-functions=4 -funroll-loops -fexpensive-optimizations \
#	-malign-double -fschedule-insns2 -mwide-multiply

all: ../exe/bent_umbrella ../exe/correlate ../exe/octopus ../exe/visualize

../exe/visualize: visualize.o
	$(CC) $(CFLAGS) -o ../exe/visualize visualize.o -lm

../exe/octopus: octopus.o
	$(CC) $(CFLAGS) -o ../exe/octopus octopus.o -lm

../exe/correlate: correlate.o
	$(CC) $(CFLAGS) -o ../exe/correlate correlate.o -lm

../exe/bent_umbrella: bent_umbrella.o circle_intersect.o
	$(CC) $(CFLAGS) -o ../exe/bent_umbrella bent_umbrella.o circle_intersect.o -lm

visualize.o: parameters.h visualize.c

octopus.o: circle_intersect.h parameters.h octopus.c

correlate.o: circle_intersect.h parameters.h correlate.c

bent_umbrella.o: parameters.h circle_intersect.h bent_umbrella.c

circle_intersect.o: parameters.h circle_intersect.h circle_intersect.c

clean:
	rm *.o
