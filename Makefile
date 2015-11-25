######################################################
#
#	Makefile
#
#	libr3d.a
#
#	See Readme.md for usage
#
######################################################

######## User options #########

# Use single-precision computations
#OPT += -DSINGLE_PRECISION

###############################

CC = gcc
CFLAGS = -Wall -I. -O3 
SRC = r3d.c r2d.c rNd.c v3d.c
DEPS = r3d.h r2d.h rNd.h v3d.h Makefile
OBJ = $(SRC:.c=.o)

all: libr3d.a

libr3d.a: $(OBJ)
	ar -rs $@ $^

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OPT)

clean:
	rm -rf libr3d.a $(OBJ) 
