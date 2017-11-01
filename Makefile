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

CC = g++
CFLAGS = -Wall -I. -O3 -shared -fPIC 
SRC = r3d.c r2d.c rNd.c v3d.c v2d.c vNd.c
DEPS = r3d.h r2d.h rNd.h v3d.h v2d.h vNd.h Makefile
OBJ = $(SRC:.c=.o)

all: libr3d.a tests

libr3d.a: $(OBJ)
	ar -rs $@ $^

tests:
	make -C tests

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OPT)

clean:
	rm -rf libr3d.a $(OBJ); make clean -C tests 
