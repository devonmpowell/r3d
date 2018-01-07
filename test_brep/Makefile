######################################################
#
#	Makefile
#
#	for r3d robustness unit testing
#
#	usage: `make`
#
######################################################

# where r3d lives
R3D_HOME = ..

# what we want to build
R3D_SRC = test_intersection.cc

# the compiler
CC = g++
#CC = icpc

# compiler flages
CFLAGS = -g -Wall -I. -std=c++11

# additional include directories
INC = -I$(R3D_HOME)

# additional link directories
LIB = -L$(R3D_HOME)

# additional link flags
LDFLAGS = -lr3d -lm

# define the object files
R3D_OBJ = $(R3D_SRC:.cc=.o)

all: test_intersection

test_intersection: $(R3D_OBJ)
	$(CC) $(LIB) $(CFLAGS) -o $@ $(R3D_OBJ) $(LDFLAGS) 

%.o: %.cc
	$(CC) $(INC) -c $(CFLAGS) -o $@ $<  

clean:
	rm -rf test_intersection  $(R3D_OBJ) 
