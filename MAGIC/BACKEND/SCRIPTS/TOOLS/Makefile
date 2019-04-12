#VPATH = ./gestures:./collision

CC=gcc

F90=ifort
F90=gfortran

FLAGS90=-C -cpp
FLAGS90=-O3 -cpp

# -std=c99 is needed for roundf and -D_GNU_SOURCE
#  to have getopt & strndup functions exported by
#  unistd.h and string.h (respectively) even with c99.
CFLAGS=-g -O3 -W -Wall -std=c99 -D_GNU_SOURCE -DDOUBLE_PREC
CFLAGS=-g -W -Wall -std=c99 -D_GNU_SOURCE -DDOUBLE_PREC

LFLAGS=
LFLAGS=-g

.SUFFIXES: .f90 .cu .mod .o

all : generate_mesh.o fillArr.so points2cells.so combine_lumen_plaques.so
	$(CC) ${LFLAGS} -o generate_mesh generate_mesh.o -lm

generate_mesh.o : generate_mesh.c
	$(CC) ${CFLAGS} -c -o generate_mesh.o generate_mesh.c

fillArr.so : fillArray.c
	$(CC) ${CFLAGS} -c fillArray.c -fPIC -g
	$(CC) ${CFLAGS} -shared -o fillArr.so fillArray.o

points2cells.so : points2cells.f90
	$(F90) -g ${FLAGS90} -c points2cells.f90 -fPIC -g
	$(F90) -g ${FLAGS90} -shared -o points2cells.so points2cells.o

combine_lumen_plaques.so : combine_lumen_plaques.f90
	$(F90) -g ${FLAGS90} -c combine_lumen_plaques.f90 -fPIC -g
	$(F90) -g ${FLAGS90} -shared -o combine_lumen_plaques.so combine_lumen_plaques.o
	$(F90) -g ${FLAGS90} -o tst combine_lumen_plaques.f90 tst.f90 

clean : 
	rm generate_mesh generate_mesh.o fillArr.so points2cells.so combine_lumen_plaques.so



%: %.o
	$(F90) $(FLAGS90) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(F90) $(FLAGS90) -c $<
