EXEDIR= ../../SHOP/
#VPATH = ./gestures:./collision
LFLAGS=
# -std=c99 is needed for roundf and -D_GNU_SOURCE
#  to have getopt & strndup functions exported by
#  unistd.h and string.h (respectively) even with c99.
CC= g++
CFLAGS= -W -Wall -std=c99 -D_GNU_SOURCE
CPPFLAGS= -W -Wall

ifdef debug
CPPFLAGS +=  -g -O0
CFLAGS +=  -g -O0
else
CPPFLAGS +=  -O
CFLAGS +=  -O3
endif

COMP =  $(CC) $(CPPFLAGS) -c


OBJECTS = genmesh.o readargs.o double3.o point.o node.o mesh.o 

all : generate_mesh.o $(OBJECTS)
	$(CC) $(CPPFLAGS) $(OBJECTS) $(LFLAGS) -o $(EXEDIR)/genmesh

nocomp:
	$(MAKE) genmesh "CPPFLAGS = $(CPPFLAGS) -DNOCOMP " 


generate_mesh.o : generate_mesh.c
	gcc ${CFLAGS} -c generate_mesh.c; gcc ${LFLAGS} -o $(EXEDIR)/generate_mesh generate_mesh.o -lm


genmesh.o:	genmesh.cpp genmesh.h point.h node.h mesh.h readargs.h
	   $(COMP) genmesh.cpp

readargs.o:	readargs.cpp readargs.h
	   $(COMP) readargs.cpp

double3.o:	double3.cpp double3.h
	   $(COMP) double3.cpp

point.o:	point.cpp point.h double3.h genmesh.h
	   $(COMP) point.cpp

node.o:		node.cpp node.h genmesh.h
	   $(COMP) node.cpp

mesh.o:		mesh.cpp mesh.h node.h point.h genmesh.h
	   $(COMP) mesh.cpp

genmesh: $(OBJECTS)
	$(CC) $(CPPFLAGS) $(OBJECTS) $(LFLAGS) -o $(EXEDIR)/genmesh

clean : 
	rm $(OBJECTS) generate_mesh.o
