.SUFFIXES:
.SUFFIXES: .c .o

MAKEFILE = Makefile

INCPATH  = ..
INCFLAGS = -Iinclude -I$(INCPATH)


CFLAGS =  -DDEBUG_ALL $(G)
CC = cc

LIBPATH = ..

INCLUDE	= $(INCPATH)/tools.h 
		
OBJ =\
	io.o \
	tools.o \
	pgon_approx.o

compile.all:   $(OBJ) $(INCLUDE) 
	ar r $(LIBPATH)/libtools.a *.o
			
%.o: %.c $(INCLUDE)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(@F:.o=.c)
