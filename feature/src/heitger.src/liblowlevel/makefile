.SUFFIXES:
.SUFFIXES: .c .o

MAKEFILE = Makefile

INCPATH  = ..
INCFLAGS = -Iinclude -I$(INCPATH)


CFLAGS =  -DDEBUG_ALL $(G)
CC = cc

LIBPATH = ..

INCLUDE	= $(INCPATH)/lowlevel.h 
		
OBJ =\
	egabor.o \
	kpt.o \
	maps.o \
	ll_processing.o \
	noise_estimation.o

compile.all:   $(OBJ) $(INCLUDE) 
	ar r $(LIBPATH)/liblowlevel.a *.o
			
%.o: %.c $(INCLUDE)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(@F:.o=.c)
