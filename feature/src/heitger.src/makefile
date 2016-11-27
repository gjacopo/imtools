.SUFFIXES:
.SUFFIXES: .c .o

MAKEFILE = Makefile

CFLOWLIBS=  -L. -llowlevel -ltools  -lm
CFLAGS   = -Bstatic

LDFLAGS  =  -lm
G        =  
CC       = cc
 
all: 	 SE Keypoints Edgethresh 
	
SE:     SE.o 
	$(CC) $(G) -o $@ $@.o $(CFLAGS)  $(CFLOWLIBS)
		
Keypoints: Keypoints.o 
	$(CC) $(G) -o $@ $@.o $(CFLAGS) $(CFLOWLIBS)

Edgethresh: Edgethresh.o 
	$(CC) $(G) -o $@ $@.o $(CFLAGS) $(CFLOWLIBS)

%.o:    %.c 
	$(CC) $(G) -c $(CFLAGS) $(INCFLAGS) $(@F:.o=.c)
