LDFLAGS	=	-lm # -g
CPPFLAGS	=	-O3 -Wall -Wextra
CC      =	gcc
CPP	=	g++
RM      =	rm
INSTALL =	/usr/bin/install
INSTDIR =	/usr/local/bin

RIPLRMLEVOBJS	= riplrmlev.o

PROG	= riplrmlev

all:		$(PROG)

riplrmlev:	$(RIPLRMLEVOBJS)
		$(CPP) $(RIPLRMLEVOBJS) $(LDFLAGS) -o riplrmlev

clean:
		$(RM) -f $(RIPLRMLEVOBJS) $(PROG)

install:
		$(INSTALL) -m 775 $(PROG) $(INSTDIR)/


riplrmlev.o: riplrmlev.cpp riplrmlev.h
