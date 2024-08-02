# What compiler and linker to use?
CC=mpicc
LD=$(CC)

# Where to install things
TPLIBDIR=/biosoft/lib/tpatterns
BINDIR=/biosoft/bin
MANDIR=/biosoft/man/man1

# If you already have PCRE 1.09 or later compiled elsewhere, you may
# wish to change the following appropriately
PCRE_DIR=./pcre-1.09
# Location of PCRE header file
PCRE_INC=$(PCRE_DIR)
# Location of compiled PCRE library
PCRE_LIBDIR=$(PCRE_DIR)
# Location of compiled PCRE libraries for dependency checking
# (you can comment this out if you have a precompiled PCRE elsewhere)
PCRE_LIBS=$(PCRE_LIBDIR)/libpcre.a

# Add any extra compilation flags
CFLAGS=-O2
LFLAGS=-O2

# Locations of some common programs
CP=/bin/cp
MKDIR=/bin/mkdir
RM=/bin/rm
EGREP=/usr/bin/egrep
LN=/sbin/ln -s

#
#### YOU SHOULD NOT NEED TO CHANGE ANYTHING BELOW HERE ####
#

PCRE_LFLAGS=$(shell pcre-config --libs)
PCRE_CFLAGS=$(shell pcre-config --cflags)

# Some makes (notably IRIX) fail to get the shell right
SHELL=/bin/sh

# Version of tpatterns?
VERSION=2.0

OBJS=compare.o fastlibs.o main.o options.o tim.o

TPOBJS=translate.o errors.o fasta.o tpglobals.o
TPLIB=libtp.a

all:	txfasta tpatterns fasta2t

# This target is for my use only
dist:	clean
	( cd .. ; gtar zcvvf tpatterns-$(VERSION).tar.gz tpatterns-$(VERSION) )
	mv ../tpatterns-$(VERSION).tar.gz ../../public_html/software/tpatterns

clean:
	$(RM) -f *.[oa] txfasta tpatterns tpatterns.1 *~

$(TPLIB):	$(TPOBJS)
	ar -cr $@ $(TPOBJS)

tpatterns:	$(TPLIB) $(OBJS)
	$(LD) $(OBJS) $(TPLIB) $(PCRE_LFLAGS) $(LFLAGS) -o $@

txfasta:	$(TPLIB) txfasta.o
	$(LD) txfasta.o $(TPLIB) $(LFLAGS) -o $@

fasta2t:	$(TPLIB) fasta2t.o
	$(LD) fasta2t.o $(TPLIB) $(LFLAGS) -o $@

install:	txfasta tpatterns tpatterns.1
	( [ -d $(TPLIBDIR) ] || $(MKDIR) -p $(TPLIBDIR) )
	$(CP) universal.txm gcg2txm.pl $(TPLIBDIR)
	$(CP) txfasta $(BINDIR)
	$(CP) tpatterns $(BINDIR)
	$(CP) tpatterns.1 $(MANDIR)
	$(LN) $(MANDIR)/tpatterns.1 $(MANDIR)/txfasta.1

tpatterns.1:	tpatterns.preman
	$(CC) -E -DTF_LIB=$(TPLIBDIR) -DTF_BIN=$(BINDIR) \
	tpatterns.preman | $(EGREP) -v '^#' > $@

# General rule for all other C source files
.c.o:	tpatterns.h tplib.h
	$(CC) -c $(CFLAGS) $(PCRE_CFLAGS) \
	-DTFVER=\"$(VERSION)\" -DTFLIB=\"$(TPLIBDIR)\" $<
