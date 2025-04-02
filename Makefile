CC = g++ -Wall -ggdb
CC = g++ -pg
CC = mpiicpx
# CC = gcc

#mpirun -np 8 ./cdhit -i -o -T -c

# default with OpenMP
# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),no)
  CCFLAGS = -DNO_OPENMP
else
  CCFLAGS = -fopenmp
endif

#LDFLAGS = -static -lz -o
#LDFLAGS = /usr/lib/x86_64-linux-gnu/libz.a -o

# default with zlib
# without zlib
# in command line:
# make zlib=no
ifeq ($(zlib),no)
  CCFLAGS += 
  LDFLAGS += -o
else
  CCFLAGS += -DWITH_ZLIB
  LDFLAGS += -lz -o
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb
else
CCFLAGS += -O2
endif

ifdef MAX_SEQ
CCFLAGS += -DMAX_SEQ=$(MAX_SEQ)
endif

# PROGS = cd-hit cd-hit-est cd-hit-2d cd-hit-est-2d cd-hit-div cd-hit-454 cd-hit-lib-test
PROGS = cd-hit

# Propagate hardening flags
CCFLAGS := $(CPPFLAGS) $(CCFLAGS) $(CXXFLAGS)

# Object files
OBJS = cdhit-common.o cdhit-utility.o cdhit-lib.o

.c++.o:
	$(CC) $(CCFLAGS) -c $<

# all: $(PROGS) $(LIBTARGET)

# clean:
# 	rm -f *.o $(PROGS) $(LIBTARGET)

all: $(PROGS)

clean:
	rm -f *.o $(PROGS)


# programs

cd-hit: cdhit-common.o cdhit-utility.o cdhit.o
	$(CC) $(CCFLAGS) cdhit.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h input_sequence.h
	$(CC) $(CCFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-utility.c++ -c

cdhit.o: cdhit.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit.c++ -c

PREFIX ?= /usr/local/bin

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
	install -m 0755 *.pl $(PREFIX);
