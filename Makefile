CC = g++ -Wall -ggdb
CC = g++ -pg
CC = mpiicpx

# default with OpenMP
# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),no)
  CCFLAGS = -DNO_OPENMP
else
  CCFLAGS = -fopenmp
endif
ifeq ($(AVX512),yes)
#   CCFLAGS += -mavx512f -mavx512vl -mavx512bw -mavx512dq -fno-vectorize 
CCFLAGS += -march=native
else
   CCFLAGS += -DNO_AVX512
endif
#LDFLAGS = -static -lz -o
#LDFLAGS = /usr/lib/x86_64-linux-gnu/libz.a -o

# default with zlib
# without zlib
# in command line:
# make zlib=no
ifeq ($(zlib),no)
  CCFLAGS += -g
  LDFLAGS += -o
else
  CCFLAGS += -DWITH_ZLIB -g
  LDFLAGS += -lz -o
endif
# CCFLAGS += -DREADY
# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb
else
CCFLAGS += -Ofast
endif

ifdef MAX_SEQ
CCFLAGS += -DMAX_SEQ=$(MAX_SEQ)
endif

PROGS = cdhit-mpi cdhit-preprocess

# Propagate hardening flags
CCFLAGS := $(CPPFLAGS) $(CCFLAGS) $(CXXFLAGS)

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)

clean:
	rm -f *.o $(PROGS)

# programs
cdhit-mpi: cdhit-common.o cdhit-utility.o cdhit-mpi.o
	$(CC) $(CCFLAGS) cdhit-mpi.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cdhit-mpi

cdhit-preprocess: cdhit-common.o cdhit-utility.o cdhit-preprocess.o
	$(CC) $(CCFLAGS) cdhit-preprocess.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cdhit-preprocess

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-utility.c++ -c

cdhit-mpi.o: cdhit-mpi.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-mpi.c++ -c

cdhit-preprocess.o: cdhit-preprocess.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-preprocess.c++ -c

PREFIX ?= /usr/local/bin

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
	install -m 0755 *.pl $(PREFIX);
