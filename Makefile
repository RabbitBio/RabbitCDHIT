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
  CCFLAGS = -fopenmp -D_REENTRANT -mcx16 -std=c++17
endif
ifeq ($(debug),yes)
  CCFLAGS += -DDEBUG
endif
ifeq ($(AVX512),yes)
#   CCFLAGS += -mavx512f -mavx512vl -mavx512bw -mavx512dq -fno-vectorize 
CCFLAGS += -march=native
else
   CCFLAGS += -DNO_AVX512 
endif

#LDFLAGS = -static -lz -o
#LDFLAGS = /usr/lib/x86_64-linux-gnu/libz.a -o


# TBB 配置 - 可以通过环境变量 TBB_ROOT 指定，或使用默认路径
TBB_ROOT ?= /opt/intel/oneapi/tbb/latest
TBB_LIB_DIR = $(TBB_ROOT)/lib/intel64/gcc4.8
TBB_INCLUDE_DIR = $(TBB_ROOT)/include

# 检查 TBB 是否存在，如果不存在则使用系统默认
ifeq ($(wildcard $(TBB_LIB_DIR)/libtbb.so),)
  # TBB 不在默认路径，尝试其他可能的位置
  TBB_LIB_DIR = $(shell if [ -d "$(TBB_ROOT)/lib/intel64" ]; then echo "$(TBB_ROOT)/lib/intel64"; \
                            elif [ -d "$(TBB_ROOT)/lib" ]; then echo "$(TBB_ROOT)/lib"; \
                            else echo ""; fi)
endif

# 添加 TBB 路径
ifneq ($(TBB_LIB_DIR),)
  LDFLAGS += -L$(TBB_LIB_DIR)
endif
ifneq ($(TBB_INCLUDE_DIR),)
  CCFLAGS += -I$(TBB_INCLUDE_DIR)
endif

# Intel oneAPI TBB 需要链接 tbb 和 tbbmalloc
LIBS = -ltbb -ltbbmalloc -latomic


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

SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

# 创建build目录
$(shell mkdir -p $(BUILD_DIR))

PROGS = cdhit-mpi cdhit-preprocess

# Propagate hardening flags
CCFLAGS := $(CPPFLAGS) $(CCFLAGS) $(CXXFLAGS)

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)

clean:
	rm -f $(BUILD_DIR)/*.o $(PROGS)

# programs
cdhit-mpi: $(BUILD_DIR)/cdhit-common.o $(BUILD_DIR)/cdhit-utility.o $(BUILD_DIR)/cdhit-mpi.o
	$(CC) $(CCFLAGS) $(BUILD_DIR)/cdhit-mpi.o $(BUILD_DIR)/cdhit-common.o $(BUILD_DIR)/cdhit-utility.o ${LIBS} $(LDFLAGS) cdhit-mpi

cdhit-preprocess: $(BUILD_DIR)/cdhit-common.o $(BUILD_DIR)/cdhit-utility.o $(BUILD_DIR)/cdhit-preprocess.o
	$(CC) $(CCFLAGS) $(BUILD_DIR)/cdhit-preprocess.o $(BUILD_DIR)/cdhit-common.o $(BUILD_DIR)/cdhit-utility.o ${LIBS} $(LDFLAGS) cdhit-preprocess

# objects
$(BUILD_DIR)/cdhit-common.o: $(SRC_DIR)/cdhit-common.c++ $(INC_DIR)/cdhit-common.h
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(SRC_DIR)/cdhit-common.c++ -c -o $@

$(BUILD_DIR)/cdhit-utility.o: $(SRC_DIR)/cdhit-utility.c++ $(INC_DIR)/cdhit-utility.h
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(SRC_DIR)/cdhit-utility.c++ -c -o $@

$(BUILD_DIR)/cdhit-mpi.o: $(SRC_DIR)/cdhit-mpi.c++ $(INC_DIR)/cdhit-utility.h
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(SRC_DIR)/cdhit-mpi.c++ -c -o $@

$(BUILD_DIR)/cdhit-preprocess.o: $(SRC_DIR)/cdhit-preprocess.c++ $(INC_DIR)/cdhit-utility.h
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(SRC_DIR)/cdhit-preprocess.c++ -c -o $@

PREFIX ?= /usr/local/bin

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
	install -m 0755 *.pl $(PREFIX);