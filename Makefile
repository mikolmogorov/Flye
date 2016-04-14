SRC := $(shell pwd)/assemble
POLISH := $(shell pwd)/polishing

export LIBBF = $(shell pwd)/libbf
export BIN_DIR = $(shell pwd)
export CXXFLAGS = -std=c++0x -I${LIBBF}
export LDFLAGS = -L${LIBBF} -lbf

.PHONY: clean all

all: 
	make all -C ${LIBBF}
	make log -C ${SRC}
	make all -C ${POLISH}
debug:
	make debug -C ${LIBBF}
	make debug -C ${SRC}
	make debug -C ${POLISH}
log:
	make all -C ${LIBBF}
	make all -C ${POLISH}
	make log -C ${SRC}
clean: 
	make clean -C ${LIBBF}
	make clean -C ${POLISH}
	make clean -C ${SRC}