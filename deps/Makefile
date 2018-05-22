CC = mpicc
CXX = mpicxx
F77 = mpif77
FC = mpif90
CFLAGS = --std=gnu11 -g -O2

CFLAGS +=-fno-common -fomit-frame-pointer

CFLAGS += -Wconversion -Wno-sign-conversion \
          -Wcast-align -Wchar-subscripts -Wall -W \
          -Wpointer-arith -Wwrite-strings -Wformat-security -pedantic \
          -Wextra -Wno-unused-parameter

CFLAGS +=-DUSE_DFLOAT_DOUBLE

# list of libraries to build
UNAME_S := $(shell uname -s)

# p4est flags
CPPFLAGS += -Ip4est/local/include
LDFLAGS += -Lp4est/local/lib
LDLIBS += -lp4est -lsc
ifeq ($(UNAME_S),Linux)
 LDFLAGS += -Wl,-rpath=$(CURDIR)/p4est/local/lib,--enable-new-dtags
endif

all:
	tar xzf downloads/p4est-*.tar.gz && mv p4est-* p4est
	cd p4est && ./configure CC=$(CC) CXX=$(CXX) F77=$(F77) FC=$(FC) --enable-mpi --enable-debug --without-blas \
           && $(MAKE) install

realclean:
	rm -rf p4est
