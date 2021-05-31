###############################################################################
#
# File      : common.mk
# Author    : Alex Stivala
# Created   : July 2008
#
#
# Definitions of compilers and compile options for all Makefiles.
# Used for compiling with MPI and POSIX threads.
# Use GNU make.
#
# set MODE=DEBUG to build with debugging and verbose printing on host,
# default is to build with optimizations on and no debug or profile
#
###############################################################################

CC           = gcc
LD           = gcc

# must ensure that MPICC uses the same compiler CC 
# e.g. for gcc must "module load openmpi-gcc" to get OpenMPI using gcc
# on tango.vpac.org
MPICC        = mpicc
MPILD        = mpicc

WARNFLAGS = -Wall
#              the following warnings are not implied by -Wall
WARNFLAGS  += -Wextra -Wfloat-equal  \
              -Wundef -Wshadow \
              -Wpointer-arith -Wcast-qual -Wcast-align\
              -Wwrite-strings \
              -Wmissing-declarations -Wunreachable-code

CDEBUG = -g -DDEBUG_CONFIG  -DDEBUG_SAMPLER  -DDEBUG_DIGRAPH -DDEBUG_ALGS -DDEBUG_SNOWBALL -DDEBUG_MEMUSAGE -DDEBUG_SIMULATE -DDEBUG_CERGM
# Do NOT use -ffast-math as we depend on IEEE handling of NaN
OPTFLAGS = -O3  #-pg
CFLAGS     = $(OPTFLAGS) $(WARNFLAGS)

# Use the Random123 library
CPPFLAGS =  -DUSE_RANDOM123 -IRandom123-1.09/include

## Use lookup table for pow()
#CPPFLAGS += -DUSE_POW_LOOKUP

# Using the uthash hash table 
# See https://troydhanson.github.io/uthash/userguide.html
# and https://github.com/troydhanson/uthash
# Enable the Bloom filter (max size 32 bits = 512 MB)
# for faster misses
# as this was found to increase speed significantly.
# Using default hash function (Jenkins) as get different
# rankings for different networks using the hashscan
# and keytests utilities so no basis to choose a generally
# best hash fucntion it seems
CPPFLAGS += -DHASH_BLOOM=32

# extra output even on optimized build for now
#CPPFLAGS += -DDEBUG_MEMUSAGE -DDEBUG_SNOWBALL

LDFLAGS =   #-pg


ifeq ($(MODE),DEBUG)
  CFLAGS = $(CDEBUG) $(WARNFLAGS)
endif

PTHREAD_CFLAGS = -pthread
PTHREAD_LDFLAGS = -pthread


LDLIBPATH  = 
LDLIBS     = -lm



# Program to build TAGS file for EMACS
MAKETAGS   = etags

MAKEDEPEND = gcc -MM -x c++ $(CPPFLAGS)

