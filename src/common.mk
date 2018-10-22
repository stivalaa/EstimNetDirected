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

CDEBUG = -g -DDEBUG_CONFIG  -DDEBUG_SAMPLER  -DDEBUG_DIGRAPH -DDEBUG_ALGS -DDEBUG_SNOWBALL -DDEBUG_MEMUSAGE
# Do NOT use -ffast-math as we depend on IEEE handling of NaN
OPTFLAGS = -O3  #-pg
CFLAGS     = $(OPTFLAGS) $(WARNFLAGS)

CPPFLAGS =  -DUSE_RANDOM123 -IRandom123-1.09/include
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

