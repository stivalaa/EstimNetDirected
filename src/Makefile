##############################################################################
#
# File      : Makefile
# Author    : Alex Stivala
# Created   : February 2009
#
#
# The default target "all" will compile and link the program.
#
# The "clean" target will remove object and executable files.
#
# set MODE=DEBUG to build with debugging and verbose printing on host,
# default is to build with optimizations on and no debug or profile
#
#
#
##############################################################################

include common.mk
-include local.mk


DEPENDFILE = .depend

COMMON_C_OBJS =  utils.o digraph.o \
                 changeStatisticsDirected.o basicSampler.o \
                 equilibriumExpectation.o configparser.o estimconfigparser.o \
                 ifdSampler.o

MPI_C_OBJS    = EstimNetDirectedMPImain.o
NONMPI_C_OBJS = EstimNetDirectedMain.o

COMMON_C_SRCS = $(COMMON_C_OBJS:.o=.c)
MPI_C_SRCS    = $(MPI_C_OBJS:.o=.c)
NONMPI_C_SRCS = $(NONMPI_C_OBJS:.o=.c)

COMMON_HASH_C_OBJS = $(COMMON_C_OBJS:.o=_hash.o)

OBJS = $(COMMON_C_OBJS) $(MPI_C_OBJS) $(NONMPI_C_OBJS) $(COMMON_HASH_C_OBJS)
SRCS = $(COMMON_C_SRCS) $(MPI_C_SRCS) $(NONMPI_C_SRCS)

all: EstimNetDirected EstimNetDirected_mpi EstimNetDirected_hashtables EstimNetDirected_mpi_hashtables

EstimNetDirected: $(COMMON_C_OBJS) $(NONMPI_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

EstimNetDirected_mpi: $(COMMON_C_OBJS) $(MPI_C_OBJS)
	$(MPILD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

# These versions use hash tables instead of arrays for two-paths lookup
EstimNetDirected_hashtables: $(COMMON_HASH_C_OBJS) $(NONMPI_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

EstimNetDirected_mpi_hashtables: $(COMMON_HASH_C_OBJS) $(MPI_C_OBJS)
	$(MPILD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

tags: $(SRCS)
	$(MAKETAGS) $(SRCS)

depend: $(SRCS)
	$(MAKEDEPEND) $(SRCS) > $(DEPENDFILE)

clean:
	rm -f $(OBJS)
	rm -f EstimNetDirected EstimNetDirected_mpi
	rm -f EstimNetDirected_hashtables EstimNetDirected_mpi_hashtables


EstimNetDirectedMPImain.o: EstimNetDirectedMPImain.c
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_HASHTABLES -c -o $@ $<



###############################################################################

-include $(DEPENDFILE)


