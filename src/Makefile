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

ESTIM_COMMON_C_OBJS =  utils.o digraph.o \
                 changeStatisticsDirected.o basicSampler.o \
                 equilibriumExpectation.o configparser.o estimconfigparser.o \
                 ifdSampler.o loadDigraph.o tntSampler.o

SIM_COMMON_C_OBJS =  utils.o digraph.o \
                 changeStatisticsDirected.o basicSampler.o \
                 configparser.o simconfigparser.o ifdSampler.o simulation.o \
                 tntSampler.o loadDigraph.o

ESTIM_MPI_C_OBJS    = EstimNetDirectedMPImain.o
ESTIM_NONMPI_C_OBJS = EstimNetDirectedMain.o

SIM_C_OBJS = SimulateERGMmain.o

ESTIM_COMMON_C_SRCS = $(ESTIM_COMMON_C_OBJS:.o=.c)
ESTIM_MPI_C_SRCS    = $(ESTIM_MPI_C_OBJS:.o=.c)
ESTIM_NONMPI_C_SRCS = $(ESTIM_NONMPI_C_OBJS:.o=.c)
SIM_C_SRCS          = $(SIM_C_OBJS:.o=.c)

ESTIM_COMMON_HASH_C_OBJS = $(ESTIM_COMMON_C_OBJS:.o=_hash.o)
SIM_COMMON_HASH_C_OBJS = $(SIM_COMMON_C_OBJS:.o=_hash.o)

ESTIM_COMMON_ARRAY_C_OBJS = $(ESTIM_COMMON_C_OBJS:.o=_array.o)
SIM_COMMON_ARRAY_C_OBJS = $(SIM_COMMON_C_OBJS:.o=_array.o)


OBJS = $(ESTIM_COMMON_C_OBJS) $(ESTIM_MPI_C_OBJS) $(ESTIM_NONMPI_C_OBJS) $(ESTIM_COMMON_HASH_C_OBJS) $(SIM_COMMON_C_OBJS) $(SIM_COMMON_HASH_C_OBJS) $(ESTIM_COMMON_ARRAY_C_OBJS) $(SIM_COMMON_ARRAY_C_OBJS) $(SIM_C_OBJS)
SRCS = $(ESTIM_COMMON_C_SRCS) $(ESTIM_MPI_C_SRCS) $(ESTIM_NONMPI_C_SRCS) $(SIM_C_SRCS)

all: EstimNetDirected EstimNetDirected_mpi EstimNetDirected_hashtables EstimNetDirected_mpi_hashtables SimulateERGM SimulateERGM_hashtables EstimNetDirected_mpi_arrays EstimNetDirected_arrays SimulateERGM_arrays

EstimNetDirected: $(ESTIM_COMMON_C_OBJS) $(ESTIM_NONMPI_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

EstimNetDirected_mpi: $(ESTIM_COMMON_C_OBJS) $(ESTIM_MPI_C_OBJS)
	$(MPILD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

SimulateERGM:  $(SIM_COMMON_C_OBJS) $(SIM_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)


# These versions use hash tables instead of arrays for two-paths lookup
EstimNetDirected_hashtables: $(ESTIM_COMMON_HASH_C_OBJS) $(ESTIM_NONMPI_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

EstimNetDirected_mpi_hashtables: $(ESTIM_COMMON_HASH_C_OBJS) $(ESTIM_MPI_C_OBJS)
	$(MPILD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

SimulateERGM_hashtables: $(SIM_COMMON_HASH_C_OBJS) $(SIM_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

# These versions use arrays for two-paths lookup
EstimNetDirected_arrays: $(ESTIM_COMMON_ARRAY_C_OBJS) $(ESTIM_NONMPI_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

EstimNetDirected_mpi_arrays: $(ESTIM_COMMON_ARRAY_C_OBJS) $(ESTIM_MPI_C_OBJS)
	$(MPILD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

SimulateERGM_arrays: $(SIM_COMMON_ARRAY_C_OBJS) $(SIM_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)





tags: $(SRCS)
	$(MAKETAGS) $(SRCS)

depend: $(SRCS)
	$(MAKEDEPEND) $(SRCS) > $(DEPENDFILE)

clean:
	rm -f $(OBJS)
	rm -f EstimNetDirected EstimNetDirected_mpi
	rm -f EstimNetDirected_hashtables EstimNetDirected_mpi_hashtables
	rm -f EstimNetDirected_arrays EstimNetDirected_mpi_arrays
	rm -f SimulateERGM SimulateERGM_arrays SimulateERGM_hashtables


EstimNetDirectedMPImain.o: EstimNetDirectedMPImain.c
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<



###############################################################################

-include $(DEPENDFILE)


