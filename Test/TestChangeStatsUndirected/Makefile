include ../../src/common.mk
-include ../../src/local.mk

SRCDIR = ../../src

CC = gcc

CFLAGS += -I$(SRCDIR)

OBJS = $(SRCDIR)/changeStatisticsGeneral.o  $(SRCDIR)/changeStatisticsUndirected.o $(SRCDIR)/graph.o $(SRCDIR)/utils.o $(SRCDIR)/loadGraph.o
HASH_OBJS = $(OBJS:.o=_hash.o)
ARRAY_OBJS = $(OBJS:.o=_array.o)

all: testChangeStatsUndirected testChangeStatsUndirected_hash testChangeStatsUndirected_array


testChangeStatsUndirected: testChangeStatsUndirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsUndirected_hash: testChangeStatsUndirectedMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsUndirected_array: testChangeStatsUndirectedMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm



clean:
	$(RM) testChangeStatsUndirected testChangeStatsUndirected_hash
	$(RM) testChangeStatsUndirectedMain.o testChangeStatsUndirectedMain_hash.o
	$(RM) testChangeStatsUndirected_array testChangeStatsUndirectedMain_array.o

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<
