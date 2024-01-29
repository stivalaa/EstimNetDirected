include ../../src/common.mk
-include ../../src/local.mk

SRCDIR = ../../src

CC = gcc

CFLAGS += -I$(SRCDIR)

OBJS = $(SRCDIR)/changeStatisticsGeneral.o  $(SRCDIR)/changeStatisticsUndirected.o $(SRCDIR)/graph.o $(SRCDIR)/utils.o $(SRCDIR)/loadGraph.o $(SRCDIR)/changeStatisticsBipartiteUndirected.o
HASH_OBJS = $(OBJS:.o=_hash.o)
ARRAY_OBJS = $(OBJS:.o=_array.o)

all: testChangeStatsUndirected testChangeStatsUndirected_hash testChangeStatsUndirected_array testStatsSumChangeStatsUndirected testStatsSumChangeStatsUndirected_hash testStatsSumChangeStatsUndirected_array


testChangeStatsUndirected: testChangeStatsUndirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsUndirected_hash: testChangeStatsUndirectedMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsUndirected_array: testChangeStatsUndirectedMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm


testStatsSumChangeStatsUndirected: testStatsSumChangeStatsUndirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testStatsSumChangeStatsUndirected_hash: testStatsSumChangeStatsUndirectedMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testStatsSumChangeStatsUndirected_array: testStatsSumChangeStatsUndirectedMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm


clean:
	$(RM) testChangeStatsUndirected testChangeStatsUndirected_hash
	$(RM) testChangeStatsUndirectedMain.o testChangeStatsUndirectedMain_hash.o
	$(RM) testChangeStatsUndirected_array testChangeStatsUndirectedMain_array.o
	$(RM) testStatsSumChangeStatsUndirected testStatsSumChangeStatsUndirected_hash
	$(RM) testStatsSumChangeStatsUndirectedMain.o testStatsSumChangeStatsUndirectedMain_hash.o
	$(RM) testStatsSumChangeStatsUndirected_array testStatsSumChangeStatsUndirectedMain_array.o


%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<
