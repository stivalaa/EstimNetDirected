include ../../src/common.mk
-include ../../src/local.mk

SRCDIR = ../../src

CC = gcc

CFLAGS += -I$(SRCDIR)

LOCALOBJS = bipartiteStats.o
OBJS = $(SRCDIR)/changeStatisticsGeneral.o  $(SRCDIR)/changeStatisticsBipartiteUndirected.o $(SRCDIR)/graph.o $(SRCDIR)/utils.o $(SRCDIR)/loadGraph.o $(SRCDIR)/changeStatisticsUndirected.o $(LOCALOBJS)
HASH_OBJS = $(OBJS:.o=_hash.o)
ARRAY_OBJS = $(OBJS:.o=_array.o)

all: testChangeStatsBipartite testChangeStatsBipartite_hash testChangeStatsBipartite_array testStatsSumChangeStatsBipartite testStatsSumChangeStatsBipartite_hash testStatsSumChangeStatsBipartite_array testDiffStatsChangeStatsBipartite testDiffStatsChangeStatsBipartite_hash testDiffStatsChangeStatsBipartite_array


testChangeStatsBipartite: testChangeStatsBipartiteMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsBipartite_hash: testChangeStatsBipartiteMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsBipartite_array: testChangeStatsBipartiteMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm



testStatsSumChangeStatsBipartite: testStatsSumChangeStatsBipartiteMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testStatsSumChangeStatsBipartite_hash: testStatsSumChangeStatsBipartiteMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testStatsSumChangeStatsBipartite_array: testStatsSumChangeStatsBipartiteMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm



testDiffStatsChangeStatsBipartite: testDiffStatsChangeStatsBipartiteMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testDiffStatsChangeStatsBipartite_hash: testDiffStatsChangeStatsBipartiteMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testDiffStatsChangeStatsBipartite_array: testDiffStatsChangeStatsBipartiteMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm




clean:
	$(RM) testChangeStatsBipartite testChangeStatsBipartite_hash
	$(RM) testChangeStatsBipartiteMain.o testChangeStatsBipartiteMain_hash.o
	$(RM) testChangeStatsBipartite_array testChangeStatsBipartiteMain_array.o
	$(RM) testStatsSumChangeStatsBipartite testStatsSumChangeStatsBipartite_hash
	$(RM) testStatsSumChangeStatsBipartiteMain.o testStatsSumChangeStatsBipartiteMain_hash.o
	$(RM) testStatsSumChangeStatsBipartite_array testStatsSumChangeStatsBipartiteMain_array.o
	$(RM) testDiffStatsChangeStatsBipartite testDiffStatsChangeStatsBipartite_hash
	$(RM) testDiffStatsChangeStatsBipartiteMain.o testDiffStatsChangeStatsBipartiteMain_hash.o
	$(RM) testDiffStatsChangeStatsBipartite_array testDiffStatsChangeStatsBipartiteMain_array.o
	$(RM) $(LOCALOBJS) $(LOCALOBJS:.o=_hash.o) $(LOCALOBJS:.o=_array.o)

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<
