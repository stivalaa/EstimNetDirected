include ../../src/common.mk
-include ../../src/local.mk

SRCDIR = ../../src

CC = gcc

CFLAGS += -I$(SRCDIR)

OBJS = $(SRCDIR)/changeStatisticsGeneral.o  $(SRCDIR)/changeStatisticsBipartiteUndirected.o $(SRCDIR)/graph.o $(SRCDIR)/utils.o $(SRCDIR)/loadGraph.o $(SRCDIR)/changeStatisticsUndirected.o
HASH_OBJS = $(OBJS:.o=_hash.o)
ARRAY_OBJS = $(OBJS:.o=_array.o)

all: testChangeStatsAttrBipartite testChangeStatsAttrBipartite_hash testChangeStatsAttrBipartite_array testBipartiteAlphaBetaChangeStats testBipartiteAlphaBetaChangeStats_hash testBipartiteAlphaBetaChangeStats_array


testChangeStatsAttrBipartite: testChangeStatsAttrBipartiteMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsAttrBipartite_hash: testChangeStatsAttrBipartiteMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsAttrBipartite_array: testChangeStatsAttrBipartiteMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm





testBipartiteAlphaBetaChangeStats: testBipartiteAlphaBetaChangeStatsMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testBipartiteAlphaBetaChangeStats_hash: testBipartiteAlphaBetaChangeStatsMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testBipartiteAlphaBetaChangeStats_array: testBipartiteAlphaBetaChangeStatsMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm



clean:
	$(RM) testChangeStatsAttrBipartite testChangeStatsAttrBipartite_hash
	$(RM) testChangeStatsAttrBipartiteMain.o testChangeStatsAttrBipartiteMain_hash.o
	$(RM) testChangeStatsAttrBipartite_array testChangeStatsAttrBipartiteMain_array.o
	$(RM) testBipartiteAlphaBetaChangeStats testBipartiteAlphaBetaChangeStats_hash
	$(RM) testBipartiteAlphaBetaChangeStatsMain.o testBipartiteAlphaBetaChangeStatsMain_hash.o
	$(RM) testBipartiteAlphaBetaChangeStats_array testBipartiteAlphaBetaChangeStatsMain_array.o


%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<

