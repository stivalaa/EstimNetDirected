include ../src/common.mk
-include ../src/local.mk

CC = gcc

CFLAGS += -I../src

OBJS =  ../src/changeStatisticsDirected.o ../src/digraph.o ../src/utils.o ../src/changeStatisticsDirected.o ../src/loadDigraph.o
HASH_OBJS = $(OBJS:.o=_hash.o)

all: testChangeStatsDirected testChangeStatsDirected_hash testSetFunctions


testChangeStatsDirected: testChangeStatsDirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsDirected_hash: testChangeStatsDirectedMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm


testSetFunctions: testSetFunctions.o  $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm


clean:
	$(RM) testChangeStatsDirected testChangeStatsDirected_hash
	$(RM) testChangeStatsDirectedMain.o testChangeStatsDirectedMain_hash.o
	$(RM) testSetFunctions

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_HASHTABLES -c -o $@ $<
