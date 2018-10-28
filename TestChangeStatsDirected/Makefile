include ../src/common.mk
-include ../src/local.mk

CC = gcc

CFLAGS += -I../src

OBJS =  ../src/changeStatisticsDirected.o ../src/digraph.o ../src/utils.o ../src/changeStatisticsDirected.o
HASH_OBJS = $(OBJS:.o=_hash.o)

all: testChangeStatsDirected testChangeStatsDirected_hash


testChangeStatsDirected: testChangeStatsDirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsDirected_hash: testChangeStatsDirectedMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	$(RM) testChangeStatsDirected testChangeStatsDirected_hash
	$(RM) testChangeStatsDirectedMain.o testChangeStatsDirectedMain_hash.o

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_HASHTABLES -c -o $@ $<
