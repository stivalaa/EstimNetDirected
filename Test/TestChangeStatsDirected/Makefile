include ../src/common.mk
-include ../src/local.mk

CC = gcc

CFLAGS += -I../src

OBJS = ../src/changeStatisticsGeneral.o  ../src/changeStatisticsDirected.o ../src/graph.o ../src/utils.o ../src/changeStatisticsDirected.o ../src/loadGraph.o
HASH_OBJS = $(OBJS:.o=_hash.o)
ARRAY_OBJS = $(OBJS:.o=_array.o)

all: testChangeStatsDirected testChangeStatsDirected_hash testSetFunctions testChangeStatsDirected_array


testChangeStatsDirected: testChangeStatsDirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsDirected_hash: testChangeStatsDirectedMain_hash.o $(HASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

testChangeStatsDirected_array: testChangeStatsDirectedMain_array.o $(ARRAY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm


testSetFunctions: testSetFunctions.o  $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm


clean:
	$(RM) testChangeStatsDirected testChangeStatsDirected_hash
	$(RM) testChangeStatsDirectedMain.o testChangeStatsDirectedMain_hash.o
	$(RM) testSetFunctions
	$(RM) testChangeStatsDirected_array testChangeStatsDirectedMain_array.o

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<
