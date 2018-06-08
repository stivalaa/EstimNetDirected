include ../src/common.mk
-include ../src/local.mk

CC = gcc

CFLAGS += -I../src

OBJS =  ../src/changeStatisticsDirected.o ../src/digraph.o ../src/utils.o ../src/changeStatisticsDirected.o

all: testChangeStatsDirected


testChangeStatsDirected: testChangeStatsDirectedMain.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	$(RM) testChangeStatsDirected
	$(RM) testChangeStatsDirectedMain.o

