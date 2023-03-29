include ../../src/common.mk
-include ../../src/local.mk

SRCDIR = ../../src

CC = gcc

CFLAGS += -I$(SRCDIR)

OBJS = $(SRCDIR)/changeStatisticsGeneral.o $(SRCDIR)/graph.o $(SRCDIR)/utils.o  $(SRCDIR)/loadGraph.o $(SRCDIR)/changeStatisticsBipartiteUndirected.o $(SRCDIR)/changeStatisticsUndirected.o
HASH_OBJS = $(OBJS:.o=_hash.o)
ARRAY_OBJS = $(OBJS:.o=_array.o)

all: testRandom


testRandom: testRandom.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm




clean:
	$(RM) testRandom
	$(RM) testRandom.o

%_hash.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -DTWOPATH_HASHTABLES -c -o $@ $<

%_array.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTWOPATH_LOOKUP -c -o $@ $<
