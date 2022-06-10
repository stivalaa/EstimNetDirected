
SUBDIRS = TestChangeStatsDirected TestChangeStatsUndirected TestSetFunctions TestChangeStatsBipartite TestChangeStatsAttrBipartite TestRandom

all:
	for dir in $(SUBDIRS); do \
		(cd $$dir && $(MAKE)) \
	done

clean:
	for dir in $(SUBDIRS); do \
		(cd $$dir && $(MAKE) clean) \
	done
