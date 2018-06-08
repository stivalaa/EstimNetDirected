#!/usr/bin/env python
# 
# run the example python script
#
import EstimNetDirectedSimpleDemo
from changeStatisticsDirected import *

EstimNetDirectedSimpleDemo.run_on_network_attr(
        'sample_statistics_n2000_directed_sim5000000.txt',
        [changeArc, changeReciprocity, changeAltInStars, changeAltOutStars],# changeAltKTrianglesT, changeAltTwoPathsTD,
        ["Arc", "Reciprocity", "AinS", "AoutS"] #, "AT-T", "A2P-TD"]
        )

