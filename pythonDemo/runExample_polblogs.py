#!/usr/bin/env python
# 
# run the example python script
#
import EstimNetDirectedSimpleDemo
from changeStatisticsDirected import *

EstimNetDirectedSimpleDemo.run_on_network_attr(
    'polblogs/polblogs_arclist.txt',
    [changeArc, changeReciprocity, changeAltInStars,
     changeAltOutStars, changeAltTwoPathsTD,
     changeAltKTrianglesT, changeAltKTrianglesC,
     changeMismatching, changeMismatchingReciprocity],
    ["Arc", "Reciprocity", "AinS", "AoutS", "A2P-TD", "AKT-T", "AKT-C",
     "Mismatching", "MismatchingReciprocity"],
    binattr_filename=None,
    catattr_filename='polblogs/polblogs_catattr.txt')

