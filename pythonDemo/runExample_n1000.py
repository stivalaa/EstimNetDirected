#!/usr/bin/env python
# 
# run the example python script
#
import EstimNetDirectedSimpleDemo
from changeStatisticsDirected import *

EstimNetDirectedSimpleDemo.run_on_network_attr(
        'sample_statistics_n1000_directed_binattr_sim620000000.txt',
        [changeArc, changeReciprocity, changeAltInStars, changeAltOutStars,
         changeAltKTrianglesT, changeAltTwoPathsTD,
         changeReceiver, changeSender, changeInteraction],
        ["Arc", "Reciprocity", "AinS", "AoutS", "AT-T", "A2P-TD",
         "Receiver", "Sender", "Interaction"],
        'binaryAttributes_50_50_n1000.txt')

