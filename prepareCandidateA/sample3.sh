#!/bin/sh
#
#
sh ../bin/mkPrmForTdds2svg.sh > plot.prm
python3 ../bin/tdds2svg.py plot.prm -up 250 -r 0.1
