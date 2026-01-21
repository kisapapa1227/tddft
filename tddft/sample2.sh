#!/bin/sh
#
#
python3 ../bin/preparePlotParemeterFile.py > plot.prm
python3 ../bin/makeTddftReport.py plot.prm
