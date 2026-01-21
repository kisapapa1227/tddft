#!/bin/sh
#
#
python3 ../bin/mkScript.py molecularBackbone.smi > runAll.sh
sh runAll.sh
