#!/bin/sh
#
#
python3 ../bin/tddftSolver.py Indigo_pcm 'C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O' 
python3 ../bin/tddftSolver.py Indigo 'C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O' -wo

