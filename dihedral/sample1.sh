#!/bin/sh
#
#
com=../uv/bin/esipt.py
tg=substance_0403_4
smiles='O=C(C1=CC=CC=C1)C2=C(NC(N(CC3=CC=CC=C3)CC4=CC=CC=C4)=O)C=CC=C2'
opt="-phase 1"
python3 $com $tg $smiles $opt
