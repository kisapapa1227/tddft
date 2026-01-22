#!/bin/sh
#
#
com=~/uv/bin/esipt.py
tg=substance_0403_4
smiles='O=C(C1=CC=CC=C1)C2=C(NC(N(CC3=CC=CC=C3)CC4=CC=CC=C4)=O)C=CC=C2'
opt="-phase 3"
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-133 -pt 1
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-113 -pt 2
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-93 -pt 3
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-73 -pt 4
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-53 -pt 5
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-33 -pt 6
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,-13 -pt 7
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,7 -pt 8
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,27 -pt 9
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,47 -pt 10
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,67 -pt 11
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,87 -pt 12
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,107 -pt 13
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,127 -pt 14
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,147 -pt 15
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,167 -pt 16
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,187 -pt 17
python3 $com $tg $smiles $opt -dihedral 1,2,9,10,207 -pt 18
