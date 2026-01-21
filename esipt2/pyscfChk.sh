#!/bin/sh
#
#
#python3 esipt.py o-cresol 'c1(C(C)(C)C)ccccc1O' -phase 1
com=../bin/esipt.py
ms=100

#v=2.8;pip install pyscf==$v;conda list > list.$v;cp AAP.txt AAP$v.txt
python3 $com AAP$v 'c1(C(=O)C)ccccc1N' -phase 2 -distance 10,19,0.98 -pt 1 -ms $ms
exit
v=2.9;pip install pyscf==$v;conda list > list.$v;cp AAP.txt AAP$v.txt
python3 $com AAP$v 'c1(C(=O)C)ccccc1N' -phase 2 -distance 10,19,0.98 -pt 1 -ms $ms
v=2.10;pip install pyscf==$v;conda list > list.$v;cp AAP.txt AAP$v.txt
python3 $com AAP$v 'c1(C(=O)C)ccccc1N' -phase 2 -distance 10,19,0.98 -pt 1 -ms $ms
v=2.11;pip install pyscf==$v;conda list > list.$v;cp AAP.txt AAP$v.txt
python3 $com AAP$v 'c1(C(=O)C)ccccc1N' -phase 2 -distance 10,19,0.98 -pt 1 -ms $ms
v=2.12;pip install pyscf==$v;conda list > list.$v;cp AAP.txt AAP$v.txt
python3 $com AAP$v 'c1(C(=O)C)ccccc1N' -phase 2 -distance 10,19,0.98 -pt 1 -ms $ms
