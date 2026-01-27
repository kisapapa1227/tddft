from rdkit import Chem
from rdkit.Chem import AllChem
import random
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
from rdkit import DataStructs
import s00Utils,sys
import itertools
# 分子リストの例

def myRep(i,n):
    return i.split(':')[0]+':'+str(n)+']'

def prep2(m2,smarts,insert):
    s=[]
    for x in smarts:
        s.append(Chem.MolToSmiles(x))

    dst=[]
    for x in itertools.permutations(s,insert):
        d=m2
        for n,i in enumerate(x):
            d=d+"."+myRep(i,n+1)
        print("molzip",d)
        dst.append(Chem.molzip(Chem.MolFromSmiles(d)))

    src=[[],[]];chk=[];cnt=2
    for i in dst:
        x=Chem.MolToSmiles(i)
        print(x)
        if not x in chk:
            chk.append(x)
            src[0].append(i)
            src[1].append(str(cnt).zfill(2))
            cnt=cnt+1

    return src

def prep1(m2,smarts):
    src=[[m2],[""]]
    for i,frag in enumerate(smarts):
        dst=[[],[]]
        for s1,s2 in zip(src[0],src[1]):
            for j,f in enumerate(frag):
                p=Chem.MolToSmiles(f)
                print(i,j,f)
                dst[0].append(s1+"."+p)
                dst[1].append(s2+str(j+2).zfill(2))
        src=dst

    dst=[[],[]]
    for s1, s2 in zip(src[0],src[1]):
        print("molzip",s1)
        dst[0].append(Chem.molzip(Chem.MolFromSmiles(s1)))
        dst[1].append(s2)
    return dst

fno=sys.argv[1].split(".")[0]
with open(sys.argv[1],"r") as fp:
    init=fp.readline().split('\n')[0]
    m1=fp.readline().split('\n')[0]
    m2=fp.readline().split('\n')[0]
    lines=fp.readlines()

smarts=[]
dummy='COCOCO'
for i,line in enumerate(lines):
    s=[]
    for l in line.split('\n')[0].split(","):
        if l=='H':
            p=Chem.MolFromSmarts('[*:'+str(i+1)+']'+dummy)
        else:
            p=Chem.MolFromSmarts('[*:'+str(i+1)+']'+l)
        s.append(p)
        print(p)
#        s.append(Chem.MolFromSmiles(Chem.MolToSmiles(p)))
        print(i,Chem.MolToSmiles(p))
    smarts.append(s)

import numpy as np

name=[];c=0;box=[];page=1
for i,frag in enumerate(smarts):
    for j,s in enumerate(frag):
        csmiles=Chem.MolToSmarts(s)
        print(csmiles)
        if dummy in csmiles:
            csmiles='None'
            box.append(Chem.MolFromSmiles('Cl'))
        else:
            box.append(s)
        name.append(str(i+1).zfill(2)+"x"+str(j+1).zfill(2)+":"+csmiles)
        c=c+1
#  name.append(str(c).zfill(2)+":"+csmiles)
        if c%(36)==0:
            img = Draw.MolsToGridImage(box, legends=name, molsPerRow=6, subImgSize=(200,150))  # Set molsPerRow and subImgSize to desired values
            s00Utils.saveAs(img,fn=fno+"_A"+str(page)+".png")
            box=[]
            name=[]
            page+=1

if len(box)>0:
     img = Draw.MolsToGridImage(box, legends=name, molsPerRow=6, subImgSize=(200,150))  # Set molsPerRow and subImgSize to desired values
     s00Utils.saveAs(img,fn=fno+"_A"+str(page)+".png")
print("num:",c)

fp=open(fno+".smi","w")
index=""
for i,frag in enumerate(smarts):
    index=index+"01"
theName=init+index
theMol=Chem.MolFromSmiles(m1)
fp.write(theName+" '"+Chem.MolToSmiles(theMol)+"'\n")

print("-----------",len(smarts),m2.count("*:"))

insert=m2.count("*:")
if len(smarts)==insert:
    src=prep1(m2,smarts)
else:
    adapter=smarts[0]
    print(smarts)
    print("<--",smarts[0])
    if len(adapter)<insert:
        print("number of smarts is less than that of insertion point")
        exit()
    else:
        src=prep2(m2,adapter,insert)
        for s1,s2 in zip(src[0],src[1]):
            print("-",s1,"-",s2)

dum=Chem.MolFromSmiles(dummy) 
cnt=1;name=[];loop=0;dst=[]
name.append(theName);dst.append(theMol)
for s1,s2 in zip(src[0],src[1]):
  x=s1.GetSubstructMatch(Chem.MolFromSmiles(dummy)) 
  if len(x)>0:
#      print(x)
      s1=AllChem.DeleteSubstructs(s1,dum)
#  new_mol=Chem.MolFromSmiles(s1)
  fp.write(init+s2+" '"+Chem.MolToSmiles(s1)+"'\n")
  name.append(init+s2)
  dst.append(s1)
  cnt+=1
  if cnt==20:
      img = Draw.MolsToGridImage(dst,legends=name,molsPerRow=4, subImgSize=(200,200))  # Set molsPerRow and subImgSize to desired values
      s00Utils.saveAs(img,fn=fno+"_"+str(loop)+".png")
      dst=[];name=[];loop+=1
      cnt=0
      if loop>10:
          break

if len(dst)>0:
  img = Draw.MolsToGridImage(dst, legends=name,molsPerRow=4, subImgSize=(200,200))  # Set molsPerRow and subImgSize to desired values
  s00Utils.saveAs(img,fn=fno+"_"+str(loop)+".png")

fp.close()
