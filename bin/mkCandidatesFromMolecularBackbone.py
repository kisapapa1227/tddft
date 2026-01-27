from rdkit import Chem
from rdkit.Chem import AllChem
import random
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
from rdkit import DataStructs
import s00Utils,sys
# 分子リストの例

fno=sys.argv[1].split(".")[0]
with open(sys.argv[1],"r") as fp:
    init=fp.readline().split('\n')[0]
    m1=fp.readline().split('\n')[0]
    m2=fp.readline().split('\n')[0]
    weight=fp.readline().split('\n')[0].split(':')[1].split(',')
    wk=fp.readline().split('\n')[0].split(':')[1].split(',')
    AtomNum=[int(wk[0]),int(wk[1])+1]
#    withAtom=fp.readline().split('\n')[0].split(':')[1].split(',')

suppl = Chem.SDMolSupplier('../bin/solubility.train.sdf')
molecules = [mol for mol in suppl if mol is not None]

molecule_list = molecules
fragment_list = []
# 分子をフラグメントに分割
for mol in molecule_list:
    core = MurckoScaffold.GetScaffoldForMol(mol)
    tmp=Chem.ReplaceCore(mol,core,labelByIndex=True)
    if tmp is not None:
      frag=Chem.GetMolFrags(tmp,asMols=True)
      fragment_list.extend(frag)

import numpy as np

new_fragment_list=[]
for frag in fragment_list:
  dummy = Chem.MolFromSmiles('[*]')
  dummy1 = Chem.MolFromSmiles('[*:1]')
  frag_ = AllChem.ReplaceSubstructs(frag, dummy, dummy1)[0]
  new_fragment_list.append(frag_)
unique_frags=np.unique([Chem.MolToSmarts(frag) for frag in new_fragment_list ]).tolist()

print(len(unique_frags))
limit=[]
for frag in unique_frags:
  m3 = Chem.MolFromSmarts(frag)
  num=m3.GetNumAtoms()
  if num>AtomNum[0] and num<AtomNum[1]:
      limit.append(frag)

unique_frags=limit

name=[];tmp_smart=[];why=[];fine_frags=[]
c=0;low_per=6;col_per=6;page=1
for frag in unique_frags:
  m3 = Chem.MolFromSmarts(frag)
  csmiles=Chem.MolToSmiles(m3)
  if csmiles in why:
      continue
  fine_frags.append(frag)
  c+=1
  why.append(csmiles)
  tmp_smart.append(m3)
#  csmiles=Chem.MolToSmarts(m3)
  name.append(str(c).zfill(3)+":"+csmiles)
#  name.append(str(c).zfill(2)+":"+csmiles)
  if c%(low_per*col_per)==0:
     img = Draw.MolsToGridImage(tmp_smart, legends=name, molsPerRow=col_per, subImgSize=(200,150))  # Set molsPerRow and subImgSize to desired values
     s00Utils.saveAs(img,fn=fno+"_A"+str(page)+".png")
     tmp_smart=[]
     name=[]
     page+=1

if len(tmp_smart)>0:
     img = Draw.MolsToGridImage(tmp_smart, legends=name, molsPerRow=col_per, subImgSize=(200,150))  # Set molsPerRow and subImgSize to desired values
     s00Utils.saveAs(img,fn=fno+"_A"+str(page)+".png")

fp=open(fno+".smi","w")
theName=init+'000'
theMol=Chem.MolFromSmiles(m1)
fp.write(theName+" '"+Chem.MolToSmiles(theMol)+"'\n")

insert=m2.count("*:")
dst=[];name=[];cnt=0;c=1;loop=1
for frag in why:
    adapter=m2
    s2=str(c).zfill(3)
    print("init",s2)
    for i in range(insert):
        adapter=adapter+"."+frag.split("[*")[0]+"[*:"+str(i+1)+"]"
    print(adapter)
    try:
        s1=Chem.molzip(Chem.MolFromSmiles(adapter))
    except:
        print("Error",adapter)
    else:
        dst.append(s1)
        fp.write(init+s2+" '"+Chem.MolToSmiles(s1)+"'\n")
        name.append(init+s2)
        cnt+=1;c=c+1
#        print(c,"<--")
    if cnt==20:
        print("name",name)
        img = Draw.MolsToGridImage(dst,legends=name,molsPerRow=4, subImgSize=(200,200))  # Set molsPerRow and subImgSize to desired values
        s00Utils.saveAs(img,fn=fno+"_"+str(loop)+".png")
        dst=[];name=[];loop+=1
        cnt=0

if len(dst)>0:
    img = Draw.MolsToGridImage(dst, legends=name,molsPerRow=4, subImgSize=(200,200))  # Set molsPerRow and subImgSize to desired values
    s00Utils.saveAs(img,fn=fno+"_"+str(loop)+".png")
fp.close()
print("num:",c,len(why))
