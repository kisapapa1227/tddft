import pandas as pd
import sys,os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pyscf
from pyscf.geomopt.geometric_solver import optimize
from pyscf import gto, scf, dft, tddft, lib
from pyscf.lib import param
import matplotlib.pyplot as plt
import dill
import datetime
import s00Utils
import math

MAX_CYCLE=20
MAX_CYCLE=50
MAX_CYCLE=200
MAX_LOOP=5
MAX_LOOP=1
NSTATES = 24

def getN(e,pe):
#
# e_min=1240.0/100.0/27.21162
#
# wavlen=1240/(mytd.e*27.21162)
#
    if len(pe)<1:
        return -1
    for n,ee in enumerate(pe):
        w1=1240.0/(e*27.21162)
        w2=1240.0/(ee*27.21162)
        rt=w1-w2
#        print("in getN",w1,w2)
#        rt=e/ee
        if rt>-5.0 and rt<5.0:
            return n
    return -1

def gaussBand(x, band, strength, stdev):
  constant=0.0003
  bandshape = constant * (strength/(1.0/stdev))*np.exp(-(((1.0/x)-(1.0/band))**2/(1.0/stdev)**2))
  return bandshape

def opt_sp_mm(mol, ff, numConfs=1000,pruneRmsThresh=0.5,numThreads=0,seed=1234):
    # 1. 1000個のコンフォマーを発生
    sm = Chem.MolToSmiles(mol)
    m = Chem.MolFromSmiles(sm)
    m_h = Chem.AddHs(m)
    cids = AllChem.EmbedMultipleConfs(m_h,
                                      numConfs=numConfs,
                                      randomSeed=seed,
                                      pruneRmsThresh=pruneRmsThresh,
                                      numThreads=numThreads)

    # 2,3. 各コンフォマーを最適化し，エネルギー計算
    energy = []
    if ff == 'uff':
        for cid in cids:
            uff = AllChem.UFFGetMoleculeForceField(m_h,
                                                   confId=cid)
            uff.Minimize()
            energy.append((uff.CalcEnergy(), cid))
    if ff == 'mmff':
        prop = AllChem.MMFFGetMoleculeProperties(m_h)
        for cid in cids:
            mmff = AllChem.MMFFGetMoleculeForceField(m_h, prop,
                                                     confId=cid)
            mmff.Minimize()
            energy.append((mmff.CalcEnergy(), cid))
    energy.sort()
    return m_h, [(i-energy[0][0],j) for i,j in energy]

class Log:
    def __init__(self,file):
        self.file=file+".log"
        with open(self.file,"w") as fp:
            fp.write("Start\n")
    def out(self,com):
        with open(self.file,"a") as fp:
            fp.write(com+"\n")
    def outs(self,pe,po):
        with open(self.file,"a") as fp:
            for ee,pp in zip(pe,po):
                fp.write(format(1240/(ee*27.21162),'.2f')+":"+format(pp,'.3f')+"\n")

def mkBlock(uff_mol,confId=0):
    mol_block0=Chem.MolToMolBlock(uff_mol,confId=confId)
    mol_block1=mol_block0.split("\n")[4:]
    mol_block2=[]
    for l in mol_block1:
        arr=l.split()
        if len(arr)>4 and not 'M' in arr[0]:
            o=[arr[3],arr[0],arr[1],arr[2]]
            mol_block2.append(o)
    mol_block3=""
    for line in mol_block2:
        mol_block3+="  ".join(line)
        mol_block3+="\n"
    return mol_block3

def elaps(tm1,tm2):
    dy1=tm1.split(' ')[0]
    hr1=tm1.split(' ')[1].split(':')

    dy2=tm2.split(' ')[0]
    hr2=tm2.split(' ')[1].split(':')

    net=(int(dy2)-int(dy1))*24*60
    net+=(int(hr2[0])-int(hr1[0]))*60
    net+=(int(hr2[1])-int(hr1[1]))

    return str(net)

def conv(ss):
    xyz='';atom='';i=0
    for s in ss:
        if s!='':
            if i<3:
#                xyz=xyz+' '+s
                xyz=xyz+' {:5.3f}'.format(float(s))
#                xyz=xyz+' {:5.3f}'.format(float(s)/param.BOHR)
            else:
                return s+xyz
            i=i+1
    return xyz

def listToStr(ss):
    ret=""
    for s in ss:
        ret=ret+s[0]+" "+str(s[1][0]*param.BOHR)+" "+\
                str(s[1][1]*param.BOHR)+" "+str(s[1][2]*param.BOHR)+"\n "

    return ret

def reshape(dd):
    r=[]
    for d in dd:
        if len(d)==3:
            r.append([int(d[0])-1,int(d[1])-1,float(d[2])])
        elif len(d)==5:
            r.append([int(d[0])-1,int(d[1])-1,int(d[2])-1,int(d[3])-1,float(d[4])])
    return r

def setDihe(dihe,fn):
    if len(dihe)<1:
        return
    ff=open(fn,"a")
    for d in dihe:
        ff.write("$set\ndihedral "+str(d[0]+1)+" "+str(d[1]+1)+" "+str(d[2]+1)+" "+str(d[3]+1)+" "+str(d[4])+"\n")
    ff.close()

def setDist(dist,mol_block,fn):
    ff=open(fn,"w")
    if len(dist)<1:
        ff.close()
        return mol_block
    mb="";
    mm=mol_block.split(";")
#    print(mm)
    for d in dist:
        print(d[0],d[1]-1)
        p0=mm[d[0]].split(' ')
        p1=mm[d[1]].split(' ')
        x0,y0,z0=float(p0[1]),float(p0[2]),float(p0[3]),
        x1,y1,z1=float(p1[1]),float(p1[2]),float(p1[3]),
        l=math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))
        x1,y1,z1=x0+(x1-x0)*d[2]/l,y0+(y1-y0)*d[2]/l,z0+(z1-z0)*d[2]/l;
        l2=math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))
#        print(f"length {l} to {l2}")
        mm[d[1]]=" ".join([p1[0],str(x1),str(y1),str(z1)])
        ff.write("$freeze\ndistance "+str(d[0]+1)+" "+str(d[1]+1)+"\n")
    ff.close()
    return ";".join(mm)
#    print(mm)

def toArray(mb):
    ss=[];s="";f=False
    for i in mb:
        if ord(i)==10:
            ss.append(s);s=""
        else:
            s=s+i
    return ss

def toBlock(ss,xyz):
    sdf="";c=0
    for i in range(4):
        sdf=sdf+ss[i]+chr(10);c=c+1
    for i in range(len(xyz)):
        l=ss[i+4].split('.')[3].split('0000')[1]
        sdf=sdf+align(xyz[i])+l+chr(10);c=c+1
    for i in range(c,len(ss)):
        sdf=sdf+ss[i]+chr(10)
    return sdf

def align(pts):
    ret=""
    for i in range(3):
        if pts[i]<0.0:
            ret=ret+"   {:.4f}".format(pts[i])
        else:
            ret=ret+"    {:.4f}".format(pts[i])
    return ret

def molBlockFromFile(fn):
    print("in molBlockFromFile : read from ",fn)
    atom=""
    with open(fn) as ff:
        for i,l in enumerate(ff):
            if i<4:
                continue
            s=l.split('\n')[0].split(' ')
            if len(s)<12:
                continue
            c=conv(s)
            if atom=="":
                atom=c
            else:
                atom=atom+";"+c

    return atom

#
# main
#
smi='c1cccc2c1c(CCN3C(=O)C(CC(=O)NCC(C)C)NC(=O)3)cN2'
basis="sto-3g"
basis="def2-qzvp"
basis="cc-pvdz"
basis="631g" # the original
basis="621g"
basis="6311g"

ms=3

if len(sys.argv)<2:
    print(sys.argv[0],"for large moleclues . i.e) complex of 3 cyclic structure")
    print("default:MAX_CYCLE:",MAX_CYCLE," MAX_LOOP", MAX_LOOP,"NSTATES",NSTATES);
    print("specify:-mc MAX_CYCLE, -ml MAX_LOOP -ns NSTATES -bs basis -wo (withoutSolvent) -ms maxstep(for optimize)")
    exit()

eps=True
opt=True
flag=False
pt="o"
dist=[];dihe=[]
for i in range(1,len(sys.argv)):
    if flag:
        flag=False
        continue
    if sys.argv[i]=='-mc':
        MAX_CYCLE=int(sys.argv[i+1])
        flag=True
    elif sys.argv[i]=='-ml':
        MAX_LOOP=int(sys.argv[i+1])
        flag=True
    elif sys.argv[i]=='-ns':
        NSTATES=int(sys.argv[i+1])
        flag=True
    elif sys.argv[i]=='-bs':
        basis=sys.argv[i+1]
        flag=True
    elif sys.argv[i]=='-wo':
        eps=False
    elif sys.argv[i]=='-no-op':
        opt=False
    elif sys.argv[i]=='-pt':
        pt=sys.argv[i+1].zfill(2)
    elif sys.argv[i]=='-ms':
        ms=int(sys.argv[i+1])
    elif sys.argv[i]=='-phase':
        phase=int(sys.argv[i+1])
    elif sys.argv[i]=='-distance':
        dist.append(sys.argv[i+1].split(","))
    elif sys.argv[i]=='-dihedral':
        dihe.append(sys.argv[i+1].split(","))

if phase > 1 and len(dist)<1 and len(dihe)<1:
    print("must specify constrict type -distance/-dihedral")
    exit()
dist=reshape(dist)
dihe=reshape(dihe)

xc = 'b3lyp'
xc = 'CAMB3LYP'
hartree_to_kcal = 627.509
smi=sys.argv[2]

print("default:MAX_CYCLE:",MAX_CYCLE," MAX_LOOP", MAX_LOOP,"NSTATES",NSTATES);

mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
m_h=Chem.AddHs(mol)

dt_start=datetime.datetime.now()
init=1235

spin=0

if phase==1:
    uff_mol, energy_cid = opt_sp_mm(mol, 'uff',seed=init)
    mol_block=mkBlock(uff_mol,confId=0)
    mol3=gto.Mole(atom=mol_block, basis=basis).build() #1
    mf = mol3.RKS(xc=xc)#2#1
    if eps:
        mf=mf.PCM()
        mf.with_solvent.method='IEF-PCM'
        mf.with_solvent.eps=4.7113 #chloroform
    mol_eq = optimize(mf, maxsteps=ms)
    gn=sys.argv[1]
    s00Utils.molBlockToFiles(smi,mol_eq,fn=gn+'.txt')
    exit(0)

if phase==2 or phase==3:
    gn=sys.argv[1];fn=gn+"_"+pt
    log=Log(fn)

    hx=gn+"_"+str(int(pt)-1)
    if os.path.isfile(hx+".txt"):
        mol_block=molBlockFromFile(hx+".txt")
    else:
        mol_block=molBlockFromFile(gn+".txt")
    mol_block=setDist(dist,mol_block,fn+"_c.txt")
    setDihe(dihe,fn+"_c.txt")

    mol3=gto.Mole(atom=mol_block, basis=basis).build() #1
    mf = mol3.RKS(xc=xc)#2#1
    if eps:
        mf=mf.PCM()
        mf.with_solvent.method='IEF-PCM'
        mf.with_solvent.eps=4.7113 #chloroform
    if opt:
        params={"constraints":fn+"_c.txt"}
        mol_eq = optimize(mf, **params,maxsteps=ms)
    else:
        mol_eq=mol3
    s00Utils.molBlockToFiles(smi,mol_eq,fn=fn+'.txt')

mf=dft.RKS(mol_eq)
mf.xc = xc
mf.kernel()
e_tot=mf.e_tot

homo_index = (mf.mo_occ > 0).nonzero()[0][-1]
lumo_index = (mf.mo_occ == 0).nonzero()[0][0]

hartree_to_ev = 27.2114
homo_energy = mf.mo_energy[homo_index] * hartree_to_ev
lumo_energy = mf.mo_energy[lumo_index] * hartree_to_ev
lumo1_energy = mf.mo_energy[lumo_index+1] * hartree_to_ev

log.out(f'HOMO Energy: {homo_energy:.6f} eV')
log.out(f'LUMO Energy: {lumo_energy:.6f} eV')
log.out(f'LUMO+1 Energy: {lumo1_energy:.6f} eV')

#for d in dir(mf):
#    print(d)
#print(mf.e_tot,mf.energy_tot())
#print(mf.energy_elec(),mf.energy_nuc())
#for k in mf.scf_summary.keys():
#    print(k,mf.scf_summary[k])


if phase==3:
    if len(dist)>0:
        s00Utils.writeEnergy(mol_block,mf,sys.argv,dist[0][2])
    else:
        s00Utils.writeEnergy(mol_block,mf,sys.argv,dihe[0][4])
    exit()

log.out("tddft statt")

mytd = tddft.TDDFT(mf)
mytd.nstates = NSTATES

#
# wavlen=1240/(mytd.e*27.21162)
#
e_max=1240.0/100.0/27.21162
e_min=1240.0/400.0/27.21162

noise=1e-4;loop=1;pc=[];pe=[];po=[]
init=False
init=True

while True:
    try:
        if init:
            x0=mytd.init_guess(mf)
            noise=1e-4
            for x in x0:
                x+=np.random.uniform(-noise,noise,x.size)
    except:
        init=False
    if not init:
        log.out("mytd.init_guess error: never mind")

    print("Start ",loop)
#    mytd.kernel(v2a=x0)
    mytd.max_cycle=MAX_CYCLE
    try:
        log.out(f"mytd.kernel: start {loop}")
        if init:
            mytd.kernel(x0=x0)
        else:
            mytd.kernel()
    except:
        log.out("mytd.kernel error -> mytd.nstates = 8")
        mytd.nstates = 8
        if init:
            mytd.kernel(x0=x0)
        else:
            mytd.kernel()

    print(mytd.converged)
#    c=mytd.converged
#    o=mytd.oscillator_strength()
    cnt=0
    for e,c,o in zip(mytd.e,mytd.converged,mytd.oscillator_strength()):
#        if e<e_min or e>e_max:
#            continue
        print(e,o,c)
        if c:
            cnt+=1
            n=getN(e,pe)
            if n<0:
                pe.append(e)
                po.append(o)
            else:
                if po[n]<o:
                    po[n]=o
    pc.append(len(pe))
    log.out(f"number of converged points in loop {loop}: {cnt}")
    log.outs(pe,po)
    for ee in pe:
        print(1240/(ee*27.21162),end=",")
    print()
    for oo in po:
        print(oo,end=",")
    print()
    print(loop,pc)
    loop+=1
    if loop>MAX_LOOP:
        break
mytd.analyze()

print("tddft completed")
dt_end=datetime.datetime.now()

with open(fn+".tdd","w") as fp:
    wavlen=1240/(np.array(pe)*27.21162)
    for f in wavlen:
        fp.write(str(f)+",")
    fp.write("\n")
    for f in po:
        fp.write(str(f)+",")
    fp.write("\n"+smi+"\nseed=")
    fp.write(f"\nbasis={basis}\nxc={xc}\nspin={spin}\n")
    fp.write(f"max_cycle={MAX_CYCLE},max_loop={MAX_LOOP},nstates={mytd.nstates}\n")
    for f in pc:
        fp.write(str(f)+",")
    fp.write(f'\nSCF energy: {e_tot} hartree\n')
    fp.write(f'HOMO Energy: {homo_energy:.6f} eV\n')
    fp.write(f'LUMO Energy: {lumo_energy:.6f} eV\n')
    fp.write(f'LUMO+1 Energy: {lumo1_energy:.6f} eV\n')
    if eps:
        fp.write("PCM : chloroform\n")
    tm1=dt_start.strftime('\n%Y%m%d %H:%M:%S')
    tm2=dt_end.strftime('%Y%m%d %H:%M:%S')
    fp.write(tm1+"-"+tm2+"-->"+elaps(tm1,tm2)+"\n")
 
if len(dist)>0:
    s00Utils.writeEnergy(mol_block,mf,sys.argv,dist[0][2])
else:
    s00Utils.writeEnergy(mol_block,mf,sys.argv,dihe[0][4])
