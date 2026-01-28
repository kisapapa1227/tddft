import pandas as pd
import sys,os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pyscf
from pyscf.geomopt.geometric_solver import optimize
from pyscf import gto, scf, dft, tddft
import matplotlib.pyplot as plt
import datetime
import s00Utils,dill

MAX_CYCLE=20
MAX_CYCLE=50
MAX_CYCLE=200

def align(pts):
    ret=""
    for i in range(3):
        if pts[i]<0.0:
            ret=ret+"   {:.4f}".format(pts[i])
        else:
            ret=ret+"    {:.4f}".format(pts[i])
    return ret

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

def toAtom(mb):
    atom=''
    for i in mb:
        x=i.split(',')
        if len(x)<13:
            continue
        if atom=='':
            atom=x[4]+' '+x[1]+' '+x[2]+' '+x[3]
        else:
            atom=atom+';'+x[4]+' '+x[1]+' '+x[2]+' '+x[3]
    print("atom",atom)
    return atom

def toAtomNg(mb):
    atom="";ss=[];s="";f=False

    print("in toAtom",mb)
    for i in mb:
        if i=='':
            continue
        print(i,type(i))
        if ord(i)==10:
            ss.append(s);s="";f=False
        else:
            if i!=' ':
                if f:
                    s=s+','
                s=s+i
                f=False
            else:
                f=True
    #
    for i in range(4,len(ss)):
        x=ss[i].split(',')
        if len(x)<10:
            return atom
        if atom=='':
            atom=x[4]+' '+x[1]+' '+x[2]+' '+x[3]
        else:
            atom=atom+';'+x[4]+' '+x[1]+' '+x[2]+' '+x[3]
    return atom

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

def toList(mb):
    ss=[];s="";f=False
    for i in mb:
        if ord(i)==10:
            ss.append(s);s="";f=False
        else:
            if i!=' ':
                if f:
                    s=s+','
                s=s+i
                f=False
            else:
                f=True
    return ss

def toArray(mb):
    ss=[];s="";f=False
    for i in mb:
        if ord(i)==10:
            ss.append(s);s=""
        else:
            s=s+i
    return ss
#
# main
#
smi='c1cccc2c1c(CCN3C(=O)C(CC(=O)NCC(C)C)NC(=O)3)cN2'
basis="sto-3g"
basis="def2-qzvp"
basis="cc-pvdz"
basis="631g"
basis="621g"
basis="6311g"
eps=True
flag=False
EPS=4.7113 #chloroform

ms=3

if len(sys.argv)<2:
    print(sys.argv[0],"for large moleclues . i.e) complex of 3 cyclic structure")
    print("default:MAX_CYCLE:",MAX_CYCLE);
    print("specify:-mc MAX_CYCLE, -bs basis -wo (withoutSolvent)")
    exit()

for i in range(3,len(sys.argv)):
    if flag:
        flag=False
        continue
    if sys.argv[i]=='-mc':
        MAX_CYCLE=int(sys.argv[i+1])
        flag=True
    elif sys.argv[i]=='-eps':
        EPS=sys.argv[i+1]
        flag=True
    elif sys.argv[i]=='-bs':
        basis=sys.argv[i+1]
        flag=True
    elif sys.argv[i]=='-wo':
        eps=False

fn=sys.argv[1];smi=sys.argv[2]
fn=fn+".erg"

print("default:MAX_CYCLE:",MAX_CYCLE,fn);

if os.path.isfile(fn):
    print(fn+" exists.")
    exit()

mol = Chem.MolFromSmiles(smi)
print("ok",smi)

dt_start=datetime.datetime.now()
loop=1
init=1235
try_init=[]

xc = 'b3lyp'
xc = 'CAMB3LYP'

log=Log(fn)

while True:
    spin=0
    uff_mol, energy_cid = opt_sp_mm(mol, 'uff',seed=init)
    try_init.append(init)
#    mol_block=mkBlock(uff_mol,confId=0)
    mol_block=mkBlock(uff_mol,confId=0)
    try:
        mol3 = gto.Mole(atom=mol_block, basis=basis).build() #1
    except:
        print("--- alternative")
        try:
            spin=1
            mol3 = gto.Mole(atom=mol_block,basis=basis,spin=spin,charge=0).build()
        except:
            print("--- another alternative")
            try:
                log.out("change spin=2")
                spin=2
                mol3 = gto.Mole(atom=mol_block,basis=basis,spin=spin,charge=0).build()
            except:
                log.out("No way")
                print("--- No way")
                exit()
    mf = mol3.RKS(xc=xc)#2#1
    if eps:
        mf=mf.PCM()
        mf.with_solvent.method='IEF-PCM'
        mf.with_solvent.eps=EPS #chloroform

    try:
        mol_eq = optimize(mf, maxsteps=ms)
        print("->",mol_eq.atom_coords(),"<-")
    except: 
        log.out("Optimizer error!!")
        print("optimizer error!!\nNever mind")
        if init > 10:
            exit()
        loop+=1
        init+=1
        continue

    rr=[]

    mf = dft.RKS(mol_eq)
    mf.xc = xc
    if eps:
        mf=mf.PCM()
        mf.with_solvent.method='IEF-PCM'
        mf.with_solvent.eps=EPS #chloroform
    mf.kernel()

    print("mf.converged<---",mf.converged,loop)
    if mf.converged or loop>5:
        if mf.converged:
            log.out("mf.converged at "+str(loop)+" loop")
        else:
            log.out("mf.converged : no way")
        try_init.append(mf.converged)
        break 
    loop+=1
    init+=1

homo_index = (mf.mo_occ > 0).nonzero()[0][-1]
lumo_index = (mf.mo_occ == 0).nonzero()[0][0]

hartree_to_ev = 27.2114
homo_energy = mf.mo_energy[homo_index] * hartree_to_ev
lumo_energy = mf.mo_energy[lumo_index] * hartree_to_ev
lumo1_energy = mf.mo_energy[lumo_index+1] * hartree_to_ev

log.out(f'HOMO Energy: {homo_energy:.6f} eV')
log.out(f'LUMO Energy: {lumo_energy:.6f} eV')
log.out(f'LUMO+1 Energy: {lumo1_energy:.6f} eV')

dt_end=datetime.datetime.now()

with open(fn,"w") as fp:
    fp.write(smi+"\nseed=")
    for l in try_init:
        fp.write(str(l)+",")
    fp.write(f"\nbasis={basis}\nxc={xc}\nspin={spin}\n")
    fp.write(f'\nHOMO Energy: {homo_energy:.6f} eV\n')
    fp.write(f'LUMO Energy: {lumo_energy:.6f} eV\n')
    fp.write(f'LUMO+1 Energy: {lumo1_energy:.6f} eV\n')
    if eps:
        fp.write(f'PCM : epsilon {EPS}\n')
    tm1=dt_start.strftime('\n%Y%m%d %H:%M:%S')
    tm2=dt_end.strftime('%Y%m%d %H:%M:%S')
    fp.write(tm1+"-"+tm2+"-->"+elaps(tm1,tm2)+"\n")

mf.mol.smiles=smi
pkl=fn.split(".")[0]+"_mf.pkl"
#with open(pkl,mode="wb") as ff:
#    dill.dump(mf,ff)

m=Chem.MolFromSmiles(smi)
mol=Chem.AddHs(m)

atom=toAtom(toList(Chem.MolToMolBlock(mol)))
mol3=gto.M(atom=atom,basis=basis,spin=spin)
emf=mol3.RKS()

emf.mol=mf.mol
emf.molBlock=toBlock(toArray(Chem.MolToMolBlock(mol)),mf.mol.atom_coords(unit="ANG"))
emf.mo_occ=mf.mo_occ
emf.mo_coeff=mf.mo_coeff
emf.mo_energy=mf.mo_energy

epkl=pkl.split('_mf.pkl')[0]+'_ext_mf.pkl'
with open(epkl,mode="wb") as ff:
    dill.dump(emf,ff)
