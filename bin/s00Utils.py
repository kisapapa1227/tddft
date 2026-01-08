import glob,os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import  DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdFingerprintGenerator

from sklearn import model_selection, metrics
import numpy as np
from pyscf import gto
from pyscf.lib import param
import xgboost as xgb
import dill

def fix(x):
    if x<0.0:
        return f"   {x:.4f}"
    return f"    {x:.4f}"

def molBlockToAtom(MB):
    ret=[]
    MB=shapeIt(MB)
    for i in range(4,len(MB)):
        x=MB[i].split()
        if len(x)>10:
            ret.append([x[3],x[0],x[1],x[2]])
#        print(i,MB[i].split())
    return ret

def molBlockToFiles(smi,eq,fn):
    ff=open(fn,"w")

    mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
    m_h=Chem.AddHs(mol)

    MB=shapeIt(Chem.MolToMolBlock(m_h))
    t=0
    coords=eq.atom_coords()*param.BOHR
    for i in range(4):
        t+=1;ff.write(MB[i]+"\n")
    for i in range(eq.natm):
        t+=1;
        x,y,z=coords[i]
        ff.write(fix(x));ff.write(fix(y));ff.write(fix(z))
        ff.write(f" {eq.atom_symbol(i)}    0  0  0  0  0  0  0  0  0  0  0  0\n")
    for i in range(t,len(MB)):
        ff.write(MB[i]+"\n")
    ff.close()

def shapeIt(blk):
    ln=[];st=''
    for i,l in enumerate(blk):
        if ord(l)!=10:
            st+=l
        else:
            ln.append(st)
            st=''
    return ln

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

def toAtom(ss):
    atom="";
    for i in range(4,len(ss)):
        x=ss[i].split(',')
        if len(x)<10:
            return atom
        if atom=='':
            atom=x[4]+' '+x[1]+' '+x[2]+' '+x[3]
        else:
            atom=atom+';'+x[4]+' '+x[1]+' '+x[2]+' '+x[3]
    return atom

def showSpectrum(result_all,filename_list,prm):

    result_=prm['result_']
    for j,res in enumerate(result_all):
        for r in res:
            r2,y_test,y_pred,idx,_=r
            for y_, y_p, i in zip(y_test,y_pred,idx):
                result_[i].append((j,y_, y_p, r2))

    for i in range(len(result_)):
        rr=result_[i]
        plt.plot([prm['bins'][j]  for j, y_, y_p, _ in rr], [y_  for j, y_, y_p, _ in rr],label="true")
        plt.plot([prm['bins'][j]  for j, y_, y_p, _ in rr], [y_p for j, y_, y_p, _ in rr],label="predict")
  #plt.xticks(bins)
        plt.legend()
        plt.title(str(i)+":"+filename_list[i])
        name,_=os.path.splitext(os.path.basename(filename_list[i]))
        filename=prm['out_path']+"/fig_"+name+".png"
        plt.savefig(filename)
        plt.cla()

def preproc(mode,prm):

    y_scale=prm['scale']
    fin="uv_"+mode+"_dataset.tsv"
    df=pd.read_csv(fin,sep="\t",index_col=0)
    os.makedirs("img_"+mode,exist_ok=True)
    print(df)

    mol_vec={}
    for i,row in df.iterrows():
        mol=Chem.MolFromSmiles(row["SMILES"])
#        bitvec=AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
        fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        bitvec = fpgen.GetFingerprint(mol)
        v = np.zeros(2048)
        DataStructs.ConvertToNumpyArray(bitvec, v)
        mol_vec[row["filename"]]=v.tolist()

    df_chem_feat = pd.DataFrame(data=mol_vec.values(),index=mol_vec.keys(),
        columns=["fp"+str(i) for i in range(2048)])

    df=pd.merge(df,df_chem_feat,left_on="filename",right_index=True,how='right')
    print(df)

    attrs_x=[]
    attrs_y=[]
    attrs_title=[]

    for target in df.columns:
        if "trans" in target or "abs" in target:
            attrs_y.append(target)
        elif "fp" in target:
            attrs_x.append(target)

    attrs_new_x=[]
    for target in df.columns:
        if target in attrs_x:
            val=df[target].std()
            if val>0.2:
                attrs_new_x.append(target)

    print(len(attrs_new_x),"/",len(attrs_x))

    X_=df[attrs_new_x].values
    Y_=df[attrs_y].values/y_scale
    filename_list=df["filename"].values
    return X_,Y_,filename_list,attrs_new_x

def experiment(X,Y,title,save_filename):
    kfold=model_selection.KFold(n_splits=5, shuffle=True, random_state=1234)
    result_list=[]
    for train_index, test_index in kfold.split(X,Y):
        X_train = X[train_index]
        y_train = Y[train_index]
        X_test  = X[test_index ]
        y_test  = Y[test_index ]
        m=~np.isnan(y_train)
        y_train=y_train[m]
        X_train=X_train[m,:]
        model=xgb.XGBRegressor()
        model.fit(X_train, y_train)

        y_pred=model.predict(X_test)
        m=~np.isnan(y_test)
        if np.sum(m)>0:
            y_test_=y_test[m]
            y_pred_=y_pred[m]
            X_test_=X_test[m,:]
            r2  = metrics.r2_score(y_test_,y_pred_)
        else:
            r2 = np.nan
        result_list.append((r2,y_test,y_pred,test_index,model))

    plt.figure()
    plt.xlim([-0.2, 1.2])
    plt.ylim([-0.2, 1.2])
    for r2,y_test,y_pred,idx,_ in result_list:
        plt.scatter(y_test,y_pred,label='R2 %0.2f' % r2)
    plt.xlabel('True')
    plt.ylabel('Prediction')
    plt.title(title)
    plt.legend(loc="lower right")
    if save_filename:
        plt.savefig(save_filename)
        plt.cla()
#    plt.show()
    return result_list

def mkRecords(data_bins,smiles,key,segment):
    bins_head=["filename","SMILES"]+[key+str(l) for l in range(segment[0],segment[1],segment[2])]
    data={item:[] for item in bins_head}
    for k,bins in data_bins.items():
        if k in smiles:
            print(k,sep='')
            data["filename"].append(k)
            print(smiles[k],sep='')
            data["SMILES"].append(smiles[k])
            print(data_bins[k])
            bins=data_bins[k]
            for i, b in enumerate(bins):
                data[bins_head[i+2]].append(np.mean(b))
    return data

def getPeaks(xs,ys,exel=False,xlim=[250,450]):
    r=[];m=0
    for i in range(1,len(ys)-1):
        if ys[i-1]<ys[i]:
            if ys[i]>ys[i+1]:
                if xs[i]<xlim[1] and xs[i]>xlim[0]:
                    r.append(xs[i])
    if exel:
        return r[::-1]
    s=""
    r=r[::-1]
    for x in r:
        if s!="":
            s=s+","
        s=s+str("{:.1f}".format(x))
    return s

def getCP(x,w,f,s):
    composite=0
    for count, peak in enumerate(w):
        this_peak=gaussBand(x,peak,f[count],s)
        composite+=this_peak
    return composite

def gaussBand(x, band, strength, stdev):
  constant=0.0003
  bandshape = constant * (strength/(1.0/stdev))*np.exp(-(((1.0/x)-(1.0/band))**2/(1.0/stdev)**2))
  return bandshape

def getNP(fn):
    with open(fn+".tdd","r") as fp:
        t1=fp.readline().split(",")# 1
        p=[]
        if len(t1)<2:
            print(fn,":No converged data:",len(t1),t1)
            return [],0,0
        for n in range(len(t1)-1):
            p.append(float(t1[n]))
        wavlen=np.array(p)

        t1=fp.readline().split(",")# 2
        p=[]
        for n in range(len(t1)-1):
            p.append(float(t1[n]))
        line=fp.readline().split("\n")
#        print(fn,"-->",line)
        smiles=line[0]
        h=False;l=False
        for line in fp:
            if 'HOMO' in line:
                h=float("".join(line.split(':')[1].split(' ')[:-1]))
            if 'LUMO' in line:
                l=float("".join(line.split(':')[1].split(' ')[:-1]))
    return wavlen,np.array(p),smiles,[h,l]

def getSmiles(filename,key):
    filename="data_info.xlsx"
    df=pd.read_excel(filename)

    smiles={}
    for k,v in zip(df["製品名"],df["SMILES"]):
        if k is not None and type(k) is str:
    ####
            path=key+"/"+k+".tsv"
            f=os.path.exists(path)
            if f:
                print(">>",path,v)
                smiles[path]=v
            else:
                print(path)
    return smiles

def standarlize(data,segment):
    data_bins={}
    for k,v in data.items():
        xs,ys=v
        bins=[[] for l in range(segment[0],segment[1],segment[2])]
        for x,y in zip(xs,ys):
            i=int((x-segment[0]))//segment[2]
            if i<len(bins):
                bins[i].append(y)
        data_bins[k]=bins
    return data_bins

def readTsv(dataPath,key):
    for filename in glob.glob(dataPath+"/*.xlsx"):
        path,_=os.path.splitext(filename)
        out_filename=path+".tsv"
        if not os.path.exists(out_filename):
            df=pd.read_excel(filename)
            print("converting ->",out_filename)
            df.to_csv(out_filename,sep="\t",index=False)
    data={}
    for filename in glob.glob(dataPath+"/*.tsv"):
        path,_=os.path.splitext(filename)
        df=pd.read_table(filename)
  ####
        x,y=[],[]
        for k,v in zip(df["波長[nm]"],df[key]):
            x.append(k)
            y.append(v)
        data[filename]=(x,y)
#    ret=sorted(data.items(),key=lambda x:x[0])
    ret={}
    for fn in sorted(data.keys()):
        ret[fn]=data[fn]
    return ret

def saveAs(img,n=None,fn=False):
    if not fn:
        if n==None:
            fn=sys.argv[0].split(".py")[0]+".png"
        else:
            fn=sys.argv[0].split(".py")[0]+n+".png"

    if not hasattr(img,"__name__"):
        if type(img)==str:# assume as SVG
            saveAsSVG(img,fn)
            return
        # assume as IPython.IP
        print("Asave as "+fn)
        img.save(fn)
        return

    print("save as "+fn)
    if img.__name__=="matplotlib.pyplot":
        img.savefig(fn)
        return
    img.save(fn)

def saveAsSVG(img,fn):
    ff=fn.split(".png")[0]+".svg"
    with open(fn,mode='w') as f:
        f.write(img)
    print("save as "+fn)

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

def toArray(mb):
    ss=[];s="";f=False
    for i in mb:
        if ord(i)==10:
            ss.append(s);s=""
        else:
            s=s+i
    return ss

def writeEnergy(mol_block,mf,argv,d=False):

    m0=Chem.MolFromSmiles(argv[2])
    sss=Chem.MolToSmiles(m0)
    m=Chem.MolFromSmiles(sss)
    mol=Chem.AddHs(m)

    mol4=gto.Mole(atom=mol_block).build() #1
    emf=mol4.RKS()

    emf.mol=mf.mol
    emf.molBlock=toBlock(toArray(Chem.MolToMolBlock(mol)),mf.mol.atom_coords(unit="ANG"))
    emf.mo_occ=mf.mo_occ
    emf.mo_coeff=mf.mo_coeff
    emf.mo_energy=mf.mo_energy

    epkl=argv[1]+'_ext_mf.pkl'
    with open(epkl,mode="wb") as ff:
        dill.dump(emf,ff)

    if d!=False:
        with open(argv[1]+".erg","a") as ff:
            ff.write(f'{d} {mf.e_tot}\n')
