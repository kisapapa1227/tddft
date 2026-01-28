import sys,os

key='tddft-git'

with open(sys.argv[1],"r") as fp:
    fp.readline()
    fni=fp.readline().split(':')[1].split("\n")[0]

wk=os.getcwd().split(key)[1].split('/')

fni='uv'+"_".join(wk)+"_"+fni

print(fni)
