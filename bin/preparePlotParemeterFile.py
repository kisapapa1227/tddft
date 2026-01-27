import sys, os
import s00Utils
import numpy as np
import math

def getRange(w,f,std):
    x=[]
    for i in range(0,13):
        x.append(i*100+1)
    xx=np.array(x)
    cp=s00Utils.getCP(xx,w,f,std)

    mx=max(cp)
#    print(mx,len(cp))

    for low in range(1,12):
        if cp[low]>mx*0.1:
            break
    low=low-1
    for up in range(11,0,-1):
        if cp[up]>mx*0.1:
            break
    up=up+1
#    print("cp------------------------->",cp)
#    print("low",low,cp[low],cp[low+1])
#    print("up",up,cp[up],cp[up+1])

    a=int(math.log(mx))-1
    xx=int((mx*1.4)/math.pow(10,a))*math.pow(10,a)
    return low*100,up*100,xx
#
# main
#
key=False;std=4000
if len(sys.argv)>1:
    for i in range(1,len(sys.argv)):
        if sys.argv[i]=='-key':
            key=sys.argv[1]
            i=i+1
        elif sys.argv[i]=='-std':
            std=float(sys.argv[1])
            i=i+1

filename=[f.name for f in os.scandir() if f.is_file()]

print("files:",end="")
init=True
for fn in filename:
    if ".tdd" in fn and not ".tdd.log" in fn:
        if key!=False and not key in fn:
            continue
        if init!=True:
            print(",",end="")
        init=False
        fx=fn.split(".tdd")[0]
        w,f,s,h=s00Utils.getNP(fx)
        r=getRange(w,f,std)
#        f=getMax("f",f)

        print(fx,end="")
if key==False:
    print("\noutput:all")
else:
    print("\noutput:plot_"+key)
#print(sys.argv[0],sys.argv[1])
print("ops:["+str(r[0])+","+str(r[1])+","+str((r[1]-r[0])+1)+"],"+'{:.2f}'.format(r[2])+","+str(std))

