import sys
import matplotlib.pyplot as plt
import s00Utils

def getErg(fn):
    with open(fn+".erg") as ff:
        xx=[];yy=[]
        for l in ff:
            s=l.split()
            xx.append(float(s[0]));yy.append(float(s[1]))
    mx=0.0
    if yy[0]<0.0:
        for y in yy:
            if mx>y:
                mx=y
    else:
        for y in yy:
            if mx<y:
                mx=y

    hv = 27.2114;hl = 627.51
    if mx<0.0:
        c=hl
    else:
        c=-hv*hl

    rx=[];ry=[]
    for x,y in zip(xx,yy):
        rx.append(x)
        ry.append((y-mx)*c)
    return rx,ry
###############################################
# main
###############################################
xaxis="Bond Length(Å)"
with open(sys.argv[1]) as ff:
    for l in ff:
        if 'files' in l:
            fns=l.split('\n')[0].split(':')[1].split(',')
        if 'label' in l:
            lbs=l.split('\n')[0].split(':')[1].split(',')
        if 'output' in l:
            png=l.split('\n')[0].split(':')[1]
        if 'axis' in l:
            xaxis=l.split('\n')[0].split(':')[1]
        if 'ops' in l:
            ops=l.split(':')[1].split("\n")[0]
            xs=ops.split("]")[0].split("[")[1].split(',')
            ym=float(ops.split("]")[1].split(',')[1])

#print(fns);print(lbs);print(png);print(xs,ym)

Figure, ax=plt.subplots()
for fn,lb in zip(fns,lbs):
    x,y=getErg(fn)
    for i,j in zip(x,y):
        print(format(i,".2f"),format(j,"6.3f"))
    ax.plot(x,y,label=lb)

ax.set_ylim([0,ym])
ax.set_xlim([float(xs[0]),float(xs[1])])
ax.legend()
plt.xlabel(xaxis)
plt.ylabel("Energy (kcal/mol)")

s00Utils.saveAs(plt,fn=png+".png")

