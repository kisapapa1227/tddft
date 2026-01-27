import sys
import math

def trans(pt,sx,scale,mg):
    pp=[]
    for ok in pt:
        x=ok[0];y=ok[1]
        pp.append([(x-sx[0][0])*scale+mg,(y-sx[0][1])*scale+mg])
    return pp

def fconv(s):
    ret=[];p=0
    for i in s:
        if i!='':
            if p!=3:
                ret.append(float(i))
                p=p+1
            else:
                ret.insert(0,i)
                return ret
    return ret

def iconv(s):
    ret=[]
    for i in s:
        if i!='':
            ret.append(int(i))
    return ret

def setPrm(ls):
    flip=False;rotate=0.0;dist=False
    for l in ls:
        i=l.split('\n')[0].split(':')
        match i[0]:
            case 'flip':
                if i[1]=='True':
                    flip=True
            case 'rotate':
                rotate=float(i[1])
            case 'distance':
                if not '-' in i[1]:
                    dist=float(i[1])

    return flip,rotate,dist

def letRotate(sx,pps,rotate):
    ox,oy=(sx[0][0]+sx[1][0])/2.0,(sx[0][1]+sx[1][1])/2.0
    pp=[]
    c=math.cos(rotate/180.0*3.14159265)
    s=math.sin(rotate/180.0*3.14159265)
    for p in pps:
        px,py=p[0]-ox,p[1]-oy
        pp.append([c*px-s*py+ox,s*px+c*py+oy])
    return pp

def letFlip(sx,pps):
    pp=[]
    ss=sx[1][0]+sx[0][0]
    for p in pps:
        pp.append([ss-p[0],p[1]])
    return pp

def setMaxMin(pt):
    x0=y0=100.0;x1=y1=-100.0
    for ok in pt:
        x=ok[0];y=ok[1]
        if x0>x:
            x0=x
        if y0>y:
            y0=y
        if x1<x:
            x1=x
        if y1<y:
            y1=y
    return [[x0,y0],[x1,y1]]

def setDist(dat,i0,i1):
    rx=dat[i0-1][1]-dat[i1-1][1]
    ry=dat[i0-1][2]-dat[i1-1][2]
    rz=dat[i0-1][3]-dat[i1-1][3]
    dst=math.sqrt(rx*rx+ry*ry+rz*rz)
    return dst

flip=False
rotate=0
dist=False
fn=sys.argv[1]
with open(fn,"r") as fp:
    fn=fp.readline().split('\n')[0]+".txt"
    nt=fp.readline().split('\n')[0].split(',')
    flip,rotate,dist=setPrm([fp.readline(),fp.readline(),fp.readline()])

pt=[]
for i in nt:
    pt.append(int(i)-1)

dat=[];cnt=[]
step=0
with open(fn,"r") as fp:
    for l in fp:
        if 'END' in l:
            break
        p=l.split('\n')[0].split(' ')
        match step:
            case 0:
                if len(p)>30:
                    step=1
            case 1:
                if len(p)<30:
                    step=2
        match step:
            case 1:
                dat.append(fconv(p))
            case 2:
                cnt.append(iconv(p))


va=[dat[pt[1]][1]-dat[pt[0]][1],\
    dat[pt[1]][2]-dat[pt[0]][2],\
    dat[pt[1]][3]-dat[pt[0]][3]]

if dist!=False:
    l=math.sqrt(va[0]*va[0]+va[1]*va[1]+va[2]*va[2])
    for i in range(3):
        dat[pt[1]][i+1]=dat[pt[0]][i+1]+va[i]/l*dist

vb=[dat[pt[2]][1]-dat[pt[0]][1],\
    dat[pt[2]][2]-dat[pt[0]][2],\
    dat[pt[2]][3]-dat[pt[0]][3]]
vn=[va[1]*vb[2]-va[2]*vb[1],va[2]*vb[0]-va[0]*vb[2],va[0]*vb[1]-va[1]*vb[0]]

#print("a.n",va[0]*vn[0]+va[1]*vn[1]+va[2]*vn[2])
#print("b.n",vb[0]*vn[0]+vb[1]*vn[1]+vb[2]*vn[2])

bt=math.atan2(-vn[0],vn[1])
#print("mx",vn[0]*math.cos(bt)+vn[1]*math.sin(bt))
my=-vn[0]*math.sin(bt)+vn[1]*math.cos(bt)
al=math.atan2(-my,vn[2])
#print("nx",my*math.cos(al)+vn[2]*math.sin(al))

sa=math.sin(al);ca=math.cos(al)
sb=math.sin(bt);cb=math.cos(bt)

i=pt[0]
hx= dat[i][1]*cb+dat[i][2]*sb
hy=-dat[i][1]*sb+dat[i][2]*cb
hz= dat[i][3]

gx= hx
gy= hy*ca+hz*sa
gz=-hy*sa+hz*ca

ss=gz

pp=[]
for i in range(len(dat)):
    hx= dat[i][1]*cb+dat[i][2]*sb
    hy=-dat[i][1]*sb+dat[i][2]*cb
    hz= dat[i][3]

    gx= hx
    gy= hy*ca+hz*sa
    gz=-hy*sa+hz*ca-ss
    pp.append([gx,gy])
    #print(dat[i],"-->",f'{gx:.4f}.,{gy:.4f},{gz:.4f}')

mg=0.1
r=0.2
mg=30.0
size=500.0

sx=setMaxMin(pp)

if flip:
    pp=letFlip(sx,pp)

if rotate!=0.0:
    pp=letRotate(sx,pp,rotate)
    sx=setMaxMin(pp)

if (sx[1][0]-sx[0][0])>(sx[1][1]-sx[0][1]):
    scale=(size-mg*2)/(sx[1][0]-sx[0][0])
else:
    scale=(size-mg*2)/(sx[1][1]-sx[0][1])

pp=trans(pp,sx,scale,mg)

#print("mg",mg)
mx=int(scale*(sx[1][0]-sx[0][0])+mg*2)
my=int(scale*(sx[1][1]-sx[0][1])+mg*2)

fp=open(fn.split(".")[0]+"_proj.svg","w")
fp.write(f'<svg width="{mx}" height="{my}" viewBox="0, 0, {mx}, {my}">\n')
r=r*scale

p=0.3
for i in cnt:
    x1=pp[i[0]-1][0];y1=pp[i[0]-1][1]
    x2=pp[i[1]-1][0];y2=pp[i[1]-1][1]
    al=math.atan2(y2-y1,x2-x1)
    xx1=x1+r*math.cos(al);yy1=y1+r*math.sin(al);
    xx2=x2-r*math.cos(al);yy2=y2-r*math.sin(al);

    dst=setDist(dat,i[0],i[1])

    match i[2]:
        case 1:d="-"
        case 2:d="="
        case 3:d="#"
    print(str(i[0]).zfill(2)+"-"+str(i[1]).zfill(2)+" : "+format(dst,'.2f')+"Å "+dat[i[0]-1][0]+d+dat[i[1]-1][0])
    match i[2]:
        case 1:
            width=scale*0.05
            fp.write(f'<line x1="{xx1}" y1="{yy1}" x2="{xx2}" y2="{yy2}" stroke="#000000" stroke-width="{width}"/>\n')
        case 2:
            width=scale*0.1
            fp.write(f'<line x1="{xx1}" y1="{yy1}" x2="{xx2}" y2="{yy2}" stroke="#000000" stroke-width="{width}"/>\n')
            width=scale*0.03
            fp.write(f'<line x1="{xx1}" y1="{yy1}" x2="{xx2}" y2="{yy2}" stroke="#ffffff" stroke-width="{width}"/>\n')
        case 3:
            width=scale*0.15
            fp.write(f'<line x1="{xx1}" y1="{yy1}" x2="{xx2}" y2="{yy2}" stroke="#000000" stroke-width="{width}"/>\n')
            width=scale*0.10
            fp.write(f'<line x1="{xx1}" y1="{yy1}" x2="{xx2}" y2="{yy2}" stroke="#ffffff" stroke-width="{width}"/>\n')
            width=scale*0.03
            fp.write(f'<line x1="{xx1}" y1="{yy1}" x2="{xx2}" y2="{yy2}" stroke="#000000" stroke-width="{width}"/>\n')
for i in range(len(dat)):
    cx=pp[i][0];cy=pp[i][1]
    fp.write(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="#7f7f7f"/>\n')
    cx=pp[i][0]-5;cy=pp[i][1]+5;st=dat[i][0]
    fp.write(f'<text x="{cx}" y="{cy}" font-family="serif" font-weiht="bold" font-size="10" fill="#ffffff">\n')
    fp.write(f'{st}\n')
    fp.write("</text>\n")
    cx=cx+10
    fp.write(f'<text x="{cx}" y="{cy}" font-family="serif" font-weiht="bold" font-size="10" fill="#000000">\n')
    fp.write(f'#{i+1}\n')
    fp.write("</text>\n")
fp.write('</svg>\n')
fp.close()

dst=setDist(dat,pt[0]+1,pt[1]+1)
print("\n"+str(pt[0]+1).zfill(2)+"-"+str(pt[1]+1).zfill(2)+" : "+format(dst,'.2f')+"Å "+dat[pt[0]][0]+" "+dat[pt[1]][0])
