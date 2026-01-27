import sys

fn=sys.argv[1]

com="tddftSolver.py"
for i in range(2,len(sys.argv)):
    if sys.argv[i]=='-c':
        com=sys.argv[i+1]

print(f"#!/bin/sh\n#\n#\nfn={fn}\nfor star in",end="")

with open(fn) as fp:
    for l in fp:
        star=l.split(" ")[0]
        print(f" {star}",end="")
print(';do\n\tx=$(cat $fn | grep $star )\necho "#!/bin/sh" > run.sh')
print('echo "#" >> run.sh\necho "#" >> run.sh\n')
print('echo python3 ../bin/'+com+' $x >> run.sh\nsh run.sh\ndone')
