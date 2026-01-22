#!/bin/sh
#
#
fn=furanDerivative.smi
for star in furan-derivative-01 furan-derivative-02 furan-derivative-03 furan-derivative-04 furan-derivative-05 furan-derivative-06 furan-derivative-07 furan-derivative-08 furan-derivative-09 furan-derivative-10 furan-derivative-11 furan-derivative-12 furan-derivative-13;do
	x=$(cat $fn | grep $star )
echo "#!/bin/sh" > run.sh
echo "#" >> run.sh
echo "#" >> run.sh

echo python3 ../bin/tddftSolver.py $x -wo >> run.sh
sh run.sh
done
