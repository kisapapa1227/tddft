#!/bin/sh
#
#

setNum(){
	if [ $1 -lt 10 ]; then
		p="0"$1
	else
		p=$1
	fi
	return $p
}
init=AAP
molp=molpp.prm
i=1

ex=$(basename $PWD)
prm=${init}${ex}.prm

sh ~/uv/bin/mkPrmForTdds2svg.sh  > $prm
python3 ~/uv/bin/tdds2svg.py $prm -num 1

echo "#!/bin/sh" > test.sh
echo "#" >> test.sh
echo "#" >> test.sh
echo -n "convert -loop 0 " >> test.sh

while [ "$i" -le 22 ] ; do
	setNum $i
	echo "yes" $p
	echo ${init}_${p} > test.prm
	tail -3 $molp >> test.prm
#	cat test.prm
	python3 ~/uv/bin/projectMol.py test.prm
	cat test.prm
	echo -n " " ${init}_${p}_proj.svg >> test.sh
	i=$(expr $i + 1)
done
	echo " " ${init}.gif >> test.sh

	sh test.sh
