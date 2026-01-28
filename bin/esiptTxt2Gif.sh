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

if [ $# -lt 1 ]; then
	molp=molpp.prm
else
	molp=$1
fi

init=$(python3 ../bin/getInit.py $PWD)

i=1
#echo "init" $init

ex=$(basename $PWD)
prm=${init}${ex}.prm
gif=${init}_${ex}_struct
echo $ex

sh ../bin/mkPrmForTdds2svg.sh  > $prm
python3 ../bin/tdds2svg.py $prm -num 1

echo "#!/bin/sh" > test.sh
echo "#" >> test.sh
echo "#" >> test.sh
echo -n "convert -loop 0 " >> test.sh

while [ "$i" -le 22 ] ; do
	setNum $i
	fn=${init}_${p}
	echo "---->" $fn
	echo $fn > test.prm
	tail -3 $molp >> test.prm
#	cat test.prm
	python3 ../bin/projectMol.py test.prm
	cat test.prm
	fn=${init}_${p}_proj.svg
	if [ ! -e $fn ]; then
		break
	fi
	echo -n " " $fn >> test.sh
	i=$(expr $i + 1)
done
	echo " " ${gif}.gif >> test.sh

	sh test.sh

sh ../bin/esipt2Gif.sh
