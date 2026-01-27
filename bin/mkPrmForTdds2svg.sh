#!/bin/sh
#
#
key=''

if [ $# -gt 0 ] ; then
	key=$1
fi
#	echo $key
echo -n "files:"

ln=$(ls ${key}*.tdd | wc | awk '{print $1}')
cn=1
for fn in $(ls ${key}*.tdd); do
	echo -n $(basename $fn ".tdd")
	if [ $cn != $ln ];then
		echo -n ','
	fi
	cn=$(expr $cn + 1)
done
if [ $# -eq 0 ] ; then
	key=all
fi
echo
echo "output:plot_"${key}
#echo "smiles:"${key}.prm
echo "ops:[100,600,501],1.2,4000"
