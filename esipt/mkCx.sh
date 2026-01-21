#!/bin/sh
#
#
id=AAP
smiles="'c1(C(=O)C)ccccc1N'"
op=3,19,
i=1
echo "#!/bin/sh"
echo "#"
echo "#"
echo "com=~/uv/bin/esipt.py"
while [ "$i" -le 22 ]; do
	n=$(expr 210 - $(expr $i \* 5))
#	n=$(expr 95 + $(expr $i \* 5))
	up=$(expr $n / 100);dwn=$(expr $n % 100)
	if [ $dwn -lt 10 ]; then
	num=$op$up".0"$dwn
	else
	num=$op$up"."$dwn
	fi
	echo python3 \$com $id $smiles -phase 3 -distance $num -pt $i
	i=$(expr $i + 1)
done
