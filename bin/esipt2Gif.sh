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


ex=$(basename $PWD)
init=$(python3 ../bin/getInit.py $PWD)
prm=${init}${ex}.prm
initB=$(python3 ../bin/getPptx.py $prm)
gif=${init}_${ex}_cb.gif
i=1
echo "#!/bin/sh" > test.sh
echo "#" >> test.sh
echo "#" >> test.sh
echo -n "convert -loop 0 " >> test.sh

while [ "$i" -le 30 ] ; do
setNum $i
p1=${initB}_${i}.png
p2a=${init}_${p}_proj.svg
p2=${init}_${p}_proj.png
p3=${init}_${p}.png

if [ -f ${init}_${p}.tdd ]; then
	convert ${p2a} -resize 40% ${p2}
	convert $p1 $p2 -gravity northeast -geometry +80+100 -compose over -composite $p3
	echo -n " " ${p3} >> test.sh
echo $p1 $p2a $p2 $p3
fi
	i=$(expr $i + 1)

#echo $p1 $p2a $p2 $p3

done
	echo " " ${gif} >> test.sh

	sh test.sh

echo ${gif}
