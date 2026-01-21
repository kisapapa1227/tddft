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
initA=AAP
initB=uv_MTG03042025_esipt_run2_plot_all
initC=xAAP
i=1
echo "#!/bin/sh" > test.sh
echo "#" >> test.sh
echo "#" >> test.sh
echo -n "convert -loop 0 " >> test.sh

while [ "$i" -le 22 ] ; do
setNum $i
p1=${initB}_${i}.png
p2a=${initA}_${p}_proj.svg
p2=${initA}_${p}_proj.png
convert ${p2a} -resize 40% ${p2}
p3=${initC}_${p}.png

convert $p1 $p2 -gravity northeast -geometry +80+100 -compose over -composite $p3

	echo -n " " ${p3} >> test.sh
	i=$(expr $i + 1)
done
	echo " " ${initC}.gif >> test.sh

	sh test.sh
