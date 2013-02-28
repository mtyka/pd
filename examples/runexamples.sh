
basedir=$(pwd)
targetdir="log"
rm -rf targetdir
mkdir targetdir
export PD_PARAM_PATH="$basedir/../param"
for j in $(find -name "run.py"); do
	i=$(echo $j | sed 's/\.\///')
	ename=$(echo $i | sed "s/run\.py//" )
	sddecname=$(echo $i | sed "s/\//_/g" | sed "s/\.py//" )
	cp ../bin/lib/python2.5/site-packages/pd/*  $ename
	cd $ename
	echo --------- Running $i in $ename -------------------------------------
	python run.py > $basedir/$targetdir/$sddecname.log
	cd $basedir

done

#sh clean.sh
