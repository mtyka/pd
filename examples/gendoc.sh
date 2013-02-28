
echo /\*!
echo @page examplelist List of Example Scripts
echo @section mainpage_examples Example Scripts
echo 
echo The example scripts can be found in examples/
echo
for j in $(find -name "run.py"); do
	i=$(echo $j | sed 's/\.\///')
	ename=$(echo $i | sed "s/run\.py//" )
	secname=$(echo $i | sed "s/\//_/g" | sed "s/\.py//" )

	## now do the Example list
	echo "  - \ref $secname" 

done



echo
echo
echo
echo
echo
echo
for j in $(find -name "run.py"); do
	i=$(echo $j | sed 's/\.\///')
	ename=$(echo $i | sed "s/run\.py//" )
	secname=$(echo $i | sed "s/\//_/g" | sed "s/\.py//" )
	echo "\page  $secname Example: example/"$ename"run.py"  

	grep "^###" $j  | sed 's/###//'

	echo "\include \""$ename"run.py\""  
	echo 
done

echo \*/
echo
