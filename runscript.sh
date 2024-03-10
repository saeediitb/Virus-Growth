#!/bin/bash
test_run()
{
  paramfolder=$1
  cd $paramfolder
  echo python ../run.py $paramfolder.json
  cd ../
}
run()
{
  paramfolder=$1
  cd $paramfolder
  python ../run.py $paramfolder.json
  cd ..
}
if [ -d $1 ]; then
	test_run $1
	echo "Enter to continue, Ctrl+c to end"
	read option
	if [[ -z $option ]]; then
	 	run $1
	else
		echo "run program denied!"
	fi 
elif [[ $1 = 'all' ]]; then
	for param in R*; do
		if [[ "$param" != "README.md" ]]; then
    		echo "$param"
    		test_run $param
		fi
	done
	echo "Enter to continue, Ctrl+c to end"
	read option
	if [[ -z $option ]]; then
	 	for param in R*; do
			if [[ "$param" != "README.md" ]]; then
				echo $param
				run $param
			fi
		done
	else
		echo "run program denied!"
	fi 
elif [[ $1 = 'allnew' ]]; then
	for param in R*; do
		cd $param
		if [[ $(ls -1 | wc -l) -eq 1 ]]; then
		    echo python ../run.py $param.json
		fi
		cd ..
	done
	echo "Enter to continue, Ctrl+c to end"
	read option
	if [[ -z $option ]]; then
	 	for param in R*; do
		    cd $param
			if [[ $(ls -1 | wc -l) -eq 1 ]]; then
			    python ../run.py $param.json
			fi
			cd ..
	    done
	else
		echo "run program denied!"
	fi 
else
    echo "no such param option"
    echo "give argv as param folder name or all"
fi