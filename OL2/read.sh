#!/bin/bash
Dirname=bash_homework
CUR_DIR=`ls $Dirname`

for i in $CUR_DIR
do
	if [ -f $Dirname/$i ];then
		echo $i >> "/home/test/linux/filename.txt"
	elif [ -d $Dirname/$i ];then
		echo $i >> "/home/test/linux/dirname.txt"
	else
		echo "Cannot recognize!"
	fi
done

exit 0