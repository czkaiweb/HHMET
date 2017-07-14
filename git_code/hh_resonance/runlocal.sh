#!/bin/bash

date
mode=$1
ystart="600"
yend="1300"
fstart="200"
fend="650"
chmod 755 ./run_script/*.sh

while [ $ystart -le $yend ]; do
	echo $ystart
	while [ $fstart -le $fend ]; do
	    echo $ystart $fstart
    
	    ./run_script/run_${1}_MY_${ystart}_MF_${fstart}.sh 
              
	    let "fstart+=50";
    	done
	let "ystart+=100";
	fstart="200"
done
