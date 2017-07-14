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
	    cp run.con run_$ystart'_'$fstart.con
#	    echo "Arguments      = $index" >> run_$index.con 
	    echo "Executable     = run_script/run_${1}_MY_${ystart}_MF_${fstart}.sh">>run_$ystart'_'$fstart.con
	    echo "Output         = log/job_${ystart}_${fstart}.out" >> run_$ystart'_'$fstart.con
	    echo "Error          = log/job_${ystart}_${fstart}.err" >> run_$ystart'_'$fstart.con
	    echo "Log            = log/job_${ystart}_${fstart}.log" >> run_$ystart'_'$fstart.con
	    echo "Queue"         >> run_$ystart'_'$fstart.con
    
	    condor_submit run_$ystart'_'$fstart.con 
              
	    mv run_$ystart'_'$fstart.con script/     
	    let "fstart+=50";
    	done
	let "ystart+=100";
	fstart="200"
done
