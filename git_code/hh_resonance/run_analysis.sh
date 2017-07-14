#! /bin/sh

cmd='delphes_analysis.C("'$1'")'
cmd2='delphes_boosted.C("'$1'")'
root -l -q -b $cmd 
root -l -q -b $cmd2
