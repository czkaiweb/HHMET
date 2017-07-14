#! /bin/sh

cmd='delphes_analysis.C("'$1'")'
root -l -q -b $cmd
