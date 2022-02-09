#!/bin/bash

for d in angle*; do
    cd $d
    name=$d
    outfile=../dataKO/$name.csv
    rm -f $outfile
    touch $outfile
    echo \# Iter Ux Uy Fx Fy >> $outfile

    for (( i=2; i<=50; i++ )); do
        dir=${name}_$i
        echo $i `tail -n 1 $dir/postProcessing/outletAverage/0/surfaceFieldValue.dat | awk '{$1=""; print $0}'` \
        `tail -n 1 $dir/postProcessing/force/0/force.dat | awk '{$1=$4=$5=$6=$7=$8=$9=$10=""; print $0}' | sed 's/(\(-\|[0-9]\)/\1/g'` >> $outfile
    done
    sed -i 's/  */ /g' $outfile
    sed -i 's/\([0-9]\) \([0-9]\|-\)/\1, \2/g' $outfile

    cd ..
done