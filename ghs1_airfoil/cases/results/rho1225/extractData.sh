#!/bin/sh

dirname=${PWD##*/}
outdir="data_$dirname"

[ ! -d $outdir ] && mkdir $outdir

for d in angle*; do
    join $d/postProcessing/outletAverage/50/surfaceFieldValue.dat $d/postProcessing/force/50/force.dat | \
    join - $d/optimisation/objective/0/*eas1 | \
    sed 's/Time//g' | sed 's/(total_x/total_x/g' | sed 's/total_z)/total_z/g' | \
    sed 's/(\(-\|[0-9]\)/\1/g' | \
    awk '{$6=$7=$8=$9=$10=$11=$12=$13=$14=$16=""; print $0}' | tail -n +24 | \
    sed 's/\s\+/, /g' | sed 's/, $//g' | sed 's/#,/# Iter,/' > $outdir/${d##*/}.csv
done

tar -cf $outdir.tar $outdir 