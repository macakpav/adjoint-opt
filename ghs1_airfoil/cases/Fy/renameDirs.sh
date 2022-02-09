#!/bin/bash

fy=(9.1 8.6 7.6 7.1 6.6 6.1)
ang=(17 19 23 25 27 29)

for (( i=0; i<6; i++ )); do
    echo angleFy_${fy[$i]} angle${ang[$i]}_Fy
    #mv angleFy_${fy[$i]} angle${ang[$i]}_Fy
done