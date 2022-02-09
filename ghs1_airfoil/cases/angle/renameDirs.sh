#!/bin/bash

uy=(45 50 62 68 47 81)
ang=(17 19 23 25 27 29)

for (( i=0; i<6; i++ )); do
    echo angleUy_${uy[$i]} angle${ang[$i]}_Uy
    #mv angleUy_${uy[$i]} angle${ang[$i]}_Uy
done