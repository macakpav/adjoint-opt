#!/bin/sh

subdir="optimisation/objective/0/angleas1"

gnuplot<<EOF
	set terminal png size 800,600
	set output 'figKonvergenceUy.png'
	set grid
	set title "Konvergence optimalizace U_{yTarget} - původní úhel 21.46°"
	set xlabel "Iterace [1]"
	set ylabel "Relativní rozdíl U_y a U_{yTarget} [1]"
    set yrange [0.00001:1]
	set logscale y
	plot "angleUy_45/$subdir" using 1:(abs((\$4-\$5)/\$5)) w lp t "- 4°", \
         "angleUy_50/$subdir" using 1:(abs((\$4-\$5)/\$5)) w lp t "- 2°", \
         "angleUy_62/$subdir" using 1:(abs((\$4-\$5)/\$5)) w lp t "+2°", \
         "angleUy_68/$subdir" using 1:(abs((\$4-\$5)/\$5)) w lp t "+4°", \
         "angleUy_74/$subdir" using 1:(abs((\$4-\$5)/\$5)) w lp t "+6°", \
         "angleUy_81/$subdir" using 1:(abs((\$4-\$5)/\$5)) w lp t "+8°", \
		0.01 with lines lt 2 notitle
EOF

subdir="postProcessing/outletAverage/50/surfaceFieldValue_50.dat"
a45=17.45819002769148
a50=19.458190027691480
a62=23.458190027691480
a68=25.4581900276915
a74=27.45819002769148
a81=29.45819002769148


gnuplot<<EOF
	set terminal png size 800,600
	set output 'figKonvergenceAlpha.png'
	set grid
	set title "Konvergence optimalizace U_{yTarget} - původní úhel 21.46°"
	set xlabel "Iterace [1]"
	set ylabel "Relativní rozdíl alpha_2 a alpha_{2Target} [1]"
    set yrange [0.00001:1]
	set logscale y
	plot "angleUy_45/$subdir" using 1:(abs((atan(\$3/\$2)*180/3.14-$a45)/$a45)) w lp t "- 4°", \
         "angleUy_50/$subdir" using 1:(abs((atan(\$3/\$2)*180/3.14-$a50)/$a50)) w lp t "- 2°", \
         "angleUy_62/$subdir" using 1:(abs((atan(\$3/\$2)*180/3.14-$a62)/$a62)) w lp t "+2°", \
         "angleUy_68/$subdir" using 1:(abs((atan(\$3/\$2)*180/3.14-$a68)/$a68)) w lp t "+4°", \
         "angleUy_74/$subdir" using 1:(abs((atan(\$3/\$2)*180/3.14-$a74)/$a74)) w lp t "+6°", \
         "angleUy_81/$subdir" using 1:(abs((atan(\$3/\$2)*180/3.14-$a81)/$a81)) w lp t "+8°", \
		0.01 with lines lt 2 notitle
EOF

# for d in ...
#     get Uy2
#     calc targetUy-Uy2
#     calc alpha2
#     get forceY