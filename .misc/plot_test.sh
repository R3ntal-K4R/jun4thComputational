math -pe
data = Import[/home/gavin/Missfit_Code/drakemissfit/"output.dat", "CSV"]
Plot3D[data[[All, {2,3,4} ]], {x,0,40}, {y,0,40}, {z,0,40}]


gnuplot -pe << EOF
set term pdf enhanced
set xlabel "x"
set ylabel "y"
set zlabel "z"
set output 'threedplot.pdf'
set title 'Particle Traveling in 3d' 
splot 'output.dat' using 2:3:4 with lines not
EOF
