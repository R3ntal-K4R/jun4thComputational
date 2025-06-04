gnuplot -pe << EOF
set term pdf enhanced
set output 'heat_map.pdf'
set view map
set dgrid3d
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
set pm3d interpolate 0,0
splot "bfield_test_file.dat" using 1:2:3 with pm3d
EOF
