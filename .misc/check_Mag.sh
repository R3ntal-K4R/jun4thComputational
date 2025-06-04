
awk '$2==110 && $3==110' Bfield_n300.dat > storage_x.dat
awk '$1==110&& $3==110' Bfield_n300.dat > storage_y.dat
awk '$1==110&& $2==110' Bfield_n300.dat > storage_z.dat






gnuplot -persist <<-EOFMarker



set title "Magnetic field along x direction @ y,z=150"
set xlabel "0.1 m in X direction"
set ylabel " Bfield"
plot "storage_x.dat" using 1:4 with lines title "Bx", "storage_x.dat" using 1:5 with lines title "By", "storage_x.dat" using 1:6 with lines title "Bz"
EOFMarker

gnuplot -persist <<-EOFMarker

set title "Magnetic field along y direction @ x,z=150"
set xlabel "0.1 m in Y direction"
set ylabel " Bfield"
plot "storage_y.dat" using 2:4 with lines title "Bx", "storage_y.dat" using 2:5 with lines title "By", "storage_y.dat" using 2:6 with lines title "Bz"
EOFMarker

gnuplot -persist <<-EOFMarker

set title "Magnetic field along z direction @ x,y=150"
set xlabel "0.1 m in Z direction"
set ylabel " Bfield"
plot "storage_z.dat" using 3:4 with lines title "Bx", "storage_z.dat" using 3:5 with lines title "By", "storage_z.dat" using 3:6 with lines title "Bz"

EOFMarker
