#WORK IN PROGRESS B-FIELD PLOTTING FILE


#sets gnuplot to output a pdf
#set term pdf enhanced
#self-explanatory
#set output 'heat_map.pdf'

set term png
set output 'heat_map.png'

#
#set view map

set view 60,30


#x and y axis labels
#set xlabel 'x [m]'
#set ylabel 'y [m]'

#somehow interpolates z values
set dgrid3d

#palette defined to be rainbow-ish
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)

#turns interpolation off because our datafile is large enough where it's not needed
set pm3d interpolate 0,0

splot "Bfield.dat" every ::1 using 1:2:3:4:5:6 with pm3d #gives pretty 2D color thing
#splot "Bfield-config1.dat" every ::1 using 1:3:4:6 with pm3d
#splot "Bfield-config1.dat" every ::1 using 1:2:3:4:5:6 with vectors #makes the fun crazy 3d plot
#splot "Bfield-config1.dat" u 1:2:3:(column(4)):(column(5)):(column(6)) with vectors #same as above
#splot "Bfield-config1.dat" u 1:3:(column(4)):(column(6)) with vectors #"all x vals undef"
#splot "Bfield-config1.dat" u 1:3:(column(4)):(column(6)) palette
# splot "Bfield-config1.dat" u 1:3:(sqrt((column(4))**2+(column(6))**2+(column(5)**2))) palette
#'every ::1' skips the first line

#Error we get withou dgrig3d
#Warning: Single isoline (scan) is not enough for a pm3d plot.
#         Hint: Missing blank lines in the data file? See 'help pm3d' and FAQ.
#       "plot.plt" line 25: all scans empty