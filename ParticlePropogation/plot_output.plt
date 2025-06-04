set terminal pdfcairo enhanced color size 12in, 8in
set xzeroaxis linetype 3 linewidth 2.5

if (ARGC < 2) {
    print "Usage: gnuplot -c plot_output.gp output.dat bfield_output.dat";
    exit;
}


#Position plots
set output 'plot_particle.pdf'
input_file = ARG1
bfield = ARG2

set multiplot layout 2, 3
# unset xtics
# unset ytics

#tic stuff to make prettier
set xtics font ",10"
set xtics rotate by -45
set xtics nomirror out
set ytics font ",10"
set ytics nomirror out
num_xtics = 4
num_ytics = 5
stats input_file u 1 name 'x' nooutput
set xtics x_max/num_xtics #this will be the same for all graphs using "output.dat"
set format x "%6.4g" # scientific notation 
set format y "%6.4g"
# set size square 0.5,0.5

# 1
set title "x vs time"
set xlabel "timestep"
set ylabel "x [m]"
# stats input_file u 3 name 'b' nooutput
# set ytics b_max/2
plot input_file u 1:3 with lines notitle

# 2
set title "y vs time"
set xlabel "timestep"
set ylabel "y [m]"
# stats input_file u 4 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:4 with lines notitle

# 3
set title "z vs time"
set xlabel "timestep"
set ylabel "z [m]"
# stats input_file u 5 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:5 with lines notitle

# 4
set title "vx vs time"
set xlabel "timestep"
set ylabel "vx [m/s]"
# stats input_file u 6 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:6 with lines notitle

# 5
set title "vy vs time"
set xlabel "timestep"
set ylabel "vy [m/s]"
# stats input_file u 7 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:7 with lines notitle

# 6
set title "vz vs time"
set xlabel "timestep"
set ylabel "vz [m/s]"
# stats input_file u 8 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:8 with lines notitle


unset multiplot

#energy loss time
set multiplot layout 2, 3

set title "Velocity magnitude vs time"
set ylabel "[m/s]"
plot input_file u 1:12 with lines notitle

set title ""

set title "Total Energy loss [eV] vs time"
# unset title
unset xlabel
unset ylabel
set xlabel "timestep"
# set ylabel "energy loss"
# set size 0.45, 0.45
# stats input_file u 13 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:13 with lines notitle

set title "Electronic Energy loss [eV] vs time"
set xlabel "timestep"
# set ylabel "electronic energy loss"
# set size 0.45, 0.45
# stats input_file u 14 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:14 with lines notitle

set title "Nuclear Energy loss [eV] vs time"
set xlabel "timestep"
# set ylabel "energy loss"
# set size 0.45, 0.45
# stats input_file u 15 name 'b' nooutput
# set ytics b_max/num_ytics
plot input_file u 1:15 with lines notitle

unset multiplot


#Bfield time
set multiplot layout 2, 3
unset xtics
unset ytics

# 1
set title "Bx vs time"
set xlabel "timestep"
set ylabel "Bx"
plot bfield u 1:2 with lines notitle

# 2
set title "By vs time"
set xlabel "timestep"
set ylabel "By"
plot bfield u 1:3 with lines notitle

# 3
set title "Bz vs time"
set xlabel "timestep"
set ylabel "Bz"
plot bfield u 1:4 with lines notitle

# 4
set autoscale
set title "Fx vs time"
set xlabel "timestep"
set ylabel "Fx"
plot bfield u 1:5 with lines notitle

# 5
set autoscale
set title "Fy vs time"
set xlabel "timestep"
set ylabel "Fy"
plot bfield u 1:6 with lines notitle

# 6
set autoscale
set title "Fz vs time"
set xlabel "timestep"
set ylabel "Fz"
plot bfield u 1:7 with lines notitle


unset multiplot
