set terminal gif animate delay 5 loop 0 optimize
set output "rot.gif"



n = 100
do for [i=1:n] {
   set view 60, i*360/n
   splot 'output.dat' every 100 u 2:3:4 w l notitle
}

set output
