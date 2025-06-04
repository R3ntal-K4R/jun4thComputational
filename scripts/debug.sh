gcc -g -lm main.c -o main
gdb main << EOF
set args Hydrogen_100ev_space.dat Bfield.dat
run
EOF
