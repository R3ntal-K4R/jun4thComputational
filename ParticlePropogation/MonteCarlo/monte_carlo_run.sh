gcc main.c -o main -lm -g
gcc monte_carlo_generator.c -o monte_carlo_generator
gcc monte_carlo_output_creator.c -o monte_carlo_output_creator

echo Enter the amount of runs wanted:  
read MAX

echo 

./monte_carlo_output_creator

COUNTER=0
while [ $COUNTER -lt $MAX ]; do
	let COUNTER=COUNTER+1
	echo Beginning run number $COUNTER
	if [ -e input.dat ]
	then
		rm input.dat
	fi
	./monte_carlo_generator
	./main Hydrogen_SRIM_assigned_spacing.dat Bfield_n400.dat
	echo
	if [ -e output.dat ]
	then
		rm output.dat
	fi
	sleep 1
done