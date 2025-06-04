#include "aux.h"

int main()
{
	FILE *fptr;
	int i;
	int lines=0;
	char file[40], chr;
	char c[1000];

	printf("Enter file name: ");
	scanf("%s", file);

	fptr=fopen(file, "r");
	chr=getc(fptr);
	while(chr!=EOF)//finds the number of lines of the file to determine the length of the variable arrays
	{
		if(chr=='\n')
		{
			lines=lines+1;
//			printf("%i\n", lines);
		}

		chr=getc(fptr);
	}
	fclose(fptr);

	printf("Program found %i lines in file %s.\n", lines, file);//test

	double* energy=malloc((lines-1)*sizeof(double));
	double* electronic=malloc((lines-1)*sizeof(double));
	double* nuclear=malloc((lines-1)*sizeof(double));
	double* total=malloc((lines-1)*sizeof(double));
	double energy_holder, electronic_holder, nuclear_holder, longitude_stragg_holder, lateral_strag_holder;

	fptr=fopen(file, "r");
	fscanf(fptr,"%[^\n]", c); //skips first line
//	printf("Just read 1st line:%s\n",c);


	for(i=0;i<lines;i++)//assigns values to variables
	{
		fscanf(fptr, "%lf %lf %lf %lf %lf", &energy_holder, &electronic_holder, &nuclear_holder, &longitude_stragg_holder, &lateral_strag_holder);
		energy[i]=energy_holder;
		electronic[i]=electronic_holder;
		nuclear[i]=nuclear_holder;
		total[i]=electronic[i]+nuclear[i];
//		printf("scanning line\n");//test
	}

	fclose(fptr);

/*	for(i=0;i<100;i++)
	{
		printf("Energy = %lf\n", energy[i]);
	}
*/
	double input_energy;//energy given by user
	double exact_energy;//interpolated energy

	while(1)//performs interpolation
	{
		printf("Enter an energy: ");
		scanf("%lf", &energy_thing);
		printf("received input energy: %lf\n", energy_thing);
//		printf("Energy[1]=%lf, Energy[0]=%lf\n", energy[1], energy[0]);
		double energy_diff=energy[1]-energy[0];
		int index=(int)energy_thing/(energy_diff);
		double delta_e=fmod(energy_thing,energy_diff);
		printf("energydiff=%lf index=%i delta_e=%lf\n", energy_diff, index, delta_e);
		exact_energy=total[index]+(delta_e/energy_diff)*(total[index+1]-total[index]);
		printf("Energy[index]=%lf Energy[index+1]=%lf\n", total[index], total[index+1]);
		printf("Exact_energy = %lf\n", exact_energy);
	}
}
