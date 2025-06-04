#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include <math.h>
#include<string.h> 

//#include "main.h"

#define PI 3.14159265
#define C 299792458
//global constants
FILE *output_file;
int run_number;
time_t t; //used for random calls

//returns mass and charge in units given to distribution file
//fluxes must all be provided in same units
//energy must be in joules (changing it to take eV would be a one line thing, very simple)


typedef struct{
	double *mass;
	int *charge;
	double *energy;
	double *flux;
	
	double fluxSum;
	int typeNumber;
} DISTRIBUTION;


DISTRIBUTION *distFromFile(const char *filename);
void generateParticle(DISTRIBUTION *d, double distance,double sx,double sy,double sz,double regulator,int choice,int current_run);


DISTRIBUTION *distFromFile(const char *filename){
	DISTRIBUTION *d = malloc(sizeof(DISTRIBUTION));
	FILE *f;
	f=fopen(filename,"r");
	char buffer[200];
	fscanf(f," %[^\n] ",buffer);
	int num;
	sscanf(buffer,"%d",&num);
	d->typeNumber=num;
	d->energy=malloc(sizeof(double)*num);
	d->mass=malloc(sizeof(double)*num);
	d->charge=malloc(sizeof(double)*num);
	d->flux=malloc(sizeof(double)*num);
	float num1;
	float num2; // for mass
	int num3;
	float num4;
	int cnt=0;
	double fluxSum=0;
	int count1 = 0;
	while(fscanf(f," %[^\n] ",buffer)!=EOF){



		sscanf(buffer,"%f %f %d %f",&num1,&num2,&num3,&num4);

		//----------------------------Input Testing----------------------------
		//printf("%s\n", buffer);
		//fscanf(f, "%s %f %d %d %f %f %f", &name, &energy, &mass, &charge, &aveflux, &minflux, &maxflux);
		//count1 ++;
		//char c1[100];
    	//sprintf(c1, "%d", count1);
    	//printf("Count: %s\n", c1);


		//char mass1[100];
    	//sprintf(mass1, "%f", num2);
    	//printf("mass: %s\n", mass1);

		//char charge1[100];
    	// sprintf(charge1, "%f", num3);
    	// printf("charge: %s\n", charge1);

		// char energy1[100];
    	// sprintf(energy1, "%f", num1);
    	// printf("energy: %s\n", energy1);

		// char flux1[100];
    	// sprintf(flux1, "%f", num4);
    	// printf("flux: %s\n", flux1);
		//--------------------------------------------------------

		d->energy[cnt]= num1;
		d->mass[cnt]= num2;
		d->charge[cnt]= num3;
		d->flux[cnt]= num4;
		fluxSum+= num4;
		cnt++;
	}
	d->fluxSum=fluxSum;
	
	fclose(f);
	return d;
}

void generateParticle(DISTRIBUTION *d, double distance,double sx,double sy,double sz, double regulator, int choice, int current_run){



	FILE *f;
	f=output_file;
	
	double draw = (double)rand() * d->fluxSum / (double)RAND_MAX;
	
	double cur=0;
	int pick = 0;
	while(1){
		cur+=d->flux[pick];
		if(draw<=cur){
			break;
		}
		pick++;
	}
	
	
	double theta = rand() * 360.0 / RAND_MAX;
	double phi = (double)rand() * 180.0 / (double)RAND_MAX;
	double x,y,z;
	// z = sz+cos(phi*PI/180.0)*distance+25;
	// y = sy+sin(phi*PI/180.0)*distance*sin(theta*PI/180.0)+25;
	// x = sx+sin(phi*PI/180.0)*distance*cos(theta*PI/180.0)+25;
	if (choice == 1) { // const x axis
		z = sz+cos(phi*PI/180.0)*distance+regulator;
		y = sy+sin(phi*PI/180.0)*distance*sin(theta*PI/180.0)+regulator;
		x = 0.1;
	} else if (choice == 2) { // const y axis
		z = sz+cos(phi*PI/180.0)*distance+regulator;
		y = 0.1;
		x = sx+sin(phi*PI/180.0)*distance*cos(theta*PI/180.0)+regulator;
	} else if (choice == 3) { // const z axis
		z = 0.1;
		y = sy+sin(phi*PI/180.0)*distance*sin(theta*PI/180.0)+regulator;
		x = sx+sin(phi*PI/180.0)*distance*cos(theta*PI/180.0)+regulator;
	} else {
		z = sz+cos(phi*PI/180.0)*distance+regulator;
		y = sy+sin(phi*PI/180.0)*distance*sin(theta*PI/180.0)+regulator;
		x = sx+sin(phi*PI/180.0)*distance*cos(theta*PI/180.0)+regulator;
	}
	
	
	fprintf(f, "%f %f %f\n",x,y,z);
	
	
	double thetaV = (double)rand() * 180.0 / (double)RAND_MAX + theta + 90.0;
	double phiV = (double)rand() * 180.0 / (double)RAND_MAX + phi + 90.0;
	double vx,vy,vz;
	double cenergy  = d->energy[pick];
//    printf("InputEnergy: %12.4e\n", cenergy);
	double cmass1 = d->mass[pick]; // WRONG ???
//    printf("cmass1: %12.4e\n", cmass1);
	double cmass = cmass1 * (1.66054* pow(10,-27));
//    printf("InputMass: %12.4e\n", cmass);
	double restMass = cmass * pow(C,2);
//    printf("RestMass: %12.4e\n", restMass);
	double beforp = (pow(cenergy,2)-(pow(cmass,2)*pow(C,4)));
//    printf("before P: %12.4e\n", beforp);
	double p = sqrt((pow(cenergy,2)-(pow(cmass,2)*pow(C,4))))/C;
//    printf("p: %12.4e\n", p);
    
	double gamma = cenergy/(cmass * pow(C,2));

//    printf("gamma: %12.4e\n", gamma);

	double v = p/(gamma*cmass);
//    printf("v: %12.4e\n", v);


	double velo2 = (1-pow((cmass*pow(C,2)/cenergy),2));//C*sqrt(1-pow((cmass*pow(C,2)/cenergy),2));
	// printf("Other Velo: %12.4e\n",velo2);
	double OutputEnergy = cmass * C * C* (gamma);
	// printf("Output Energy: %12.4e\n",OutputEnergy);

	vz = cos(phi*PI/180.0)*v;
	vy = sin(phi*PI/180.0)*v*sin(theta*PI/180.0);
	vx = sin(phi*PI/180.0)*v*cos(theta*PI/180.0);
	
	fprintf(f, "%f %f %f\n",vx,vy,vz);

	double max_t = 10/v;
	max_t = fmin(max_t,1e-3);
	
	if (current_run == run_number-1) {
	fprintf(f,"%d %12.3e\n",d->charge[pick],d->mass[pick]);  // fprintf(f,"%f %f %f\n",d->mass[pick],d->charge[pick],d->energy[pick]);
	fprintf(f,"1.0e-11 	%4.3e \n", max_t);
	fprintf(f,"8000 \n5\n1.000\n1\n0\n");
	fprintf(f,"//x y z\n");
	fprintf(f,"//vx vy vz\n");
	fprintf(f,"//charge, mass (in number of subatomic particles)\n");
	fprintf(f,"//dt, tmax\n");
	fprintf(f,"//stepout\n");
	fprintf(f,"//radius of bubble\n");
	fprintf(f,"//B-field multiplier\n");
	fprintf(f,"//switch to change stepout when inside bubble\n");
	fprintf(f,"//straggling switch (1 for straggling, else for no straggling)\n");
	fprintf(f,"//Note: USE SMALLER dt FOR ELECTRONS!\n");
	}
}



void main(int argc, char **argv){
	printf("running main\n");

	if (argc < 3) {
        printf("Error: arguments must be <radiation file> <number of runs>.\n");
        exit(1);
    }

    const int inputname_size = 256;
    char inputname[inputname_size];
    strcpy(inputname, argv[1]); // radiation data file
    DISTRIBUTION *d = distFromFile(inputname);
	run_number = atoi(argv[2]); // number of runs to perform
    double distance = 20; //Noah's value for 500 cubed field
	double distance_small = 10; // Value for 200 cubed field
    double sx = 1;
    double sy = 1;
    double sz = 1;
	double regulator = 25; //Noah's value for 500 cubed field
	double regulator_small = 10; // Value for 200 cubed field
	int choice_x = 1;
	srand((unsigned) time(&t)); //seed random generator

	const char filename[] = "../../input.dat";
	output_file = fopen(filename, "w"); //open the file
	fprintf(output_file,"%d\n",run_number); //write the initial run number
	printf("%d\n",run_number);

	int i;
	for (i=0; i<run_number; i++) {
		/// !!
		generateParticle(d, distance_small, sx, sy, sz, regulator_small, choice_x, i); //change regulator in this call to account for different B-field size
		/// !!
	}


	fclose(output_file); //close file at the end of the program
    }



