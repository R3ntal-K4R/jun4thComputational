#include <auxiliary.h>
#include <generator.h>

//returns mass and charge in units given to distribution file
//fluxes must all be provided in same units
//energy must be in joules (changing it to take eV would be a one line thing, very simple)

DISTRIBUTION *distFromFile(const char *filename){
	DISTRIBUTION *d = malloc(sizeof(DISTRIBUTION));
	FILE *f;
	f=fopen(filename,"r");
	if (f == NULL) {
      printf("Error: Could not find flux file %s. \n", filename);
      exit(1);
  }
	char buffer[200];
	fscanf(f," %[^\n] ",buffer);
	int num1;
	sscanf(buffer,"%d",&num1);
	d->typeNumber=num1;
	d->mass=malloc(sizeof(double)*num1);
	d->charge=malloc(sizeof(double)*num1);
	d->energy=malloc(sizeof(double)*num1);
	d->flux=malloc(sizeof(double)*num1);

	int num2;
	int num3;
	int num4;
	int cnt=0;
	double fluxSum=0;
	while(fscanf(f," %[^\n] ",buffer)!=EOF){
		sscanf(buffer,"%d %d %d %d",&num1,&num2,&num3,&num4);
		d->mass[cnt]= num1;
		d->charge[cnt]= num2;
		d->energy[cnt]= num3;
		d->flux[cnt]= num4;
		fluxSum+= num4;
		cnt++;
	}
	d->fluxSum=fluxSum;

	fclose(f);
	return d;
}

particle generateParticle(const char *filename){
	// FILE *f;
	// f=fopen(filename,"w");

	DISTRIBUTION *d = distFromFile(filename);

	time_t t;
	srand((unsigned) time(&t));
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

	particle p;

	//need to randomize these
	double sx = 0; double sy = 0; double sz = 0;
	double distance = 1;

	double theta = rand() * 360.0 / RAND_MAX;
	double phi = (double)rand() * 180.0 / (double)RAND_MAX;
	// double x,y,z;
	p.z0 = sz+cos(phi*M_PI/180.0)*distance;
	p.y0 = sy+sin(phi*M_PI/180.0)*distance*sin(theta*M_PI/180.0);
	p.x0 = sx+sin(phi*M_PI/180.0)*distance*cos(theta*M_PI/180.0);

	// fprintf(f, "%f %f %f\n",x,y,z);

	double thetaV = (double)rand() * 180.0 / (double)RAND_MAX + theta + 90.0;
	double phiV = (double)rand() * 180.0 / (double)RAND_MAX + phi + 90.0;
	// double vx,vy,vz;
	#ifdef RELATIVISTICV
	double v = sqrt((pow(d->energy[pick],2)-pow(d->mass[pick],2)*pow(SoL,4))/pow(SoL,2))/d->mass[pick];
	#else
	double v = sqrt(d->energy[pick]*2/d->mass[pick]);
	#endif
	p.vz0 = cos(phi*M_PI/180.0)*v;
	p.vy0 = sin(phi*M_PI/180.0)*v*sin(theta*M_PI/180.0);
	p.vx0 = sin(phi*M_PI/180.0)*v*cos(theta*M_PI/180.0);

	printf("\n%lf \t %lf\n", d->mass[pick], d->mass[pick]);
	p.mass = d->mass[pick];
	p.charge = d->charge[pick];

	return p;

	// fprintf(f, "%f %f %f\n",vx,vy,vz);
	//
	// fprintf(f,"%f %f %f\n",d->mass[pick],d->charge[pick],d->energy[pick]);

	// fclose(f);
}
