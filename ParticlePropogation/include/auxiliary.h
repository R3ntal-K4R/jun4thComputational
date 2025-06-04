#ifndef AUX_INCLUDED
#define AUX_INCLUDED
#ifndef LIBS_INCLUDED
  #define LIBS_INCLUDED
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <string.h>
  #include <assert.h>
  #include <stdbool.h>
  #include <float.h>
  //#include <sys/stat.h> //for directory manipulation
#endif

// define constants
#define SoL 2.997925e8 //speed of light [m/s]
#define Mu0 1.25663706212e-6 //permeability of vacuum [H/m]
#define Epsilon0 8.8541878128e-12 //permitivity of free space [F/m]
#define J2eV 6.242e18 //conversion factor for joules to electron volts
#define Kg2amu 9.223e18 //conversion factor for kilograms to atomic mass units
#define AU2Coul 1.602177e-19 //conversion factor for atomic units (charge) to coulombs
#define SubAtomicNumber2Kg 1.672622e-27 //conversion factor for number of protons/neutrons to kilograms
#define Planck 6.62607004e-34 //planck constant [m^2 kg/s]
#define FineStructure 1/137.036 //fine structure constant [unitless]
#define ZeroTolerance 1e-15
#define DIM_MAX 501 //maximum number of elements in each direction for magnetic field
#define ENABLE_ELOSS //allow for energy loss
#define ENABLE_BFIELD //turn on magnetic field
#define PROPAGATE //propagate particle
#define RELATIVISTICV //calculate velocity relativistically in generator
#define M_PI 3.14159265359 //pi for when usin -std=c99 compiler flag

//define new types: vector, particle, gas, eloss, spacecraft
typedef struct vector {
  double x, y, z;
} vector;

typedef struct fresnel_vectors {
  vector Tangent, Normal, Binormal;
} fresnel_vectors;

typedef struct particles {
  int charge;
  double x, y, z;
  double vx, vy, vz;
  double prev_vx, prev_vy, prev_vz;
  double x0, y0, z0;
  double vx0, vy0, vz0;
  double mass;
  double t, tmax;
  double ax, ay, az;
  double Fx, Fy, Fz;
  double dx, dt, dr;
  double totaldist, elosstotaldist;
  double totaleloss, kinetic_energy;
  double electronic, nuclear; //referring to eloss values, used to validate eloss
  double curr_energy, curr_vel, prev_vel;
  double threshold;
  int stepout, step;
  clock_t tic, toc;
  double bubbleRadius;
  bool insideBubble;
  bool bubbleSwitch;
  int straggle_switch;
  bool run_status; // used to determine if the run is over
} particle;

typedef struct gas {
  double density, atomic_number, atomic_mass;
  double excite_potential;
} gas;

typedef struct spacecrafts {
  int n, dim, bfile_offset;
  double deltaB, B_multiplier;
  double dims, cntr, dimz;
  double (*BFx)[DIM_MAX][DIM_MAX], (*BFy)[DIM_MAX][DIM_MAX], (*BFz)[DIM_MAX][DIM_MAX];
  fpos_t bfile_position;
  double height, radius;
} spacecraft;

typedef struct eloss_arrs {
  double * energy;
  double * electronic;
  double * nuclear;
  double * range;
  double * lateral;
  double * longitudinal;
  double * total;
  int currindex;
  double eloss_multiplier;
} eloss;

typedef struct{
	double *mass;
	double *charge;
	double *energy;
	double *flux;

	double fluxSum;
	int typeNumber;
} DISTRIBUTION;

//check if particle is inside spacecraft
bool isinside(particle p, spacecraft s);

//check if particle is inside bubble for eloss
bool InBubble(particle part, spacecraft space);

//complete program and print statistics
void done(particle p, spacecraft s);

//check if particle is within bounds of magnetic field
bool checkbounds(particle p, spacecraft s);

//log particle output
void logParticle(FILE* outputFile, particle p, double stepeloss);

//log magnetic field that particle is interacting with
void logBField(FILE* outputFile, particle p, vector bfield);

//create a terminal progress bar
void progressBar(double progress);

//calculate dot product of two vectors
double dot(vector a, vector b);

//calculate cross product of two vectors
vector cross(vector a, vector b);

//choose random straggle ranges
vector box_muller(double sigma_lat, double sigma_long);

#endif
