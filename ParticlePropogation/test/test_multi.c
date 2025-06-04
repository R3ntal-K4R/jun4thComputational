#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
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
  bool bubbleswitch;
  int straggle_switch;
//   int run_status; // used to determine if the run is over
} particle;
double gamma1(particle p) {
  double v_magnitude2 = ((p.vx * p.vx)+(p.vy * p.vy)+(p.vz * p.vz));
  double c2 = SoL*SoL;
  double Beta2 = v_magnitude2 / (c2);
  double g = 1.0e0 / pow((1.0e0 - Beta2), 0.5);
  return g;
}
double rmag(double x, double y, double z) {
    return pow(x * x + y * y + z * z, .5);
}

int main(int argc, char ** argv) {
    //create a new particle
    int total_runs, temp;
    FILE *particle_file;
    particle_file=fopen(argv[3],"r");
    fscanf(particle_file, "%d", &total_runs);

    // struct particle *particle_array;
    // //define particle array based on number of runs from input file
    // particle_array=malloc(sizeof(particle)*total_runs);
    particle p_array[total_runs]; 
    int i;
    for (i=0; i<total_runs; i++){
        fscanf(particle_file, "%lf %lf %lf", &p_array[i].x0, &p_array[i].y0, &p_array[i].z0);
        fscanf(particle_file, "%lf %lf %lf", &p_array[i].vx0, &p_array[i].vy0, &p_array[5].vz0);
    }
    fscanf(particle_file, "%d %lf", &p_array[0].charge, &p_array[0].mass);
    fscanf(particle_file, "%d", &p_array[0].stepout);
    fscanf(particle_file, "%lf", &p_array[0].bubbleRadius);
    fscanf(particle_file, "%d", &temp);
    p_array[0].bubbleswitch = temp;
    fscanf(particle_file, "%d", &p_array[0].straggle_switch);
    fclose(particle_file);
    // p_array[0].run_status = 1;
    if(p_array[0].mass == 0.0 && p_array[0].charge == -1){//electron rest mass is 1836 times smaller than proton
        p_array[0].mass = 0.0005446623; //1/1836
    }

    for (i=0; i<total_runs; i++){
        p_array[i].charge = p_array[0].charge;
        p_array[i].mass = p_array[0].mass;
        p_array[i].stepout = p_array[0].stepout;
        p_array[i].bubbleRadius = p_array[0].bubbleRadius;
        p_array[i].bubbleswitch = p_array[0].bubbleswitch;
        p_array[i].straggle_switch = p_array[0].straggle_switch;
        // p_array[i].run_status = 1;
        p_array[i].mass = p_array[i].mass * SubAtomicNumber2Kg;
        //set initial conditions
        p_array[i].x = p_array[i].x0;
        p_array[i].y = p_array[i].y0;
        p_array[i].z = p_array[i].z0;
        p_array[i].vx = p_array[i].vx0;
        p_array[i].vy = p_array[i].vy0;
        p_array[i].vz = p_array[i].vz0;
        p_array[i].prev_vx = p_array[i].vx0;
        p_array[i].prev_vy = p_array[i].vy0;
        p_array[i].prev_vz = p_array[i].vz0;
        p_array[i].prev_vel = rmag(p_array[i].vx0, p_array[i].vy0, p_array[i].vz0);
        if ((rmag(p_array[i].vx, p_array[i].vy, p_array[i].vz) / SoL) >= 1.0){
            printf("\non %d, Input velocity magnitude exceeds c!\n",i);
            printf("Beta is %10.10e\n",(rmag(p_array[i].vx, p_array[i].vy, p_array[i].vz) / SoL) );
            exit(1);
        }
        p_array[i].step = 0;
        p_array[i].totaleloss = 0;
        p_array[i].kinetic_energy = (gamma1(p_array[i]) - 1) * p_array[i].mass * SoL * SoL;
        p_array[i].insideBubble = false;
    }

    //TESTING
    for(i=0;i<total_runs;i++){
        printf("%d particle initial: %lf %lf %lf \n", i, p_array[i].x, p_array[i].y, p_array[i].z);
    }
}