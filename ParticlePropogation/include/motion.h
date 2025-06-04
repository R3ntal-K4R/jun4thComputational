#ifndef INCLUDE_MOTION
#define INCLUDE_MOTION

//find the gamma of the particle
double gamma1(particle p);

//add components in quadrature
double rmag(double x, double y, double z);

//calculates acceleration of particle in x direction from x-boost
double Axx(particle p);

//calculates acceleration of particle in x direction from y-boost
double Axy(particle p);

//calculates acceleration of particle in x direction from z-boost
double Axz(particle p);

//calculates acceleration of particle in y direction from x-boost
double Ayx(particle p);

//calculates acceleration of particle in y direction from y-boost
double Ayy(particle p);

//calculates acceleration of particle in y direction from z-boost
double Ayz(particle p);

//calculates acceleration of particle in z direction from x-boost
double Azx(particle p);

//calculates acceleration of particle in z direction from y-boost
double Azy(particle p);

//calculates acceleration of particle in z direction from z-boost
double Azz(particle p);

//calculate the acceleration of the particle
vector calcAcceleration(particle p);

//update position, velocity, and acceleration of particle
particle update_particle(particle p, eloss e, spacecraft s);

//straggle particle velocity with Fresnel transforms
vector straggle_velocity(vector v, vector prev_v, eloss e);

//calculate velocity in local particle frame
fresnel_vectors fresnel(vector v, vector prev_v);

#endif
