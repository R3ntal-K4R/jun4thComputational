#ifndef INCLUDE_ENERGYLOSS
#define INCLUDE_ENERGYLOSS

//calculate the shell correction for the bethe equation
// double shell_correction(particle p, gas g);

//calculate the density correction to the bethe equation
double density_correction(particle p, gas g);

//calculate the bloch correctoin to the bethe-bloch equation
double bloch_correction(particle p);

//unsure
// double w_max(particle p, gas g);

//calculate the maximum kinetic energy transfer of an incident electron
double Tmax(particle p);

//calculate bethe equation energy loss
double bethe(particle p, gas g);

//calculate beta of particle (v/c)
double betta(particle p);

//interpolate the energy loss file to get energy loss (either elec. or nuclear)
double interp_eloss(eloss e, double currE, int selector);

//calculate the energy loss from synchrotron radiation
vector synchrotron(particle p, vector Bfield);

//calculate the energy lost at a specific time step
vector calcEloss(particle p, gas g, eloss e, double lost_energy);

#endif
