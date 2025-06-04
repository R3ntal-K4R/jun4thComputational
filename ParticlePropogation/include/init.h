#ifndef INCLUDE_INIT
#define INCLUDE_INIT

//"constructor" - create a new particle
particle initializeParticle(const char *flux_file);

//"constructor" - create a new gas
gas initGas(double density, double a_number, double a_mass, double potential);

//read energy loss SRIM file
unsigned long int InitElossFile(char * filename);

//fill and allocate energy loss arrays
void AllocateElossArrays(char * filename, unsigned long int lines, eloss e);

//"constructor" - calls other eloss functions (above) to create new eloss
eloss initializeEloss(particle p, char* elossFile);

//"constuctor" - create new spacecraft
spacecraft initializeBfield(char* bfieldfile);

#endif
