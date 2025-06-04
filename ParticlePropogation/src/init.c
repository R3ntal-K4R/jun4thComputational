#include <auxiliary.h>
#include <init.h>
#include <motion.h>
#include <generator.h>

/*creates a new particle by filling in values from the input file
  inputs: none
  outputs: particle struct

  reads in the input file and fills the necessary variables with the corresponding values
*/
particle initializeParticle(const char *flux_file){

  // particle p = generateParticle(flux_file);
  //read input file
  int total_runs;
  particle p;
  int temp,straggle;
  FILE *ptfilein;
	ptfilein=fopen(flux_file,"r");
  fscanf(ptfilein, "%d", &total_runs);
  fscanf(ptfilein, "%lf %lf %lf", &p.x0, &p.y0, &p.z0);
  fscanf(ptfilein, "%lf %lf %lf", &p.vx0, &p.vy0, &p.vz0);
  fscanf(ptfilein, "%d %lf", &p.charge, &p.mass);
  fscanf(ptfilein, "%lf %lf", &p.dt, &p.tmax);
  fscanf(ptfilein, "%d", &p.stepout);
  fscanf(ptfilein, "%lf", &p.bubbleRadius);
  fscanf(ptfilein, "%d", &temp);
  fscanf(ptfilein, "%d", &straggle);
  fclose(ptfilein);
  //convert mass to Kg from # subatomic particles
  // p.dt = 1.0e-13;
  // p.tmax = 1.0e-7;
  // p.stepout = 500;

  if(p.mass == 0.0 && p.charge == -1){//electron rest mass is 1836 times smaller than proton
    p.mass = 0.0005446623; //1/1836
  }
  p.mass = p.mass * SubAtomicNumber2Kg;
  //set initial conditions
  p.x = p.x0;
  p.y = p.y0;
  p.z = p.z0;
  p.vx = p.vx0;
  p.vy = p.vy0;
  p.vz = p.vz0;
  p.prev_vx = p.vx0;
  p.prev_vy = p.vy0;
  p.prev_vz = p.vz0;
  p.prev_vel = rmag(p.vx0, p.vy0, p.vz0);
  if ((rmag(p.vx, p.vy, p.vz) / SoL) >= 1.0){
      printf("\nInput velocity magnitude exceeds c!\n");
      printf("Beta is %10.10e\n",(rmag(p.vx, p.vy, p.vz) / SoL) );
      exit(1);
  }
  p.step = 0;
  p.totaleloss = 0;
  p.kinetic_energy = (gamma1(p) - 1) * p.mass * SoL * SoL;
  p.insideBubble = false;
  p.bubbleSwitch = temp;
  p.straggle_switch = straggle;
  // p.run_status = 1;
  return p;
}

/*creates a new gas based on passed parameters
 inputs: bethe constant, gas density, atomic number, atomic mass, excitation potential, charge
 outputs: gas g

 want to just be able to specify the name and fill the rest automatically eventually
*/
gas initGas(double density, double a_number, double a_mass, double potential){
  gas g;
  g.density = density; // in g/cm^3
  g.atomic_number = a_number;
  g.atomic_mass = a_mass;
  g.excite_potential = potential;
  return g;
}

/*gets the number of lines from the SRIM file
  inputs: SRIM file name
  outputs: number of lines in file

  can maybe move into AllocateElossArrays function
*/
// popen doesn't exist in c99, and petridis doesn't want us to use it
// unsigned long int InitElossFile(char * filename){
//   unsigned long int lines;
//   char eloss_wc_cmd[100] = "wc -l ";
//   strcat(eloss_wc_cmd,filename);
//   FILE * eloss_wc_file = popen(eloss_wc_cmd,"r");
//   fscanf(eloss_wc_file,"%ld",&lines);
//   //printf("%ld",lines);
//   return lines;
// }

unsigned long int InitElossFile(char *filename) {
    unsigned long int lines = 0;
    char buffer[256];  // Adjust the size of the buffer as needed
    FILE *eloss_file = fopen(filename, "r");

    if (eloss_file == NULL) {
        // Handle file opening error, for example, print an error message and return an error code.
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 0;  // Returning 0 lines as an error indicator
    }

    while (fgets(buffer, sizeof(buffer), eloss_file) != NULL) {
        lines++;
    }

    fclose(eloss_file);
    // printf("%ld",lines);
    return lines;
}


/*allocates and fills energy loss arrays
  inputs: SRIM file name, number of lines in file, eloss struct e
  output: none

  should maybe add this into initializeEloss instead of a separate function. may have problems with accurately accessing values otherwise
*/
void AllocateElossArrays(char * filename, unsigned long int lines, eloss e){
  double energy_holder, electronic_holder, nuclear_holder;
  double range, lateral, longitudinal;
  FILE * elossfile = fopen(filename, "r");
  char chr;
  printf("ElossFile In Allocation:\t %s  \n",filename);
  //fill arrays with corresponding values
  while((chr = getc(elossfile)) != '\n'){}
    for(int i=0; i<lines; i++){
      fscanf(elossfile, "%lf %lf %lf %lf %lf %lf", &energy_holder, &electronic_holder, &nuclear_holder,\
            &range, &lateral, &longitudinal);//
      e.energy[i] = energy_holder;
      e.electronic[i] = electronic_holder;
      e.nuclear[i] = nuclear_holder;
      e.range[i] = range;
      e.lateral[i] = lateral;
      e.longitudinal[i] = longitudinal;
      //printf("\t %f \t %f \t %f \t %f \t %f \t %f \n",energy_holder,electronic_holder,nuclear_holder,range,lateral,longitudinal);
      e.total[i] = electronic_holder+nuclear_holder;
    }
  fclose(elossfile);
}

/*creates a new eloss struct and fills values
  inputs: particle p, SRIM file name
  outputs: eloss struct e

  may need to combine all eloss functions into one to precent memory access issues
*/
eloss initializeEloss(particle p, char* elossFile){
  eloss e;
  if (p.charge == -1){ //SRIM doesn't deal with electrons
    e.energy=malloc(sizeof(double));
    e.electronic=malloc(sizeof(double));
    e.nuclear=malloc(sizeof(double));
    e.total=malloc((sizeof(double)));
  } else {
    printf("ElossFile\t %s  \n",elossFile);
    unsigned long int lines = InitElossFile(elossFile);
    e.energy=malloc(lines*sizeof(double));
    e.electronic=malloc(lines*sizeof(double));
    e.nuclear=malloc(lines*sizeof(double));
    e.range=malloc(lines*sizeof(double));
    e.lateral=malloc(lines*sizeof(double));
    e.longitudinal=malloc(lines*sizeof(double));
    e.total=malloc(lines*sizeof(double));

    e.eloss_multiplier = 1.0;
    AllocateElossArrays(elossFile, lines, e);
    double initial_energy = (gamma1(p)-1)*p.mass*SoL*SoL*J2eV;// multiplied by 6.242e18 to get energy in [eV] for searching the file If ev *J2eV
    e.currindex = lines-1;
    while (initial_energy < e.energy[e.currindex]){ //Set pointer to Eloss array to point just before the current energy
        e.currindex--;
    }
    if (initial_energy<e.energy[0]){
      char c1[100];
      sprintf(c1, "%f", e.energy[0]);
      printf("Min Energy: %s\n", c1);
      char c2[100];
      sprintf(c2, "%f", initial_energy);
      printf("initial_energy: %s\n", c2);
      printf("Error: array out of bounds - input energy too low\n");
      exit(1);
    } else if (initial_energy>e.energy[lines-1]) {
      printf("Energy: %12.4e\n", initial_energy);
      printf("Max Energy: %12.4e\n", e.energy[lines-1]);
      printf("Error: array out of bounds - input energy too high\n");
      exit(1);
    }
  }
  return e;
}

/*creates a new spacecraft struct to hold magnetic field information
  inputs: name of magnetic field file
  outputs: spacecraft struct

  may change what is included in the spacecraft struct for ease of understanding
*/
spacecraft initializeBfield(char* bfieldfile){
  FILE * bfile;
  bfile = fopen(bfieldfile,"r");
  spacecraft s;
  // s.height = 10;
  // s.radius = 14.1421;
  s.B_multiplier = 1; //not changed until particle initialization
  fscanf(bfile, "%d",&s.n);
  s.dim = s.n-1;
  if (s.n > DIM_MAX){ //make sure magnetic field isn't too big
    printf("Need more allowed values for B-field...\n ");
    printf("Current max: %d\n ", DIM_MAX);
    printf("Actual: %d\n ", s.dim);
    printf("Change dim_max in aux.h ");
    exit(1);
  }
  fscanf(bfile, "%lf %lf %lf %lf",&s.deltaB,&s.dims,&s.cntr,&s.dimz);
  // s.deltaB = 0.1;

  //allocate memory of magnetic field arrays
  double static (* BFx)[DIM_MAX][DIM_MAX];
  double static (* BFy)[DIM_MAX][DIM_MAX];
  double static (* BFz)[DIM_MAX][DIM_MAX];
  BFx = malloc((DIM_MAX) * sizeof(*BFx));
  BFy = malloc((DIM_MAX) * sizeof(*BFy));
  BFz = malloc((DIM_MAX) * sizeof(*BFz));

  int i,j,k;
  double currBx,currBy,currBz;
  fgetpos(bfile, &s.bfile_position);
  fscanf(bfile,"%d",&s.bfile_offset);
  if(s.bfile_offset != 0){
    printf("The indices of %s must start with 0!\n",bfieldfile);
    exit(1);
  }
  fsetpos(bfile, &s.bfile_position);
  //fill magnetic field arrays
  while(fscanf(bfile,"%d %d  %d %lf %lf %lf",&i,&j,&k,&currBx,&currBy,&currBz) != EOF){
    BFx[i][j][k] = currBx;
    BFy[i][j][k] = currBy;
    BFz[i][j][k] = currBz;
  }
  s.BFx = BFx; s.BFy = BFy; s.BFz = BFz;//do not unallocate. the struct references these values
  fclose(bfile);
  printf("Dim Max: %d\n",DIM_MAX);
  return s;
}
