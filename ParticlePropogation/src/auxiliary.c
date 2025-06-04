#include <auxiliary.h>
#include <motion.h>

/*checks if particle is inside the spacecraft
  inputs: particle p, spacecraft s
  output: boolean - if particle is in spacecraft habitable zone
*/
bool isinside(particle p, spacecraft s){
  double centerpoint = s.cntr * s.deltaB;
  // printf("Here ");
  // printf(" %lf <= %lf , %lf <= %lf", sqrt((p.x-centerpoint)*(p.x-centerpoint)+(p.y-centerpoint)*(p.y-centerpoint)), (s.dims * s.deltaB), fabs(p.z - centerpoint), (s.dimz * s.deltaB));
  // printf(" %d \n", (sqrt((p.x-centerpoint)*(p.x-centerpoint)+(p.y-centerpoint)*(p.y-centerpoint)) <= (s.dims * s.deltaB) && fabs(p.z - centerpoint) <= (s.dimz * s.deltaB)));
  return (sqrt((p.x-centerpoint)*(p.x-centerpoint)+(p.y-centerpoint)*(p.y-centerpoint)) <= (s.dims * s.deltaB) && fabs(p.z - centerpoint) <= (s.dimz * s.deltaB));
}

/*checks if particle is inside the bubbles
  inputs: particle p, spacecraft space
  output: boolean - if particle is in either top or bottom bubble
*/
bool InBubble(particle part, spacecraft space){
  double centerpoint = space.cntr * space.deltaB;
  double top = space.deltaB * (space.cntr + space.dimz) + part.bubbleRadius; // z-coordinate of top bubble
  double bot = space.deltaB * (space.cntr - space.dimz) - part.bubbleRadius; // z-coordinate of bottom bubble
  // printf("top: %6.3e\n", top);
  // printf("bot: %6.3e\n", bot);
  // exit(1);
  // we need the distance between the center of the bubble and the particle to be less than the radius
  return (rmag(centerpoint - part.x, centerpoint - part.y, top - part.z) <= part.bubbleRadius || \
        rmag(centerpoint - part.x, centerpoint - part.y, bot - part.z) <= part.bubbleRadius);
}

/*completes run and prints run statistics
  inputs: particle p, spacecraft s
  output: none

  prints total runtime, total energy lost, tota distance traveled, and prints initial conditions to 'monte_carlo_output.dat' (commented out)
*/
void done(particle p, spacecraft s){
  clock_t toc = clock();
  double runtime = (double)(toc - p.tic) / CLOCKS_PER_SEC; //convert clocks elapsed to seconds elapsed
  printf("\n\nRuntime:\t%lf\n", runtime);
  // printf("Total eloss (eV)= %4.15f\n", p.totaleloss*J2eV);
  // printf("eloss/KE = %4.15f\n", p.totaleloss / p.kinetic_energy);
  printf("Initial KE (eV)= %4.15f\n", p.kinetic_energy*J2eV);
  // printf("Total Eloss (eV)= %4.15f\n", p.totaleloss*J2eV);
  printf("Final KE (eV)= %4.15f\n", ((gamma1(p) - 1) * p.mass * SoL * SoL)*J2eV);
  printf("Total distance traveled with eloss  = %4.15fm\n", p.elosstotaldist);
  printf("Total distance traveled without eloss  = %4.15fm\n", p.totaldist-p.elosstotaldist);
  printf("Total distance traveled  = %4.15fm\n", p.totaldist);

  //print initial conditions
  // FILE* montecarlo_out = fopen("monte_carlo_output.dat", "a");
  // double beta = (rmag(p.vx0, p.vy0, p.vz0) / SoL);
  // fprintf(montecarlo_out, "%d         %d          %2.5e      %2.5e      %2.5e      %1.5e     %1.5e     %1.5e      %1.5e       %4.4f\n",
  // isinside(p, s), p.charge, p.x0, p.y0, p.z0, p.vx0, p.vy0, p.vz0, beta, (p.totaleloss / p.kinetic_energy));
  // fclose(montecarlo_out);

  //exit program
  printf("Run finished. Live long and prosper! \n");
  // exit(0);
  p.run_status = false; // run is finished for this particle
}

/*checks if particle is in bounds of magnetic field, exits if not
  inputs: particle p, spacecraft s
  output: none

  prints magnetic field dimension and final coordinates of particle if particle exits magnetic field area
*/
bool checkbounds(particle p, spacecraft s){
  char reason[25];
  if(p.x <= 0.0 || p.y <= 0.0 || p.z <= 0.0 || p.x >= s.dim*s.deltaB || p.y >= s.dim*s.deltaB || p.z >= s.dim*s.deltaB){
    if (p.x <= 0.0) {
        strcpy(reason,"p.x <= 0.0");
    } else if (p.y <= 0.0) {
        strcpy(reason,"p.y <= 0.0");
    } else if (p.z <= 0.0) {
        strcpy(reason,"p.z <= 0.0");
    } else if (p.x >= s.dim*s.deltaB) {
        strcpy(reason,"p.x >= s.dim*s.deltaB =");
    } else if (p.y >= s.dim*s.deltaB) {
        strcpy(reason,"p.y >= s.dim*s.deltaB =");
    } else {
        strcpy(reason,"p.z >= s.dim*s.deltaB =");
    }
    printf("Particle has exited the region of interest (Where magnetic field is known).\n Last coordinates (x,y,z) are %14.4e %14.4e %14.4e\n Forces (Fz,Fy,Fz) are %14.4e %14.4e %14.4e\n %s\n",p.x,p.y,p.z,p.Fx,p.Fy,p.Fz,reason);
    printf("dim*delta = %10.10e\n", s.dim*s.deltaB);
    // printf("\nExiting program\n");
    p.run_status = false;
    printf("r: %d\n",p.run_status);
    // done(p, s);
    return false;
  }
  return true;
}

/*log information about the particle to a specified file
  inputs: file pointer outputFile, particle p, spacecraft s
  output: none
*/
void logParticle(FILE* outputFile, particle p, double stepeloss){
  fprintf(outputFile, "%12d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", p.step, p.t, p.x, p.y, p.z, p.vx, p.vy, p.vz, p.ax, p.ay, p.az, rmag(p.vx, p.vy, p.vz), stepeloss*J2eV*p.dr, p.electronic*p.dr, p.nuclear*p.dr);
  // printf("%12d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", p.step, p.t, p.x, p.y, p.z, p.vx, p.vy, p.vz, p.ax, p.ay, p.az, rmag(p.vx, p.vy, p.vz), stepeloss*J2eV);
}

/*log information about the magnetic field the particle interacts with
  inputs: file pointer outputFile, particle p, vector struct of x,y,z components of bfield at the called step
  output: none
*/
void logBField(FILE* outputFile, particle p, vector bfield){
  fprintf(outputFile, "%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n", p.t, bfield.x, bfield.y, bfield.z, p.Fx, p.Fy, p.Fz);
}

/*create a progress bar in the terminal to show completion
  input: percentage complete in decimal
  output: none
*/
void progressBar(double progress){
  int width = 50; //total progress bar width
  int pos = floor(width * progress); //position of completed bar
  printf("["); //print beginning of progress bar
  int i;
  for (i = 0; i < width; ++i){ //print = for each completed step and > for final
    if (i < pos) printf("=");
    else if (i == pos) printf(">");
    else printf(" ");
  }
  printf("] %d %%\r", (int) floor(progress * 100)); //cap progress bar and print percentage
  fflush(stdout); //force terminal printing
  return;
}

/*calculate the dot product of two vectors
  inputs: vector a, vector b
  output: dot product of vectors a and b
*/
double dot(vector a, vector b){
  return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

/*calculate the cross product of two vectors
  inputs: vector a, vector b
  output: cross product of vectors a and b
*/
vector cross(vector a, vector b){
  vector returnVec;
  returnVec.x = a.y * b.z - a.z * b.y;
  returnVec.y = a.z * b.x - a.x * b.z;
  returnVec.z = a.x * b.y - a.y * b.x;
  return returnVec;
}

/*calculate a normal distribution of theta and phi values using Box-Muller method
  inputs: none
  output: vector with random theta as x-value and random phi as y-value
*/
vector box_muller(double sigma_lat, double sigma_long){
  double div, u, v, w, x;
  double mu = 0;
  div = RAND_MAX;
  u = (rand() / div);
  v = (rand() / div);
  w = (rand() / div);
  x = (rand() / div);
  // printf("random numbers: %12.6e %12.6e\n",u,v);
  vector returnVec;
  returnVec.x = (sigma_lat * sqrt(-2 * log(u)) * cos(2 * M_PI * v) + mu);
  returnVec.y = (sigma_long * sqrt(-2 * log(u)) * sin(2 * M_PI * v) + mu);
  returnVec.z = w;
  // returnVec.z = sqrt((sqrt(-2 * log(w)) * sin(2 * M_PI * x) + mu)+(sqrt(-2 * log(u)) * cos(2 * M_PI * v) + mu));
  return returnVec;
}
