#include <auxiliary.h>
#include <init.h>
#include <magnetic.h>
#include <energy_loss.h>
#include <motion.h>
#include <particleTrajectory.h>
#include <errno.h>

eloss e; //run once at the beginning, use in all calls
gas g;

int main(int argc, char ** argv) {

  // make sure there are enough input files
  if(argc != 4){
      printf("Usage:\t %s elossfile.dat magneticfile.dat fluxfile.dat \n",argv[0]);
      exit(1);
  }

#ifdef ENABLE_BFIELD
  //create a new spacecraft/magnetic field
  printf("Initializing BField\n");
  spacecraft s = initializeBfield(argv[2]);
#endif // END ENABLE_BFIELD

//INITIALIZE PARTICLES
  printf("Initializing particles\n");
  int total_runs, temp, straggle, i;
  FILE *particle_file;
  particle_file=fopen(argv[3],"r");
  fscanf(particle_file, "%d", &total_runs);


  //define particle array based on number of runs from input file
  particle p_array[total_runs]; 
  for (i=0; i<total_runs; i++){
      fscanf(particle_file, "%lf %lf %lf", &p_array[i].x0, &p_array[i].y0, &p_array[i].z0);
      fscanf(particle_file, "%lf %lf %lf", &p_array[i].vx0, &p_array[i].vy0, &p_array[i].vz0);
  }
  //load first particle's info so we can use it in the following particles' definitions
  fscanf(particle_file, "%d %lf", &p_array[0].charge, &p_array[0].mass);
  fscanf(particle_file, "%lf %lf", &p_array[0].dt, &p_array[0].tmax);
  fscanf(particle_file, "%d", &p_array[0].stepout);
  fscanf(particle_file, "%lf", &p_array[0].bubbleRadius);
  fscanf(particle_file,"%lf", &s.B_multiplier);
  fscanf(particle_file, "%d", &temp);
  fscanf(particle_file, "%d", &straggle);
  fclose(particle_file);
  p_array[0].bubbleSwitch = temp;
  p_array[0].straggle_switch = straggle;
  p_array[0].run_status = true;
  if(p_array[0].mass == 0.0 && p_array[0].charge == -1){//electron rest mass is 1836 times smaller than proton
      p_array[0].mass = 0.0005446623; //1/1836
  }
  p_array[0].mass = p_array[0].mass * SubAtomicNumber2Kg;

  for (i=0; i<total_runs; i++){
      p_array[i].charge = p_array[0].charge;
      p_array[i].mass = p_array[0].mass;
      p_array[i].dt = p_array[0].dt;
      p_array[i].tmax = p_array[0].tmax;
      p_array[i].stepout = p_array[0].stepout;
      p_array[i].bubbleRadius = p_array[0].bubbleRadius;
      p_array[i].bubbleSwitch = p_array[0].bubbleSwitch;
      p_array[i].straggle_switch = p_array[0].straggle_switch;
      p_array[i].run_status = true;
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
      p_array[i].ax = 0; p_array[i].ay = 0; p_array[i].az = 0; //initial acceleration, needed to fix seg fault
      if ((rmag(p_array[i].vx, p_array[i].vy, p_array[i].vz) / SoL) >= 1.0){
          printf("\non %d, Input velocity magnitude exceeds c!\n",i);
          printf("Beta is %10.10e\n",(rmag(p_array[i].vx, p_array[i].vy, p_array[i].vz) / SoL) );
          exit(1);
      }
      p_array[i].step = 0;
      p_array[i].totaleloss = 0;
      p_array[i].totaldist = 0; p_array[i].elosstotaldist = 0;
      p_array[i].kinetic_energy = (gamma1(p_array[i]) - 1) * p_array[i].mass * SoL * SoL;
      p_array[i].insideBubble = false;
  }
//END INITIALIZATION OF PARTICLES


  //create a new gas
  // (<density [g/cm^3]>, <atomic number>, <molecular mass [g/mol]>, <mean excitation potential>)
  // want to have this dynamic eventually
  g = initGas(1.13796e-1, 7, 28.0134, 82); // N2

// #ifdef ENABLE_ELOSS //Should probably always do this regardless, as otherwise we can't turn off eloss
  //create new eloss arrays if enabled
  printf("Arguments:\t %s , \t %s ,\t %s  \n",argv[0],argv[1],argv[2]);
  e = initializeEloss(p_array[0], argv[1]); //run once at the beginning, use in all calls
// #endif // END ENABLE_ELOSS

  for (i=0;i<total_runs;i++) {
    propagate(p_array[i],i+1, s);
  }

}

void propagate(particle p, int current_run_number, spacecraft s) {
  #ifdef PROPAGATE
  //print initial conditions and information about the run
  printf("\n||||| RUN #%d |||||\n",current_run_number);
  long long int maxsteps = (long long int) (p.tmax/p.dt);
  printf("Max Steps = %lld\n", maxsteps);
  printf("%12.4e %12.4e %12.4e \n", p.x0, p.y0, p.z0);
  printf("%12.4e %12.4e %12.4e \n", p.vx0, p.vy0, p.vz0);
  printf(" %d %12.4e \n", p.charge, p.mass);
  printf("%12.4e  %12.4e\n", p.dt, p.tmax);
  printf("Outputting every %d steps\n", p.stepout);
  printf("Estimated output file size: %lld kB\n",4*40*(maxsteps)/(p.stepout*1000));
  printf("Beta = %4.4f\n", (rmag(p.vx, p.vy, p.vz) / SoL));
  double kinetic_energy = p.mass * SoL * SoL * (gamma1(p) - 1);  //JOULES
  printf("KE = %8.15f (eV) %4.10e (J)\n", kinetic_energy*J2eV, kinetic_energy);

//more printout testing
  printf("accel init %12.6e %12.6e %12.6e\n",p.ax,p.ay,p.az);
  printf("B-field Multiplier: %4.3e\n", s.B_multiplier);

  int original_stepout;
  FILE *ptfileout; // output file
  char fileout_name[256];  // Assuming the file name won't be longer than 255 characters
  snprintf(fileout_name, sizeof(fileout_name), "./output/%d_output.dat", current_run_number);
  // printf("%s\n",fileout_name);

  ptfileout = fopen(fileout_name, "w");
  if(!ptfileout){
   //handle the error
   printf("something went wrong: %s\n", strerror(errno));
   exit(1);
  }
  // printf("Opened output file successfully\n");

  fprintf(ptfileout,"     timestep     time              x             y              z              vx             vy           \
  vz             ax             ay             az         v_magnitude      ELoss        Electronic      Nuclear\n");
  // printf("printed intial line to output file\n");

  #ifdef ENABLE_BFIELD
  //create output file for magnetic field logging if enabled
  char bfield_out_name[256];  // Assuming the file name won't be longer than 255 characters
  snprintf(bfield_out_name, sizeof(bfield_out_name), "./output/%d_bfield_output.dat", current_run_number);
  // printf("%s\n", bfield_out_name);
  FILE *bfieldout = fopen(bfield_out_name, "w");
  if(!bfieldout){
   //handle the error
   printf("something went wrong: %s\n", strerror(errno));
   exit(1);
  }
  // printf("Opened bfield output file successfully\n");

  fprintf(bfieldout,"       time             b_x              b_y              b_z              F_x              F_y              F_z\n");
  #endif // END ENABLE_BFIELD
  // printf("Got past bfield output in run\n");

  //define variables for propagation
  p.tic = clock();
  bool prev_inside_bubble = false;
  double loss_energy, stepeloss;
  vector magneticField;
  p.curr_vel = rmag(p.vx, p.vy, p.vz);
  printf("Starting propagation...\n");
  p.run_status  = true;

  //log initial position for plotting convenience
  logParticle(ptfileout, p, 0.0);
  for (p.t = 0.0; p.t < p.tmax; p.t += p.dt){

    // printf("run_status: %d\n",p.run_status);
    if (!checkbounds(p,s) || !p.run_status) {
      break; //end current run
    }
    //create/update progress bar
    InBubble(p,s);
    if(isinside(p, s)){
      printf("\nParticle has entered the spacecraft\n");
      // done(p, s);
      break;
    }
    progressBar(p.t/p.tmax);
    p.dr = rmag(p.vx, p.vy, p.vz)*p.dt; //distance traveled this step
#ifdef ENABLE_BFIELD
#ifdef ENABLE_ELOSS
    if (InBubble(p, s)){
      if (!prev_inside_bubble && p.bubbleSwitch){
        original_stepout = p.stepout;
        p.stepout=1; //change stepout for more analysis
      }
      if (!prev_inside_bubble) printf("\nParticle has entered the bubble.\n");
      prev_inside_bubble = true;
      // do eloss propagation
      p.curr_energy = (gamma1(p)-1)*p.mass*SoL*SoL*J2eV; //in eV because of J2eV conversion multiplier
      // printf("curr_energy: %12.6e\n",p.curr_energy);
       while(p.curr_energy < e.energy[e.currindex]){  //Update energy pointer
         e.currindex--;
       }
      magneticField = MagField(p, s);
      // energy loss is different for electrons then for everything else. 
      // loss_energy = 0; // for testing
      loss_energy = interp_eloss(e, p.curr_energy, 0) * e.eloss_multiplier; //eV, for everything but electrons
      p.electronic = interp_eloss(e, p.curr_energy, 1); //eV, electronic eloss for logging, not for electrons
      // printf("elec %12.6e\n",p.electronic);
      p.nuclear = interp_eloss(e, p.curr_energy, 2); //eV, nuclear eloss for logging, not for electrons
      // printf("interped eloss: %12.6e\n",loss_energy);
      vector energy_lost = calcEloss(p, g, e, loss_energy); //joules
      // printf("energy_lost before synch x,y,z %12.6e %12.6e %12.6e \n",energy_lost.x,energy_lost.y,energy_lost.z);
      vector synch_rad = synchrotron(p, magneticField); //joules I think
      // synch_rad.x = 0; synch_rad.y = 0; synch_rad.z = 0; // testing no Eloss
      // printf("synch_rad x,y,z %12.6e %12.6e %12.6e \n",synch_rad.x,synch_rad.y,synch_rad.z);
      energy_lost.x += synch_rad.x * p.dr;
      energy_lost.y += synch_rad.y * p.dr;
      energy_lost.z += synch_rad.z * p.dr;
      stepeloss = rmag(energy_lost.x, energy_lost.y, energy_lost.z);
      // printf("stepeloss: %12.6e \n",stepeloss);
      p.totaleloss += stepeloss;//in Joules
      // printf("total Eloss post synch: %12.6e %12.6e %12.6e \n%12.6e\n",energy_lost.x,energy_lost.y,energy_lost.z,stepeloss);
      //calculate forces on particle
      p.Fx = (p.charge * AU2Coul * ((p.vy * magneticField.z)-(p.vz * magneticField.y))) + energy_lost.x;
      p.Fy = (p.charge * AU2Coul * ((p.vz * magneticField.x)-(p.vx * magneticField.z))) + energy_lost.y;
      p.Fz = (p.charge * AU2Coul * ((p.vx * magneticField.y)-(p.vy * magneticField.x))) + energy_lost.z;
      // printf("FORCE: %12.6e %12.6e %12.6e\n",p.Fx,p.Fy,p.Fz);
      p.elosstotaldist += p.dr;
    } else {
      if (prev_inside_bubble) {
        printf("\nParticle has left the bubble.\n");
        // printf("part: %f %f %f\n", p.x, p.y, p.z);
        p.stepout = original_stepout;
      }
      prev_inside_bubble = false;
      // do free propagation
      magneticField = MagField(p, s);
      vector energy_lost = synchrotron(p, magneticField);
      //calculate forces on particle (no energy loss)
      p.Fx = (p.charge * AU2Coul * ((p.vy * magneticField.z)-(p.vz * magneticField.y))) + energy_lost.x * p.dr;
      p.Fy = (p.charge * AU2Coul * ((p.vz * magneticField.x)-(p.vx * magneticField.z))) + energy_lost.y * p.dr;
      p.Fz = (p.charge * AU2Coul * ((p.vx * magneticField.y)-(p.vy * magneticField.x))) + energy_lost.z * p.dr;
      // printf("ELOSS %12.6e %12.6e %12.6e\n",energy_lost.x,energy_lost.y,energy_lost.z);
      // printf("FORCE: %12.6e %12.6e %12.6e\n",p.Fx,p.Fy,p.Fz);

      //update eloss if particle is electron
      stepeloss = rmag(energy_lost.x*p.dr, energy_lost.y*p.dr, energy_lost.z*p.dr);
      if(p.charge == -1 && p.mass < SubAtomicNumber2Kg){ // synchrotron eloss for electrons
        p.totaleloss += stepeloss;
        p.elosstotaldist += p.dr;
      }
    }
    if(p.step % p.stepout == 0){
      logBField(bfieldout, p, magneticField);
    }
#else // bfield ONLY
    stepeloss = 0.0;
    vector magneticField = MagField(p, s);
    //calculate forces on particle (no energy loss)
    p.Fx = (p.charge * AU2Coul * ((p.vy * magneticField.z)-(p.vz * magneticField.y)));
    p.Fy = (p.charge * AU2Coul * ((p.vz * magneticField.x)-(p.vx * magneticField.z)));
    p.Fz = (p.charge * AU2Coul * ((p.vx * magneticField.y)-(p.vy * magneticField.x)));
    logBField(bfieldout, p, magneticField);
#endif //end ENABLE_ELOSS
#else // ifndef ENABLE_BFIELD
#ifdef ENABLE_ELOSS
    p.curr_energy = (gamma1(p)-1)*p.mass*SoL*SoL*J2eV; //Energy calculated is in eV for search and interpolation
    while(p.curr_energy < e.energy[e.currindex]){
      e.currindex--;
    }
    loss_energy = interp_eloss(e, p.curr_energy, 0) * e.eloss_multiplier;
    energy_lost = calcEloss(p, g, e, loss_energy);
    stepeloss = rmag(energy_lost.x*p.dr, energy_lost.y*p.dr, energy_lost.z*p.dr);
    p.totaleloss += stepeloss;
    p.elosstotaldist += p.dr;
    //calculate forces on particle (no magnetic field)
    p.Fx = energy_lost.x;
    p.Fx = energy_lost.y;
    p.Fx = energy_lost.z;
#else // ifndef ENABLE_ELOSS
    stepeloss = 0.0;
    p.Fx = 0;
    p.Fy = 0;
    p.Fz = 0;
#endif // end ENABLE_ELOSS
#endif // end ENABLE_BFIELD
    //update particle and log if necessary
    // printf("to the update %d\n",p.step);
    p = update_particle(p, e, s);
    if (p.step % p.stepout == 0){
      logParticle(ptfileout, p, stepeloss);
    }
  } //end time loop
  fclose(ptfileout);
#endif
  done(p, s);
  printf("Finished with run %d\n",current_run_number);
}
