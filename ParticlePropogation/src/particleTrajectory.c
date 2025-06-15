#include <auxiliary.h>
#include <init.h>
#include <magnetic.h>
#include <energy_loss.h>
#include <motion.h>
#include <particleTrajectory.h>
#include <errno.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <sys/stat.h> 
#include <sys/types.h> 


eloss e; //run once at the beginning, use in all calls
gas g;

void propagate(particle p, int current_run_number, spacecraft s, const char* output_dir_path);

int main(int argc, char ** argv) {

  // make sure there are enough input files
  if(argc != 4){
      printf("Usage:\t %s elossfile.dat magneticfile.dat fluxfile.dat \n",argv[0]);
      exit(1); // 
  }

  // Get simulation name from user
  char simulation_name[256];
  char output_dir[512];

  printf("Enter the name for this simulation: ");
  if (scanf("%255s", simulation_name) != 1) {
      fprintf(stderr, "Error reading simulation name.\n");
      exit(1);
  }
  // Clear the input buffer after scanf
  int c;
  while ((c = getchar()) != '\n' && c != EOF);

  // Create the path for the new directory
  snprintf(output_dir, sizeof(output_dir), "outputs/%s", simulation_name);

  // Create the simulation-specific output directory
  // For POSIX systems (Linux, macOS):
  if (mkdir(output_dir, 0777) == -1) {
      if (errno != EEXIST) { // Don't treat "directory already exists" as a fatal error
          perror("Error creating simulation directory");
          // exit(1); // You might want to exit if creation fails for other reasons
      }
  }

  printf("Output will be saved in: %s\n", output_dir);

  // Create the log directory
  char log_dir_path[512];
  snprintf(log_dir_path, sizeof(log_dir_path), "%s/log", output_dir);
  if (mkdir(log_dir_path, 0777) == -1 && errno != EEXIST) {
      perror("Error creating log directory");
      exit(1);
  }

  // Redirect stdout to a file in the log directory
  char log_file_path[512];
  snprintf(log_file_path, sizeof(log_file_path), "%s/terminal_output.log", log_dir_path);
  if (freopen(log_file_path, "w", stdout) == NULL) {
      fprintf(stderr, "Error redirecting stdout to file: %s\n", strerror(errno));
      exit(1);
  }
  setbuf(stdout, NULL); // Disable buffering


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
      p_array[0].mass = 0.0005446623; 
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
      p_array[i].ax = 0; p_array[i].ay = 0; p_array[i].az = 0; // initial acceleration, need to fix seg fault
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
  g = initGas(1.13796e-1, 7, 28.0134, 82); 

  printf("Arguments:\t %s , \t %s ,\t %s  \n",argv[0],argv[1],argv[2]); 
  e = initializeEloss(p_array[0], argv[1]); 

  for (i=0;i<total_runs;i++) {
    // Pass the output_dir to the propagate function
    propagate(p_array[i], i+1, s, output_dir); 
  }
  return 0; 
}

void propagate(particle p, int current_run_number, spacecraft s, const char* output_dir_path) {
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
  double kinetic_energy = p.mass * SoL * SoL * (gamma1(p) - 1); 
  printf("KE = %8.15f (eV) %4.10e (J)\n", kinetic_energy*J2eV, kinetic_energy); 

  printf("accel init %12.6e %12.6e %12.6e\n",p.ax,p.ay,p.az); 
  printf("B-field Multiplier: %4.3e\n", s.B_multiplier); 

  int original_stepout = p.stepout;
  FILE *ptfileout; // output file
  char fileout_name[512]; 
  snprintf(fileout_name, sizeof(fileout_name), "%s/%d_output.dat", output_dir_path, current_run_number); 
  ptfileout = fopen(fileout_name, "w"); 
  if(!ptfileout){
   printf("something went wrong opening %s: %s\n", fileout_name, strerror(errno)); 
   exit(1); 
  }

  fprintf(ptfileout,"     timestep     time              x             y              z              vx             vy           \
  vz             ax             ay             az         v_magnitude      ELoss        Electronic      Nuclear\n");
  // printf("printed intial line to output file\n");

  #ifdef ENABLE_BFIELD
  char bfield_out_name[512]; 
  snprintf(bfield_out_name, sizeof(bfield_out_name), "%s/%d_bfield_output.dat", output_dir_path, current_run_number);  
  FILE *bfieldout = fopen(bfield_out_name, "w"); 
  if(!bfieldout){
   printf("something went wrong opening %s: %s\n", bfield_out_name, strerror(errno)); 
   exit(1); 
  }
  fprintf(bfieldout,"       time             b_x              b_y              b_z              F_x              F_y              F_z\n");
  #endif // END ENABLE_BFIELD
  // printf("Got past bfield output in run\n");

  p.tic = clock(); 
  bool prev_inside_bubble = false; 
  double loss_energy, stepeloss;
  vector magneticField;
  p.curr_vel = rmag(p.vx, p.vy, p.vz);
  printf("Starting propagation...\n");
  p.run_status  = true; 
  logParticle(ptfileout, p, 0.0); 
  for (p.t = 0.0; p.t < p.tmax; p.t += p.dt){ 

    if (!checkbounds(p,s) || !p.run_status) { 
      break; 
    }
    InBubble(p,s); 
    if(isinside(p, s)){ 
      printf("\nParticle has entered the spacecraft\n"); 
      break; 
    }
    progressBar(p.t/p.tmax); 
    p.dr = rmag(p.vx, p.vy, p.vz)*p.dt; 
#ifdef ENABLE_BFIELD
#ifdef ENABLE_ELOSS
    if (InBubble(p, s)){ 
      if (!prev_inside_bubble && p.bubbleSwitch){ 
        original_stepout = p.stepout; 
        p.stepout=1;  
      }
      if (!prev_inside_bubble) printf("\nParticle has entered the bubble.\n"); 
      prev_inside_bubble = true; 
      p.curr_energy = (gamma1(p)-1)*p.mass*SoL*SoL*J2eV; 
      while(p.curr_energy < e.energy[e.currindex]){  
         e.currindex--; 
      }
      magneticField = MagField(p, s); 
      loss_energy = interp_eloss(e, p.curr_energy, 0) * e.eloss_multiplier; 
      p.electronic = interp_eloss(e, p.curr_energy, 1); 
      p.nuclear = interp_eloss(e, p.curr_energy, 2); 
      vector energy_lost = calcEloss(p, g, e, loss_energy); 
      vector synch_rad = synchrotron(p, magneticField); 
      synch_rad.z = 0; 
      energy_lost.x += synch_rad.x * p.dr; 
      energy_lost.y += synch_rad.y * p.dr; 
      energy_lost.z += synch_rad.z * p.dr; 
      stepeloss = rmag(energy_lost.x, energy_lost.y, energy_lost.z); 
      p.totaleloss += stepeloss; 
      p.Fx = (p.charge * AU2Coul * ((p.vy * magneticField.z)-(p.vz * magneticField.y))) + energy_lost.x; 
      p.Fy = (p.charge * AU2Coul * ((p.vz * magneticField.x)-(p.vx * magneticField.z))) + energy_lost.y;  
      p.Fz = (p.charge * AU2Coul * ((p.vx * magneticField.y)-(p.vy * magneticField.x))) + energy_lost.z;  
      p.elosstotaldist += p.dr; 
    } else {
      if (prev_inside_bubble) { 
        printf("\nParticle has left the bubble.\n"); 
        p.stepout = original_stepout; 
      }
      prev_inside_bubble = false; 
      magneticField = MagField(p, s); 
      vector energy_lost = synchrotron(p, magneticField); 
      p.Fx = (p.charge * AU2Coul * ((p.vy * magneticField.z)-(p.vz * magneticField.y))) + energy_lost.x * p.dr; 
      p.Fy = (p.charge * AU2Coul * ((p.vz * magneticField.x)-(p.vx * magneticField.z))) + energy_lost.y * p.dr;  
      p.Fz = (p.charge * AU2Coul * ((p.vx * magneticField.y)-(p.vy * magneticField.x))) + energy_lost.z * p.dr; 
      stepeloss = rmag(energy_lost.x*p.dr, energy_lost.y*p.dr, energy_lost.z*p.dr); 
      if(p.charge == -1 && p.mass < SubAtomicNumber2Kg){ 
        p.totaleloss += stepeloss; 
        p.elosstotaldist += p.dr; 
      }
    }
    if(p.step % p.stepout == 0){ 
      logBField(bfieldout, p, magneticField); 
    }
#else 
    stepeloss = 0.0; // 
    vector magneticField = MagField(p, s); // 
    p.Fx = (p.charge * AU2Coul * ((p.vy * magneticField.z)-(p.vz * magneticField.y))); // 
    p.Fy = (p.charge * AU2Coul * ((p.vz * magneticField.x)-(p.vx * magneticField.z))); // 
    p.Fz = (p.charge * AU2Coul * ((p.vx * magneticField.y)-(p.vy * magneticField.x))); // 
    logBField(bfieldout, p, magneticField); // 
#endif 
#else 
#ifdef ENABLE_ELOSS
    p.curr_energy = (gamma1(p)-1)*p.mass*SoL*SoL*J2eV; // 
    while(p.curr_energy < e.energy[e.currindex]){ // 
      e.currindex--; // 
    }
    loss_energy = interp_eloss(e, p.curr_energy, 0) * e.eloss_multiplier; // 
    vector energy_lost = calcEloss(p, g, e, loss_energy); // 
    stepeloss = rmag(energy_lost.x*p.dr, energy_lost.y*p.dr, energy_lost.z*p.dr); // 
    p.totaleloss += stepeloss; // 
    p.elosstotaldist += p.dr; // 
    p.Fx = energy_lost.x; // 
    // p.Fy = energy_lost.y; // Typo in original, assuming Fx, Fy, Fz intended
    // p.Fz = energy_lost.z; // Typo in original
    p.Fy = energy_lost.y; // Corrected from original Fx = energy_lost.y 
    p.Fz = energy_lost.z; // Corrected from original Fx = energy_lost.z 
#else 
    stepeloss = 0.0; // 
    p.Fx = 0; // 
    p.Fy = 0; // 
    p.Fz = 0; // 
#endif 
#endif 
    p = update_particle(p, e, s); 
    if (p.step % p.stepout == 0){ 
      logParticle(ptfileout, p, stepeloss); 
    }
  } 
  fclose(ptfileout);
  #ifdef ENABLE_BFIELD
  fclose(bfieldout); // Close bfieldout only if it was opened
  #endif
#endif
  done(p, s);
  printf("Finished with run %d\n",current_run_number);
}
