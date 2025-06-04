#include <auxiliary.h>
#include <motion.h>

/*calculates the gamma of the particle
  inputs: particle p
  outputs: gamma of particle
*/
double gamma1(particle p) {
  double v_magnitude2 = ((p.vx * p.vx)+(p.vy * p.vy)+(p.vz * p.vz));
  double c2 = SoL*SoL;
  double Beta2 = v_magnitude2 / (c2);
  double g = 1.0e0 / pow((1.0e0 - Beta2), 0.5);
  return g;
}

/*add inputs in quadrature
  inputs: three doubles to be added in quadrature
  outputs: double
*/
double rmag(double x, double y, double z) {
    return pow(x * x + y * y + z * z, .5);
}

/*calculates acceleration of particle in x direction from x-boost
  inputs: particle p
  output: Axx value for particle
*/
double Axx(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*((p.vx) * p.vx)+(gamma1(p) * p.mass);
}

/*calculates acceleration of particle in x direction from y-boost
  inputs: particle p
  output: Axy value for particle
*/
double Axy(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vx * p.vy);
}

/*calculates acceleration of particle in x direction from z-boost
  inputs: particle p
  output: Axz value for particle
*/
double Axz(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vx * p.vz);
}

/*calculates acceleration of particle in y direction from x-boost
  inputs: particle p
  output: Ayx value for particle
*/
double Ayx(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vy * p.vx);
}

/*calculates acceleration of particle in y direction from y-boost
  inputs: particle p
  output: Ayy value for particle
*/
double Ayy(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vy * p.vy)+(gamma1(p) * p.mass);
}

/*calculates acceleration of particle in y direction from z-boost
  inputs: particle p
  output: Ayz value for particle
*/
double Ayz(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vy * p.vz);
}

/*calculates acceleration of particle in z direction from x-boost
  inputs: particle p
  output: Azx value for particle
*/
double Azx(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vz * p.vx);
}

/*calculates acceleration of particle in z direction from y-boost
  inputs: particle p
  output: Azy value for particle
*/
double Azy(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vz * p.vy);
}

/*calculates acceleration of particle in z direction from z-boost
  inputs: particle p
  output: Azz value for particle
*/
double Azz(particle p) {
    return ((pow(gamma1(p), 3.0) * p.mass) / (SoL * SoL))*(p.vz * p.vz)+(gamma1(p) * p.mass);
}

/*calculate the acceleration of the particle from forces and relativistic effects
  inputs: particle p
  outputs: x,y,z component of acceleration
*/
vector calcAcceleration(particle p){
  // denom may be giving a nan error
  // printf("part v inserted to acceleration calc: %10e\n%10e %10e %10e\n",rmag(p.vx,p.vy,p.vz),p.vx,p.vy,p.vz);
  double ax,ay,az;
  double denom = Axx(p) * Ayy(p) * Azz(p) - Axx(p) * Ayz(p) * Azy(p) - Axy(p) * Ayx(p) * Azz(p) + Axy(p) * Ayz(p) * Azx(p) + Axz(p) * Ayx(p) * Azy(p) - Axz(p) * Ayy(p) * Azx(p);
  if (denom != 0) { // DBL_EPSILON = 2.220446e-16
    ax = (p.Fx * Ayy(p) * Azz(p) - p.Fx * Ayz(p) * Azy(p) - p.Fy * Axy(p) * Azz(p) + p.Fy * Axz(p) * Azy(p) + p.Fz * Axy(p) * Ayz(p) - p.Fz * Axz(p) * Ayy(p))
    / (Axx(p) * Ayy(p) * Azz(p) - Axx(p) * Ayz(p) * Azy(p) - Axy(p) * Ayx(p) * Azz(p) + Axy(p) * Ayz(p) * Azx(p) + Axz(p) * Ayx(p) * Azy(p) - Axz(p) * Ayy(p) * Azx(p));
    ay = (p.Fy * Axx(p) * Azz(p) - p.Fz * Axx(p) * Ayz(p) - p.Fx * Ayx(p) * Azz(p) + p.Fx * Ayz(p) * Azx(p) + p.Fz * Ayx(p) * Axz(p) - p.Fy * Axz(p) * Azx(p))
    / (Axx(p) * Ayy(p) * Azz(p) - Axx(p) * Ayz(p) * Azy(p) - Axy(p) * Ayx(p) * Azz(p) + Axy(p) * Ayz(p) * Azx(p) + Axz(p) * Ayx(p) * Azy(p) - Axz(p) * Ayy(p) * Azx(p));
    az = (p.Fx * Ayx(p) * Azy(p) - p.Fx * Ayy(p) * Azx(p) - p.Fy * Axx(p) * Azy(p) + p.Fy * Axy(p) * Azx(p) + p.Fz * Axx(p) * Ayy(p) - p.Fz * Axy(p) * Ayx(p))
    / (Axx(p) * Ayy(p) * Azz(p) - Axx(p) * Ayz(p) * Azy(p) - Axy(p) * Ayx(p) * Azz(p) + Axy(p) * Ayz(p) * Azx(p) + Axz(p) * Ayx(p) * Azy(p) - Axz(p) * Ayy(p) * Azx(p));
  } else {
    ax = 0.0;
    ay = 0.0;
    az = 0.0;
  }
  vector returnVec;
  returnVec.x = ax;
  returnVec.y = ay;
  returnVec.z = az;
  // printf("Calc Accel return %10e \n %10e %10e %10e \n",rmag(ax,ay,az),ax,ay,az);
  return returnVec;
}

/*updates the position, velocity, and acceleration of particle
  inputs: particle p at previous step
  outputs: update particle struct (technically a new struct)
*/
particle update_particle(particle p, eloss e, spacecraft s){
  // printf("update_particle start\n");
  // printf("Particle v before updating: %10e\n%10e %10e %10e\n",rmag(p.vx,p.vy,p.vz),p.vx,p.vy,p.vz);
  double magnitude = rmag(p.vx, p.vy, p.vz);
  double new_x,new_y,new_z; //for updating dr
  bool status = InBubble(p,s); // called once, used twice
  if (magnitude/SoL <= 1e-4){
    printf("\nVelocity of particle too low\n");
    // done(p,s);
    p.run_status = false; // end current run
    // printf("run status: %d\n", p.run_status);
  }
  if (p.straggle_switch == 1 && status){
    // printf("Straggled on step: %4i\n", p.step);
    vector vel, straggled_vel, prev_vel;
    vel.x = p.vx; vel.y = p.vy; vel.z = p.vz;
    prev_vel.x = p.prev_vx; prev_vel.y = p.prev_vy; prev_vel.z = p.prev_vz;
    straggled_vel = straggle_velocity(vel, prev_vel, e);
    p.vx = straggled_vel.x;
    p.vy = straggled_vel.y;
    p.vz = straggled_vel.z;
    if (fabs(magnitude - rmag(p.vx, p.vy, p.vz)) > 1e-7){
      printf("The velocity has changed by %6.3e during straggling\n", fabs(magnitude - rmag(p.vx, p.vy, p.vz)));
    }
    // exit(1);
  }
  p.prev_vel = p.curr_vel;
  p.prev_vx = p.vx;
  p.prev_vy = p.vy;
  p.prev_vz = p.vz;

  p.totaldist += p.dr;
  vector curr_a = calcAcceleration(p);
  // printf("curr_a about to be inserted %10e \n %10e %10e %10e \n",rmag(curr_a.x,curr_a.y,curr_a.z),curr_a.x,curr_a.y,curr_a.z); //a is wrong rn 7/25/2023
  p.x = p.x + p.vx * p.dt + 0.5 * curr_a.x * p.dt * p.dt;
  p.y = p.y + p.vy * p.dt + 0.5 * curr_a.y * p.dt * p.dt;
  p.z = p.z + p.vz * p.dt + 0.5 * curr_a.z * p.dt * p.dt;
  p.vx = p.vx + curr_a.x * p.dt;
  p.vy = p.vy + curr_a.y * p.dt;
  p.vz = p.vz + curr_a.z * p.dt;
  p.curr_vel = rmag(p.vx, p.vy, p.vz);
  p.ax = curr_a.x;
  p.ay = curr_a.y;
  p.az = curr_a.z;
  p.step++;
  // printf("update_particle finish\n");
  // printf("Particle v after updating: %10e\n%10e %10e %10e\n",rmag(p.vx,p.vy,p.vz),p.vx,p.vy,p.vz);
  return p;
}

/*straggle particle velocity using transformation to and from Fresnel system
  inputs: vector of current velocity, vector of previous step velocity, eloss arrays
  output: vector of new straggles velocity

*/
vector straggle_velocity(vector v, vector prev_v, eloss e){

  vector gaussian = box_muller(e.lateral[e.currindex],e.longitudinal[e.currindex]);
  // printf("lat, long: %12.6e %12.6e\n",e.lateral[e.currindex],e.longitudinal[e.currindex]);
  // printf("gaussian: %12.6e %12.6e %12.6e\n",gaussian.x,gaussian.y,gaussian.z);
  double phi = gaussian.z * 2 * M_PI;
  double theta = atan2(gaussian.x, e.range[e.currindex] + gaussian.y);

  double mag = rmag(v.x, v.y, v.z);

  vector local_velocity, lab_velocity;
  local_velocity.x = mag * sin(theta) * cos(phi);   
  local_velocity.y = mag * sin(theta) * sin(phi);
  local_velocity.z = mag * cos(theta);

  // printf("phi: %12.6e theta: %12.6e\n",phi, theta);
  // printf("Magnitude: %8.3e\n", mag);
  printf("%8.3e\n %12.6e %12.6e %12.6e\n", rmag(local_velocity.x, local_velocity.y, local_velocity.z),local_velocity.x, local_velocity.y, local_velocity.z);
  // printf("%f %f %f", local_velocity.x, local_velocity.y, local_velocity.z);

  fresnel_vectors fresnel_vec = fresnel(v, prev_v);
  // printf("Fresnel Magnitudes:\n");
  // printf("Tangent: %8.3e\n", rmag(fresnel_vec.Tangent.x, fresnel_vec.Tangent.y, fresnel_vec.Tangent.z));
  // printf("%8.3e %8.3e %8.3e\n", fresnel_vec.Tangent.x, fresnel_vec.Tangent.y, fresnel_vec.Tangent.z);
  // printf("Normal: %8.3e\n", rmag(fresnel_vec.Normal.x, fresnel_vec.Normal.y, fresnel_vec.Normal.z));
  // printf("%8.3e %8.3e %8.3e\n", fresnel_vec.Normal.x, fresnel_vec.Normal.y, fresnel_vec.Normal.z);
  // printf("Binormal: %8.3e\n", rmag(fresnel_vec.Binormal.x, fresnel_vec.Binormal.y, fresnel_vec.Binormal.z));
  // printf("%8.3e %8.3e %8.3e\n", fresnel_vec.Binormal.x, fresnel_vec.Binormal.y, fresnel_vec.Binormal.z);

  
  lab_velocity.x = fresnel_vec.Tangent.x * local_velocity.x + fresnel_vec.Normal.x * local_velocity.y + fresnel_vec.Binormal.x * local_velocity.z;
  lab_velocity.y = fresnel_vec.Tangent.y * local_velocity.x + fresnel_vec.Normal.y * local_velocity.y + fresnel_vec.Binormal.y * local_velocity.z;
  lab_velocity.z = fresnel_vec.Tangent.z * local_velocity.x + fresnel_vec.Normal.z * local_velocity.y + fresnel_vec.Binormal.z * local_velocity.z;
  // lab_velocity.x = fresnel_vec.Tangent.x * local_velocity.x + fresnel_vec.Tangent.y * local_velocity.y + fresnel_vec.Tangent.z * local_velocity.z;
  // lab_velocity.y = fresnel_vec.Normal.x * local_velocity.x + fresnel_vec.Normal.y * local_velocity.y + fresnel_vec.Normal.z * local_velocity.z;
  // lab_velocity.z = fresnel_vec.Binormal.x * local_velocity.x + fresnel_vec.Binormal.y * local_velocity.y + fresnel_vec.Binormal.z * local_velocity.z;

  // printf("Lab velocity: %8.3e\n %12.6e %12.6e %12.6e\n", rmag(lab_velocity.x, lab_velocity.y, lab_velocity.z), lab_velocity.x, lab_velocity.y, lab_velocity.z);
  return lab_velocity;
}

/*transform to the Fresnel coordinate system
  inputs: vector of current velocity, vector of previous step velocity
  output: velocity vector in local Fresnel system
*/
fresnel_vectors fresnel(vector v, vector prev_v){

  vector Tangent, Tangent_prev, Normal, Binormal;
  double mag = rmag(v.x, v.y,v.z);
  // printf("MAG: %8.3e\n", mag);
  if (mag >= ZeroTolerance) {
    Tangent.x = v.x / mag;
    Tangent.y = v.y / mag;
    Tangent.z = v.z / mag;
  } else {
    Tangent.x = 0.0;
    Tangent.y = 0.0;
    Tangent.z = 0.0;
  }

  mag = rmag(prev_v.x, prev_v.y, prev_v.z);
  // printf("MAG: %8.3e\n", mag);
  if (mag >= ZeroTolerance) {
    Tangent_prev.x = prev_v.x / mag;
    Tangent_prev.y = prev_v.y / mag;
    Tangent_prev.z = prev_v.z / mag;
  } else {
    Tangent_prev.x = 0.0;
    Tangent_prev.y = 0.0;
    Tangent_prev.z = 0.0;
  }

  mag = rmag(Tangent.x - Tangent_prev.x, Tangent.y - Tangent_prev.y, Tangent.z - Tangent_prev.z);
  // printf("MAG: %16.14e\n", mag);
  if (mag >= ZeroTolerance) {
    Normal.x = (Tangent.x - Tangent_prev.x)/mag;
    Normal.y = (Tangent.y - Tangent_prev.y)/mag;
    Normal.z = (Tangent.z - Tangent_prev.z)/mag;
  } else {
    Normal.x = 0.0;
    Normal.y = 0.0;
    Normal.z = 0.0;
    // printf("%8.3e %8.3e %8.3e\n", Normal.x, Normal.y, Normal.z);
  }

  Binormal = cross(Tangent, Normal);

  // mag = rmag(Tangent.x, Tangent.y, Tangent.z);
  // printf("TANGENT MAG: %8.3e\n", mag);
  // mag = rmag(Normal.x, Normal.y, Normal.z);
  // printf("NORMAL MAG: %8.3e\n", mag);
  // mag = rmag(Binormal.x, Binormal.y, Binormal.z);
  // printf("BINORMAL MAG: %8.3e\n", mag);

  fresnel_vectors returnVec;
  returnVec.Tangent = Tangent;
  returnVec.Normal = Normal;
  returnVec.Binormal = Binormal;

  return returnVec;
}

/*straggles the current velocity of the particle
  inputs: velocity vector v
  outputs: straggled velocity vector
*/
// vector straggle_velocity(vector v){
//   vector z_vec, cross_vec;
//   double R_mat[3][3], I_mat[3][3], K_mat[3][3], K2_mat[3][3], R_phi[3][3];
//   double vel[3], vel_prime[3];
//   z_vec.x = 0; z_vec.y = 0;
//   z_vec.z = rmag(v.x,v.y,v.z);
//
//   for(int i = 0; i <= 2; i++){
// 		for(int j = 0; j <= 2; j++){
// 			K_mat[i][j] = 0.0;
// 			R_mat[i][j] = 0.0;
// 			R_phi[i][j] = 0.0;
// 			I_mat[i][j] = 0.0;
//       I_mat[i][i] = 1.0;
// 			K2_mat[i][j] = 0.0;
//       }
// 		}
//
//   double angle = atan(v.y/v.x);
//   R_phi[0][0] = cos(angle);
// 	R_phi[0][1] = -sin(angle);
// 	R_phi[1][0] = sin(angle);
// 	R_phi[1][1] = cos(angle);
// 	R_phi[2][2] = 1.0;
//
//   cross_vec = cross(z_vec, v);
//   double cross_mag = rmag(cross_vec.x, cross_vec.y, cross_vec.z);
//   K_mat[0][2] = cross_vec.x / cross_mag;
//   K_mat[1][2] = cross_vec.y / cross_mag;
//   K_mat[2][0] = -cross_vec.x / cross_mag;
//   K_mat[2][1] = -cross_vec.y / cross_mag;
//
//   vector angles = box_muller();
//   double theta = angles.x; double phi = angles.y;
//   vel[0] = sin(theta)*cos(phi);
//   vel[1] = sin(theta)*sin(phi);
//   vel[2] = cos(theta);
//
//   //calculate K_squared
// 	for(int i = 0; i <= 2; i++){
// 		for(int j = 0; j <= 2; j++){
// 			for(int k = 0; k <= 2; k++){
// 				K2_mat[i][j] += K_mat[i][k] * K_mat[k][j];
// 			}
// 		}
// 	}
//
//   double eta = v.z / rmag(v.x, v.y, v.z);
//   for(int i = 0; i <= 2; i++){
// 		for(int j = 0; j <= 2; j++){
// 			R_mat[i][j] = I_mat[i][j] + sqrt(1 - eta * eta) * K_mat[i][j] + K2_mat[i][j] * (1 - eta);
// 		}
// 	}
//
// 	for(int i = 0; i <= 2; i++){
// 		for(int j = 0; j <= 2; j++){
// 			vel_prime[i] += R_phi[j][i] * R_mat[j][i] * vel[j] * z_vec.z;
// 		}
// 	}
//
//   double vel_mag = rmag(vel_prime[0], vel_prime[1], vel_prime[2]);
//
//   z_vec.x = vel_prime[0] / vel_mag;
//   z_vec.y = vel_prime[1] / vel_mag;
//   z_vec.z = vel_prime[2] / vel_mag;
//   return z_vec;
// }
