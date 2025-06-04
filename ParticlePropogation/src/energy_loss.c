#include <auxiliary.h>
#include <motion.h>
#include <energy_loss.h>

/*compute the shell correction for the bethe equation
  inputs: particle p, gas g
  output: shell correction coefficient
*/
// double shell_correction(particle p, gas g) {
//     double n = gamma1(p) * rmag(p.vx, p.vy, p.vz) / SoL;
//     return (0.422377 * pow(n, -2) + 0.0304043 * pow(n, -4) - 0.00038106 * pow(n, -6))*.000001 * g.excite_potential * g.excite_potential + (3.850190 * pow(n, -2) - 0.1667989 * pow(n, -4) + 0.00157955 * pow(M_PI, -6))*.000000001 * pow(g.excite_potential, 3);
// }

/*compute the density correction for the bethe equation
  inputs: particle p, gas g
  output: density correction coefficient
*/
// double density_correction(particle p, gas g) {
//     double x = log10(rmag(p.vx, p.vy, p.vz) * gamma1(p) / SoL);
//     double vp = pow((6.022e23 * g.density * g.atomic_number) / (g.atomic_mass), .5);
//     double C = -(2 * log((g.excite_potential) / (Planck * vp)) + 1);
//     if (p.x < 1.738) {
//         return 0;
//     } else if (p.x > 1.738 && p.x < 4.13) {
//         return 4.6052 * x + C + 0.1534 * pow((1.738 - x), 3.21);
//     } else {
//         return 4.6052 * x + C;
//     }
//
// }

/*compute the maximum angular frequency the particle can have in the gas?
  inputs: particle p, gas g
  output: maximum w for the particle

  need more information about this
*/
// double w_max(particle p, gas g) {
//   return (4 * p.mass * g.atomic_mass/Kg2amu) / pow(g.atomic_mass/Kg2amu + p.mass, 2);
// }

// double bethe(particle p, gas g) {
//     double coeff = 0.1535 * p.charge * p.charge / pow(betta(p),2) * (g.atomic_number / g.atomic_mass);
//     double bracket = log((2*m*DELTA*pow(betta(p),2)*pow(gamma1(p),2))/I) - 2 * pow(betta(p),2);
//     bracket += 27.631 - delta;
//
//
//     double result_1 = (g.bethe_constant * g.density * g.atomic_number * (g.charge_incident / g.atomic_mass)) / pow(rmag(p.vx, p.vy, p.vz), 2); /*all values outside brackets*/
//     double result_2 = (log((2 * p.mass * pow(gamma1(p), 2) * pow(rmag(p.vx, p.vy, p.vz), 2) * w_max(p, g)) / pow(g.excite_potential, 2))) - 2 \
//     * pow(rmag(p.vx, p.vy, p.vz), 2) / (SoL * SoL) - density_correction(p, g) - 2 * shell_correction(p, g) / g.atomic_number;
//     return result_2*result_1;
// }

/*compute the density correction for the bethe equation
  inputs: particle p, gas g
  output: density correction factor

  from: https://pdg.lbl.gov/2007/reviews/passagerpp.pdf

  constants are for nitrogen, taken from:
  https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir83-2785.pdf
*/
double density_correction(particle p, gas g){
  double x0 = 1.7378;
  double x1 = 4.1323;
  double density_C = -10.5400;
  double x = log10(betta(p) * gamma1(p));
  double a = 0.15349;
  double m = 3.2125;
  if (x >= x1) {
    return 2 * log(10) * x - density_C;
  } else if (x < x1 && x >= x0) {
    return 2 * log(10) * x - density_C + a * pow(x1 - x, m);
  } else {
    return 0; //non-conductive
  }
}

/*compute the second order Bloch (shell) correction to Bethe-Bloch equation
  inputs: particle p
  output: value of the Bloch correction term

  from ICRU Report 49 (eq. 2.6-2.7)
*/
double bloch_correction(particle p){
  double y = p.charge * FineStructure / betta(p);
  double correction = 1.042 - 0.8549 * pow(y,2) + 0.343 * pow(y,4);
  return -1 * pow(y,2) * (1.20206 - pow(y,2)*correction);
}

/*Compute maximum kinetic energy imparted to a free electron in a single collision
  inputs: particle p
  output: maximum kinetic energy for an incident electron

  from: https://pdg.lbl.gov/2007/reviews/passagerpp.pdf

  maybe include the difference between incident and scattered electron from:
  https://hallaweb.jlab.org/12GeV/experiment/E12-07-108/Weekly_meetings/12Oct2015/thir_12/Ioni_energy_loss.pdf
*/
double Tmax(particle p){
  double particleMass = p.mass * 931.49432; // needs to be in MeV
  double numerator = 2 * 0.510998918 * pow(betta(p),2) * pow(gamma1(p),2);
  double denominator = 1 + 2 * 0.510998918 * gamma1(p) / particleMass + pow(0.510998918 / particleMass, 2);
  return numerator / denominator;
}

/*compute the bethe-bloch equation value for energy loss in a gas
  inputs: particle p, gas g
  output: value of the bethe-bloch equation for a given velocity

  This is only called for electrons

  from: https://pdg.lbl.gov/2007/reviews/passagerpp.pdf
*/
double bethe(particle p, gas g){
  double atomicMass = g.atomic_mass; // needs to be in g/mol
  double KA_ratio = 0.307075 / (atomicMass); // MeV cm^2/g
  double dedx = KA_ratio * pow(p.charge,2) * g.atomic_number / pow(betta(p),2);
  double bracket = 0.5 * log(2 * 0.510998918 * pow(betta(p),2) * pow(gamma1(p),2) * Tmax(p) / pow(g.excite_potential,2));
  bracket -= pow(betta(p),2) - density_correction(p, g)/2 + bloch_correction(p);
  bracket *= dedx; // MeV cm^2 / g
  return bracket * g.density * 1.60218e-13 / 10000; // Joules / m
}

/*compute beta (v/c) for the particle
  inputs: paricle p
  output: beta value for the particle

  maybe move this to motion.h/c?
*/
double betta(particle p) {
    return rmag(p.vx, p.vy, p.vz) / SoL;
}

/*interpolate the energy loss value from the values in the SRIM file
  inputs: energy loss arrays e, current energy of particle currE
  output: interpolated energy lost for a given energy
*/
double interp_eloss(eloss e, double currE, int selector){
    double differnce, low_weight, low_diff, high_weight, high_diff;
    // printf("next energy: %12.6e\n current energy: %12.6e\n",e.energy[e.currindex+1],e.energy[e.currindex]);
    differnce = e.energy[e.currindex+1]-e.energy[e.currindex];
    low_diff = currE - e.energy[e.currindex];
    high_diff = e.energy[e.currindex+1] - currE;
    high_weight = low_diff/differnce;
    low_weight = high_diff/differnce;
        if (selector == 1) { 
          //electronic eloss
        return low_weight*e.electronic[e.currindex]+high_weight*e.electronic[e.currindex+1];
    } else if (selector == 2) { 
        //nuclear eloss
        return low_weight*e.nuclear[e.currindex]+high_weight*e.nuclear[e.currindex+1];
    }
    //force from both components
    return low_weight*e.total[e.currindex]+high_weight*e.total[e.currindex+1];
}

/*calculate the energy lost due to synchrotron radiation in x,y,z direction
  inputs: particle p, vector of the magnetic field Bfield
  output: vector of energy lost in each direction
*/
vector synchrotron(particle p, vector Bfield){
  // printf("synchrotron start\n");
  vector beta, beta_dot;
  beta.x = p.vx / SoL; beta.y = p.vy / SoL; beta.z = p.vz / SoL;
  beta_dot.x = p.ax / SoL; beta_dot.y = p.ay / SoL; beta_dot.z = p.az / SoL;
  double coeff = p.charge * p.charge * AU2Coul * AU2Coul;
  coeff = coeff / (6 * M_PI * Epsilon0 * SoL);
  coeff = coeff * pow(gamma1(p), 6);
  // printf("coeff synch %12.6e\n",coeff);
  double energy_lost = dot(beta_dot, beta_dot) - dot(cross(beta, beta_dot), cross(beta, beta_dot));
  // printf("synch before coeff %12.6e\n",energy_lost);
  energy_lost = -1 * energy_lost * coeff;
  // printf("synchrotron %12.6e\n",energy_lost);
  vector returnVec;
  returnVec.x = energy_lost * p.vx / rmag(p.vx, p.vy, p.vz);
  returnVec.y = energy_lost * p.vy / rmag(p.vx, p.vy, p.vz);
  returnVec.z = energy_lost * p.vz / rmag(p.vx, p.vy, p.vz);
  // printf("synchrotron end\n");
  return returnVec;
}

/*calculate the energy lost at a particular time step in the x,y,z direction
  inputs: particle p, gas g, eloss e, magnetic field vector Bfield, energy lost (from interp eloss)
  output: vector of the energy lost in each direction

  maybe move interp_eloss call inside this function?
*/
vector calcEloss(particle p, gas g, eloss e, double lost_energy) {
  // lost_energy inputted as eV from interp_eloss function
  double eloss_x, eloss_y, eloss_z, v_mag;
  v_mag = rmag(p.vx,p.vy,p.vz);
  //calculate energy loss for an electron using bethe-bloch equation for beta > 0.15
  if (p.charge == -1 && p.mass < 1 && betta(p) > 0.15) {
    double bethe_step = bethe(p, g); //
    // printf("bethe_eloss: %12.6e\n", bethe_step);
    eloss_x = bethe_step * p.vx / v_mag;
    eloss_y = bethe_step * p.vy / v_mag;
    eloss_z = bethe_step * p.vz / v_mag;
    // printf("eloss: %12.6e, %12.6e, %12.6e\n", eloss_x,eloss_y,eloss_z);
  } else {
    eloss_x = -(1.602176e-19)*lost_energy*p.vx/v_mag;
    eloss_y = -(1.602176e-19)*lost_energy*p.vy/v_mag;
    eloss_z = -(1.602176e-19)*lost_energy*p.vz/v_mag; //converted components from eV; to joules
    // eloss_x = -lost_energy*p.vx/v_mag;
    // eloss_y = -lost_energy*p.vy/v_mag;
    // eloss_z = -lost_energy*p.vz/v_mag;
  }
  vector returnVec;
  returnVec.x = eloss_x;
  returnVec.y = eloss_y;
  returnVec.z = eloss_z;
  // returnVec.x = 0; // testing no eloss
  // returnVec.y = 0;
  // returnVec.z = 0;

  return returnVec;
}
