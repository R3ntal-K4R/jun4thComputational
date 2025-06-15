#ifndef GENERATOR_ORIGIN_HIGH_ENERGY_H
#define GENERATOR_ORIGIN_HIGH_ENERGY_H

typedef struct {
    double x, y, z;          // Position coordinates
    double direction_x, direction_y, direction_z; // Direction vector
    double energy;           // Energy of the particle
    double velocity_x, velocity_y, velocity_z; // Velocity vector
} Particle;

int main_origin_high_energy();

#endif
