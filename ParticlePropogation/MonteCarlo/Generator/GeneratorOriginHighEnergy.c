#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "GeneratorOriginHighEnergy.h"

// Define the boundaries of the space
double x_min = -5.0, x_max = 5.0;
double y_min = -5.0, y_max = 5.0;
double z_min = -5.0, z_max = 5.0;

// Define the target point
double target_x = 25.0;
double target_y = 25.0;
double target_z = 35.0;

// Function to generate a random number between 0 and 1
double rand_double() {
    return rand() / (double)RAND_MAX;
}

// Function to generate a random position on the border
void generate_position_border(Particle *particle) {
    // Randomly choose a face of the cube
    int face = rand() % 6;

    switch (face) {
        case 0: // x = x_min
            particle->x = x_min;
            particle->y = y_min + rand_double() * (y_max - y_min);
            particle->z = z_min + rand_double() * (z_max - z_min);
            break;
        case 1: // x = x_max
            particle->x = x_max;
            particle->y = y_min + rand_double() * (y_max - y_min);
            particle->z = z_min + rand_double() * (z_max - z_min);
            break;
        case 2: // y = y_min
            particle->x = x_min + rand_double() * (x_max - x_min);
            particle->y = y_min;
            particle->z = z_min + rand_double() * (z_max - z_min);
            break;
        case 3: // y = y_max
            particle->x = x_min + rand_double() * (x_max - x_min);
            particle->y = y_max;
            particle->z = z_min + rand_double() * (z_max - z_min);
            break;
        case 4: // z = z_min
            particle->x = x_min + rand_double() * (x_max - x_min);
            particle->y = y_min + rand_double() * (y_max - y_min);
            particle->z = z_min;
            break;
        case 5: // z = z_max
            particle->x = x_min + rand_double() * (x_max - x_min);
            particle->y = y_min + rand_double() * (y_max - y_min);
            particle->z = z_max;
            break;
    }
}

// Function to generate a direction vector towards the target
void generate_direction_to_target(Particle *particle) {
    // Calculate the direction vector towards the target (target_x, target_y, target_z)
    double direction_x = target_x - particle->x;
    double direction_y = target_y - particle->y;
    double direction_z = target_z - particle->z;

    // Normalize the direction vector
    double magnitude = sqrt(direction_x * direction_x + direction_y * direction_y + direction_z * direction_z);
    particle->direction_x = direction_x / magnitude;
    particle->direction_y = direction_y / magnitude;
    particle->direction_z = direction_z / magnitude;
}

// Function to generate a high energy value and calculate velocity
void generate_energy_and_velocity(Particle *particle) {
    // Define the energy range
    double energy_min = 10000.0;  // Minimum energy (eV)
    double energy_max = 20000.0; // Maximum energy (eV)

    // Generate a random energy value within the defined range
    particle->energy = energy_min + rand_double() * (energy_max - energy_min);

    // Assuming the particle is an electron
    double mass_electron = 510998.9461; // Electron mass in eV/c^2
    double speed_of_light = 299792458;  // Speed of light in m/s

    // Calculate the kinetic energy
    double kinetic_energy = particle->energy;

    // Calculate the velocity magnitude
    double velocity_magnitude = speed_of_light * sqrt(2 * kinetic_energy / mass_electron);

    // Calculate the velocity components
    particle->velocity_x = velocity_magnitude * particle->direction_x;
    particle->velocity_y = velocity_magnitude * particle->direction_y;
    particle->velocity_z = velocity_magnitude * particle->direction_z;
}

// Function to generate a particle with properties: border position, direction to center, high energy
Particle generate_particle_border_to_target() {
    Particle particle;

    generate_position_border(&particle);
    generate_direction_to_target(&particle);
    generate_energy_and_velocity(&particle);

    return particle;
}

// Function to write particle data to a file in the specified format
void write_particle_data(Particle particle, FILE *fp) {
    fprintf(fp, "%lf %lf %lf\n", particle.x, particle.y, particle.z);
    fprintf(fp, "%lf %lf %lf\n", particle.velocity_x, particle.velocity_y, particle.velocity_z);
}

int main_origin_high_energy() {
    // Seed the random number generator
    srand(time(NULL));

    int num_particles = 250; // Number of particles to generate

    // Open a file to write the particle data
    FILE *fp = fopen("ParticlePropogation/MonteCarlo/Generator/input.dat", "w");
    if (fp == NULL) {
        perror("Error opening file");
        return 1; // Indicate an error occurred
    }

    // Write the number of particles
    fprintf(fp, "%d\n", num_particles);

    // Generate and write data for each particle
    for (int i = 0; i < num_particles; i++) {
        Particle particle = generate_particle_border_to_target();
        write_particle_data(particle, fp);
    }

    // Write the other parameters
    fprintf(fp, "26    5.585e+01\n");
    fprintf(fp, "1.0e-11 	6.557e+15\n");
    fprintf(fp, "8000\n");
    fprintf(fp, "5\n");
    fprintf(fp, "1.000\n");
    fprintf(fp, "1\n");
    fprintf(fp, "0\n");

    // Close the file
    fclose(fp);

    printf("Generated %d particles from the border angled to (25, 25, 35) with high energy and wrote data to input.dat\n", num_particles);

    return 0; // Indicate successful execution
}
