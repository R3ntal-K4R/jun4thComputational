#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>

double* box_muller(double* rand_theta, double* rand_phi)
{
	double u, v;
	double min, max, range, div;

	assert(rand_theta);
	assert(rand_phi);

	min = 0.0;
	max = 1.0;

	range = max - min;
	div = RAND_MAX / range;

	u = min + (rand() / div);
	v = min + (rand() / div);

	//printf("%lf    %lf\n", u, v);

	*rand_theta = M_PI * sqrt(-2 * log(u)) * cos(2 * M_PI * v);
	*rand_phi = 2 * M_PI * sqrt(-2 * log(u)) * sin(2 * M_PI * v);
}

double* vector_rotation(double theta, double phi, double vx, double vy, double vz, double* vx_prime, double* vy_prime, double* vz_prime)
{
	double v_mag, vx_dir, vy_dir, vz_dir;
	double R[3][3], I[3][3], matrix[3][3], K[3][3], K_squared[3][3], R_phi[3][3];
	double cross[3], v_matrix[3], v_matrix_prime[3];
	double cross_mag, rho, phi_prime, straggling_dist, eta, rotate;

	assert(vx_prime);
	assert(vy_prime);
	assert(vz_prime);

	v_mag = sqrt(vx * vx + vy * vy + vz * vz);
	vx_dir = vx / v_mag;
	vy_dir = vy / v_mag;
	vz_dir = vz / v_mag;

	//put zeroes in matrices
	for(int i = 0; i <= 2; i++)
	{
		for(int j = 0; j <= 2; j++)
		{
			K[i][j] = 0.0;
			R[i][j] = 0.0;
			R_phi[i][j] = 0.0;
			I[i][j] = 0.0;
			K_squared[i][j] = 0.0;
		}
	}

	rotate = atan(vy/vx);
	R_phi[0][0] = cos(rotate);
	R_phi[0][1] = -sin(rotate);
	R_phi[1][0] = sin(rotate);
	R_phi[1][1] = cos(rotate);
	R_phi[2][2] = 1.0;

	//cross product of velocity vector and z-aligned vector
	cross[0] = -vy * v_mag;
	cross[1] = -v_mag * vx;
	cross[2] = 0.0;
	cross_mag = sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

	//assigning k matrix of rodrigues formula from the dot product
	K[0][2] = cross[0] / cross_mag;
	K[1][2] = cross[1] / cross_mag;
	K[2][0] = -cross[0] / cross_mag;
	K[2][1] = -cross[1] / cross_mag;

	//make sure not dividing by 0
	if(cross[0] == 0.0 && cross[1] == 0.0)
	{
		K[0][2] = 0.0;
		K[1][2] = 0.0;
		K[2][0] = 0.0;
		K[2][1] = 0.0;
	}

	//dot product equals eta
	eta = (vz * v_mag) / (v_mag * v_mag);

	//applying rotation of theta and phi
	v_matrix[0] = sin(theta) * cos(phi);
	v_matrix[1] = sin(theta) * sin(phi);
	v_matrix[2] = cos(theta);

	//declaring identity matrix
	I[0][0] = 1.0;
	I[1][1] = 1.0;
	I[2][2] = 1.0;

	//calculate K_squared
	for(int i = 0; i <= 2; i++)
	{
		for(int j = 0; j <= 2; j++)
		{
			for(int k = 0; k <= 2; k++)
			{
				K_squared[i][j] += K[i][k] * K[k][j];
			}
		}
	}

	for(int i = 0; i <= 2; i++)
	{
		for(int j = 0; j <= 2; j++)
		{
			R[i][j] = I[i][j] + sqrt(1 - eta * eta) * K[i][j] + K_squared[i][j] * (1 - eta);
		}
	}

	for(int i = 0; i <= 2; i++)
	{
		for(int j = 0; j <= 2; j++)
		{
			v_matrix_prime[i] += R_phi[j][i] * R[j][i] * v_matrix[j] * v_mag;
		}
	}

	*vx_prime = v_matrix_prime[0];
	*vy_prime = v_matrix_prime[1];
	*vz_prime = v_matrix_prime[2];

}

void main()
{
	FILE *test_angular_straggling = fopen("test_angular_straggling.dat", "w");
	int i;
	double rand_theta, rand_phi, pi;
	double r, theta, phi;
	double x0, y0, z0, vx0, vy0, vz0;
	double x, y, z, vx, vy, vz, vx_prime, vy_prime, vz_prime, v_mag;
	double t, i_max, t_final, dt;
	srand(time(NULL));

	pi = acos(-1);

	//initial conditions in m and m/s
	t = 0.0;
	x0 = 1.0;
	y0 = 1.0;
	z0 = 1.0;
	vx0 = 1.0;
	vy0 = 1.0;
	vz0 = 1.0;
	dt = 1.0;
	t_final = 4;

	x = x0;
	y = y0;
	z = z0;
	vx = vx0;
	vy = vy0;
	vz = vz0;
	v_mag = sqrt(vx * vx + vy * vy + vz * vz);

	fprintf(test_angular_straggling, "# t, x, y, z, vx, vy, vz, v_mag, rand_theta, rand_phi\n");
	fprintf(test_angular_straggling, "%lf %lf %lf %lf %lf %lf %lf %lf\n", t, x, y, z, vx, vy, vz, v_mag);

	i_max = t_final / dt;

	for(i = 0; i <= i_max; i++){
		//get random angles
		box_muller(&rand_theta, &rand_phi);
		rand_theta = rand_theta * pi;
		rand_phi = rand_phi * pi;

		/*
		test angles
		rand_theta = 0.0;
		rand_phi = M_PI/2;
		*/

		vector_rotation(rand_theta, rand_phi, vx, vy, vz, &vx_prime, &vy_prime, &vz_prime);

		vx = vx_prime;
		vy = vy_prime;
		vz = vz_prime;

		t = t + dt;

		x = x + vx * dt;
		y = y + vy * dt;
		z = z + vz * dt;

		v_mag = sqrt(vx * vx + vy * vy + vz * vz);

		fprintf(test_angular_straggling, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t, x, y, z, vx, vy, vz, v_mag, rand_theta, rand_phi);
	}

	fclose(test_angular_straggling);
}