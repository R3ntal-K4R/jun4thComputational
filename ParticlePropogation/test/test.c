  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <string.h>
  #include <assert.h>
  #include <stdbool.h>
  #include <float.h>

typedef struct vector {
  double x, y, z;
} vector;

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
  returnVec.z = sqrt(pow(2,(sqrt(-2 * log(w)) * sin(2 * M_PI * x) + mu))+pow(2,(sqrt(-2 * log(u)) * cos(2 * M_PI * x) + mu)));
  // printf("%12.6e %12.6e %12.6e \n", returnVec.z, w, x);
  return returnVec;
}


int main() {
    int i;
    double j;
    double div = RAND_MAX;
    FILE *test_BM = fopen("test_BM.dat", "w");
    fprintf(test_BM,"run    muller x     muller y      phi     theta             gaussian.z \n");
    // srand(time(0));
    printf("%12.6e\n", div);
    for (i = 0; i < 100000; i++) {
        j = rand() / div;
        // Code to be executed in each iteration
        vector gaussian = box_muller(1,1); //testing with 1 meter lat and long
        double phi = gaussian.z * 2 * M_PI;
        double theta = atan2(gaussian.x, 10.0 + gaussian.y); //range of 10 m
        fprintf(test_BM,"%d %12.6e %12.6e %12.6e %12.6e %12.6e\n",i,gaussian.x, gaussian.y, phi, theta, gaussian.z);
    }

    return 0;
}