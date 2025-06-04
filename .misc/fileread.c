#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

int randint(int n) {
  if ((n - 1) == RAND_MAX) {
    return rand();
  } else {
    // Supporting larger values for n would requires an even more
    // elaborate implementation that combines multiple calls to rand()
    assert (n <= RAND_MAX);

    // Chop off all of the values that would cause skew...
    int end = RAND_MAX / n; // truncate skew
    assert (end > 0);
    end *= n;

    // ... and ignore results from rand() that fall above that limit.
    // (Worst case the loop condition should succeed 50% of the time,
    // so we can expect to bail out of this loop pretty quickly.)
    int r;
    while ((r = rand()) >= end);

    return r % n;
  }}


int main () {
   FILE *fp;
   char c;
   int pos=5;
   double val;
   int value;
   int index;

   fp = fopen("Bfield.dat","r");
   int i;
   double time;
for(i =0; i<45; i++){
    index = randint(108000000);

    clock_t tic = clock();
    fseek(fp, i , SEEK_SET);

    clock_t toc = clock();
    fscanf(fp,"%lf",&val);

    printf("num at position %d is %lf\n",i,val);
    double runtime = (double)(toc - tic) / CLOCKS_PER_SEC;
    time += runtime;
    printf("runtime = %lf\n",runtime);
}
   double avg = time/i;
   printf("Avg readtime = %lf\n",avg);
   fclose(fp);
   
   return(0);
}