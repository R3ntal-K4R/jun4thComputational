#include "generator.c"
#include <stdio.h>
void main(){
    printf("running main");
    const char inputname[] = "Test.txt";
    DISTRIBUTION *d = distFromFile(inputname);
    const char filename[] = "oupt.txt";
    double distance = 1;
    double sx = 1;
    double sy = 1;
    double sz = 1;


    generateParticle(d, filename, distance, sx, sy, sz);}