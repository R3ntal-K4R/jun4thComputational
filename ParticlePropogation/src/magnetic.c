#include <auxiliary.h>
#include <magnetic.h>

/*computes the magnetic field the particle interacts with at a given point
  inputs: particle p, spacecraft s
  outputs: vector of x,y,z, components of magnetic field felt by particle
*/
vector MagField(particle p, spacecraft s) {
    // checkbounds(p, s);
    //find index position
    // printf("%f %f %f\n",p.x,p.y,p.z);
    double i2 = floor(p.x/s.deltaB);
    double j2 = floor(p.y/s.deltaB);
    double k2 = floor(p.z/s.deltaB);
    // printf("%f %f %f\n",i2,j2,k2);
    //difference from index position
    double deltax = p.x - i2 * s.deltaB;
    double deltay = p.y - j2 * s.deltaB;
    double deltaz = p.z - k2 * s.deltaB;
    int i = i2;
    int j = j2;
    int k = k2;
    //interpolate magnetic field for particle position
    // Keegan thinks we shouldn't have /deltaB in these terms.
    // printf("bx\n");
    // printf("%d %d %d\n",i,j,k);
    double Bx = s.BFx[i][j][k] + (deltax/s.deltaB)*(s.BFx[i+1][j][k]-s.BFx[i][j][k]) + (deltay/s.deltaB)*(s.BFx[i][j+1][k]-s.BFx[i][j][k]) + (deltaz/s.deltaB)*(s.BFx[i][j][k+1]-s.BFx[i][j][k]);
    // printf("end bx\n");
    double By = s.BFy[i][j][k] + (deltax/s.deltaB)*(s.BFy[i+1][j][k]-s.BFy[i][j][k]) + (deltay/s.deltaB)*(s.BFy[i][j+1][k]-s.BFy[i][j][k]) + (deltaz/s.deltaB)*(s.BFy[i][j][k+1]-s.BFy[i][j][k]);
    double Bz = s.BFz[i][j][k] + (deltax/s.deltaB)*(s.BFz[i+1][j][k]-s.BFz[i][j][k]) + (deltay/s.deltaB)*(s.BFz[i][j+1][k]-s.BFz[i][j][k]) + (deltaz/s.deltaB)*(s.BFz[i][j][k+1]-s.BFz[i][j][k]);
    vector returnVec;
    //multiply by scale factor
    returnVec.x = Bx*s.B_multiplier;
    returnVec.y = By*s.B_multiplier;
    returnVec.z = Bz*s.B_multiplier;
    return returnVec;
}
