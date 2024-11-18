double temp1[2];
double temp2[2];
double stemp[2];
int isevilequal(double z1[2], double z2[2]){
  if(z1[0] == z2[0] && z1[1] == z2[1]) return 1;
  else return 0;
}
int isevilnull(double z1[2]){
  if(fabs(z1[0]) < pow(10,-20) && fabs(z1[1]) < pow(10,-20)){
    return 1;
  }
  else return 0;
}
void evilequal(double z1[2], double z2[2]){
    z1[0] = z2[0];
    z1[1] = z2[1];
}

double evilnormsq(double z1[2]){
    return z1[0]*z1[0] +(z1[1]*z1[1]);
}
double *evilmult(double z1[2],double z2[2]){
    //double *temp = (double *)malloc(2*sizeof(double));
    stemp[0] = z1[0]*z2[0] - (z1[1]*z2[1]);
    stemp[1] = z1[0]*z2[1] + (z1[1]*z2[0]);
    return stemp;
    //free(temp);
}
double *eviladd(double z1[2],double z2[2]){
    //double *temp = (double *)malloc(2*sizeof(double));
    temp1[0] = z1[0]+z2[0];
    temp1[1] = z1[1]+z2[1];
    return temp1;
    //free(temp);
}
double *evilsub(double z1[2],double z2[2]){
    //double *temp = (double *)malloc(2*sizeof(double));
    temp1[0] = z1[0]-z2[0];
    temp1[1] = z1[1]-z2[1];
    return temp1;
    //free(temp);
}
double *evildivi(double z1[2], double z2[2]){
    //double *temp = (double *)malloc(2*sizeof(double));
    if(z1[0]==1 && z2[0]==0){
        temp1[0] = z2[0];
        temp1[1] = z2[1];
        return temp1;
    }
    double z2n = evilnormsq(z2);
    temp1[0] = (z1[0]*z2[0] + (z1[1]*z2[1]))/z2n;
    temp1[1] = (z1[1]*z2[0] - (z1[0]*z2[1]))/z2n;
    return temp1;
    //free(temp);
}
double *evilcon(double z1[2]){
    //double *temp = (double *)malloc(2*sizeof(double));
    temp2[0] = z1[0];
    temp2[1] = -z1[1]; 
    return temp2;
    //free(temp);
}
double *evilscale(double a[2],double k){
  //double *temp = (double *)malloc(2*sizeof(double));
  temp1[0] = a[0]*k;
  temp1[1] = a[1]*k;
  return temp1;
  //free(temp);
} 
double *evilsqrt(double z1[2]) {
    double magnitude = sqrt(sqrt(evilnormsq(z1)));
    double angle = atan2(z1[1], z1[0]) / 2; 
    stemp[0] = magnitude * cos(angle);
    stemp[1] = magnitude * sin(angle);
    return stemp;
}
