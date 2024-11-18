#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>       
#include <omp.h>
#include <unistd.h>     
#include "evil.h"
#include <string.h>
#define n 100
#define epsilon 1e-4
#define delta 1e-4
#include <stdint.h>
// temp vars;
double t[n][n][2] = {0};
double t1[n][n][2] = {0};
double t2[n][n][2] = {0};
double qi[n][2];
double qj[n][2]; 
double tt[2];
double eig[n][2];
double null[2] = {0,0};
double b[2]={0}, c[2]={0};
double roots[2][2] = {0};
double currpd[n] = {0};
double currapd[n-1] = {0};
double currbpd[n-1] = {0};
double prevpd[n] = {0};
double prevapd[n-1] = {0};
double prevbpd[n-1] = {0};
double diffpd[n] = {0};
double diffapd[n-1] = {0};
double diffbpd[n-1] = {0};

//    USE ulimit -s unlimited in UNIX BASED SYSTEMS FOR n > 800
//    NOT A RECOMMENDED SOLUTION BUT SOMETHING BETTER THAN NOTHING

// funcs
void nuller(int N,double a[N]){
  for(int i=0;i<N;i++){
    a[i] = 0;
  }
}
void matmult(double a[n][n][2],double b[n][n][2], double p[n][n][2]){
        for(int i=0;i<n;i++){
          for(int j=0;j<n;j++){
            nuller(2,t[i][j]);
          }
        }
        int i,j,k;double temp[2];
        //#pragma omp parallel for private(i,j,k) shared(a,b,t)
	  for(k=0;k<n;k++){
		  for(i=0;i<n;i++){
			  for(j=0;j<n;j++){
			    
				  evilequal(t[i][j],eviladd(t[i][j],evilmult(a[i][k],(b[k][j]))));
				  //evilequal(t[i][j],temp);
			  }
		  }
	  }
	for(int i=0;i<n;i++){
	  for(int j=0;j<n;j++){
	    p[i][j][0] = t[i][j][0];
	    p[i][j][1] = t[i][j][1];
	  }
	}
}

void matscale(double a[n][n][2], double k){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			evilequal(a[i][j],evilscale(a[i][j],k));
		}
	}
}
void transmat(double a[n][n][2], double at[n][n][2]){
        //double temp[n][n][2];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			t[i][j][0] = a[j][i][0];
			t[i][j][1] = a[j][i][1];
		}
	}
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			at[i][j][0] = t[i][j][0];
			at[i][j][1] = t[i][j][1];
		}
	}
}
void matadd(double a[n][n][2], double b[n][n][2], double c[n][n][2]){
        //double temp[n][n][2];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			t1[i][j][0] = a[i][j][0]+ b[i][j][0];
			t1[i][j][1] = a[i][j][1]+ b[i][j][1];
                        c[i][j][0] = t1[i][j][0];
			c[i][j][1] = t1[i][j][1];
		}
	}/*
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			c[i][j][0] = t[i][j][0];
			c[i][j][1] = t[i][j][1];
		}
	}*/
}
void matsub(double a[n][n][2], double b[n][n][2], double c[n][n][2]){
        //double temp[n][n][2];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			t2[i][j][0] = a[i][j][0]- b[i][j][0];
			t2[i][j][1] = a[i][j][1]- b[i][j][1];
			c[i][j][0] = t2[i][j][0];
			c[i][j][1] = t2[i][j][1];
		}
	}/*
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			c[i][j][0] = t2[i][j][0];
			c[i][j][1] = t2[i][j][1];
		}
	}*/
}
void eye(double a[n][n][2], double k[2]){
	for(int i=0; i<n;i++){
		for(int j=0;j<n;j++){
			a[i][j][0] = (i==j)?k[0]:0;
			a[i][j][1] = (i==j)?k[1]:0;
		}
	}
}
void ithcoltoarr(int i, double qi[n][2], double q[n][n][2]){
	for(int j=0; j<n;j++){
		q[j][i][0] = qi[j][0];
		q[j][i][1] = qi[j][1];
		//evilequal(q[j][i],qi[j])
        }
}

void arrtoithcol(int i, double qi[n][2], double q[n][n][2]){
	for(int j=0;j<n;j++){
		qi[j][1] = q[j][i][1];
		qi[j][0] = q[j][i][0];
	}
}
double *inner(double a[n][2], double b[n][2]){
	//double *sum= (double *)malloc(2*sizeof(double));
	tt[0] = 0; tt[1] = 0;
	for(int i=0;i<n;i++){
		evilequal(tt,eviladd(tt,evilmult(evilcon(a[i]),b[i])));
	}
	return tt;
}
double norm(double a[n][2]){
	double sum=0;
	for(int i=0;i<n;i++){
		sum += (a[i][1]*a[i][1]) + a[i][0]*a[i][0];
	}
	return sqrt(sum);
}
void qr(double a[n][n][2], double q[n][n][2], double r[n][n][2]){
	for(int i=0; i<n; i++){
		arrtoithcol(i, qi, a);
		for(int j=0; j<i; j++){
			arrtoithcol(j, qj, q);
			evilequal(r[j][i],evilscale(inner(qj, qi),1));
			for (int k = 0; k < n; k++) {
                            evilequal(tt, evilsub(qi[k], evilmult(r[j][i], qj[k])));
                            evilequal(qi[k], tt);
                        }
		}
		r[i][i][0] = norm(qi);
		if(r[i][i][0]==0){
		  for(int j=0;j<n;j++){
		    evilequal(qi[j],evilscale(qi[j],0));
		  }
		}
		else{
		  for(int k=0; k<n; k++){
			  evilequal(qi[k],evilscale(qi[k], 1/r[i][i][0]));
		  }
		}
		ithcoltoarr(i, qi, q);
	}
}
void pdsub(void){
  for(int i=0;i<n;i++){
    if(i<n-1){
      diffbpd[i] = currbpd[i] - prevbpd[i];
      diffapd[i] = currapd[i] - prevapd[i];
      diffpd[i] = currpd[i] - prevpd[i];
    }
    else diffpd[i] = currpd[i] - prevpd[i];
  }
}

void extractpd(double A[n][n][2]){
  for(int i=0;i<n;i++){
    if(i<n-1){
      currbpd[i] = sqrt(evilnormsq(A[i+1][i]));
      currapd[i] = sqrt(evilnormsq(A[i][i+1]));
      currpd[i] = sqrt(evilnormsq(A[i][i+1]));
    } 
    else currpd[i] = sqrt(evilnormsq(A[i][i+1]));
  }
}

void updatepd(void){
  for(int i=0;i<n;i++){
    if(i<n-1){
      prevbpd[i] = currbpd[i];
      prevapd[i] = currapd[i];
      prevpd[i] = currpd[i];
    } 
    else prevpd[i] = currpd[i];
  }
}
void printmat(double A[n][n][2]){
  for(int i=0;i<n;i++){
    for(int j=0; j<n;j++){
      printf("(%lf , %lf) ",A[i][j][0],A[i][j][1]);
    }
    printf("||\n");
    printf("\n");
  }
}
void printvec(double A[n][2]){
    for(int j=0; j<n;j++){
      printf("(%e , %e)\n",A[j][0],A[j][1]);
    }
}
int isHessenberg(double A[n][n][2]) {
    for (int i = 2; i < n; i++) {
        for (int j = 0; j < i - 1; j++) {
            if (A[i][j][0] != 0 || A[i][j][1] != 0) {
                return 0;
            }
        }
    }
    return 1;
}
int isConverged(double A[n][n][2]){
  extractpd(A);
  pdsub();
  updatepd();
  for(int i=0;i<n;i++){
    if(i<n-1){
      if(fabs(diffpd[i]) > epsilon){
      //|| fabs(diffapd[i]) > epsilon || fabs(diffbpd[i]) > epsilon){
          //printf("0");
          return 0;
      }
    }
  }
  return 1;
}

void deflate(double A[n][n][2]){
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(fabs(sqrt(evilnormsq(A[i][j]))) < delta) evilequal(A[i][j], null);
      else continue;
    }
  }		
}
int schur(double A[n][n][2]){
  double Q[n][n][2]={0}, R[n][n][2]={0}; int i= 0;
  while(i<10000000){
    qr(A,Q,R);
    matmult(R,Q,A);
    deflate(A);
    if(isConverged(A)){
      return i+1;
    }
    else i++;
  }
  return i;
}
int schurshift(double A[n][n][2]){
  double Q[n][n][2]={0}, R[n][n][2]={0}, U[n][n][2]; int i=0;
  while(i<10000){
    //transmat(A,A);
    eye(U,A[n-1][n-1]);
    matsub(A,U,A);
    qr(A,Q,R);
    matmult(R,Q,A);
    matadd(A,U,A);
    //deflate(A);
    if(isConverged(A) || isHessenberg(A)){
      return i+1;
    }
    else i++;
  }
  return i;
}

void root(void){
  //printf("%lf\n",dis[0]);
  evilequal(roots[0], eviladd(b,evilsqrt(evilsub(evilmult(b,b),evilscale(c,4)))));
  evilequal(roots[1], evilsub(b,evilsqrt(evilsub(evilmult(b,b),evilscale(c,4)))));
  evilequal(roots[0],evilscale(roots[0], (double)1/2));
  evilequal(roots[1],evilscale(roots[1], (double)1/2));
}

void eigen(double a[n][n][2]){
  int i=0;
  while(i<n){
    if(isevilnull(a[i+1][i]) || isevilnull(a[i][i+1])){
      evilequal(eig[i], a[i][i]);i++;
      continue;
    }
    else{
      evilequal(tt, evilmult(a[i][i], a[i+1][i+1]));
      evilequal(c, evilsub(tt,evilmult(a[i+1][i], a[i][i+1]))); tt[0] = 0; tt[1] = 0;
      //printf("%lf\n",c[0]);
      evilequal(b, eviladd(a[i][i], a[i+1][i+1]));
      root();
      evilequal(eig[i], roots[0]);evilequal(eig[i+1], roots[1]);
      i+=2;
    }
  }
}
int main(void){
  //double A[n][n][2] = {{{4,0},{1,0},{1,0},{7,4}},{{8,1},{3,3},{4,2},{3,4}},{{7,2},{5,3},{6,1},{2,9}},{{1,3},{4,7},{8,0},{2,3}}};
  //double A[n][n][2] = {{{1,0},{1,0},{0,0}},{{1,0},{0,0},{1,0}},{{0,0},{1,0},{1,0}}};
  //double A[2][2][2] = {{{1,0},{1,0}},{{-4,0},{-2,0}}};
  
  FILE *out = fopen("out.dat", "a");
  double A[n][n][2];
  srand(time(NULL));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j][0] = (double)rand()/10;
            A[i][j][1] = (double)rand()/10;
        }
    }
  //printmat(A);
  int k=0;
  double start_time = omp_get_wtime();
  if(n>2){
    k = schurshift(A);
  }
  printmat(A);
  eigen(A);
  double end_time = omp_get_wtime(); 
  double time_taken = end_time - start_time;
  printvec(eig);
  fprintf(out, "%d %lf\n", n, time_taken);
  fclose(out);
  printf("schur() took %f seconds and %d iterations to execute \n", time_taken,k);
}
