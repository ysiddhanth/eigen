#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>       
#include <unistd.h>     
#include "evil.h"
#define n 800
#define epsilon pow(10,-4)
#define delta pow(10,-15)

// temp vars;
  double ***t;
  double ***t1;
  double ***t2;
  double **qi;
  double **qj;
  double **eig;
  double tt[2];
  double null[2] = {0,0};
  double b[2]={0}, c[2]={0};
  double roots[2][2] = {0};
  double *currpd;
  double *prevpd;
  double *diffpd;
  double *currbpd;
  double *prevbpd;
  double *diffbpd;
  double *currapd;
  double *prevapd;
  double *diffapd;
//func def
double ***createMat(void);
double **createArr(int N);
double *createRArr(int N);
void nuller(int N,double *a);
void matmult(double ***a,double ***b, double ***p);
void matscale(double ***a, double k);
void transmat(double ***a, double ***at);
void matadd(double ***a, double ***b, double ***c);
void matsub(double ***a, double ***b, double ***c);
void eye(double ***a, double *k);
void ithcoltoarr(int i, double **qi, double ***q);
void arrtoithcol(int i, double **qi, double ***q);
double *inner(double **a, double **b);
double norm(double **a);
void qr(double ***a, double ***q, double ***r);
void pdsub(void);
void extractpd(double ***A);
void updatepd(void);
void printmat(double ***A);
void printvec(double **A);
int isConverged(double ***A);
void deflate(double ***A);
int schur(double ***A);
void root(void);
void eigen(double ***a);
// main
int main(void){
  // temp vars;
  t = createMat();
  t1 = createMat();
  t2 = createMat();
  qi = createArr(n);
  qj = createArr(n);
  eig = createArr(n);
  currpd = createRArr(n);
  prevpd = createRArr(n);
  diffpd = createRArr(n);
  currbpd = createRArr(n-1);
  prevbpd = createRArr(n-1);
  diffbpd = createRArr(n-1);
  currapd = createRArr(n-1);
  prevapd = createRArr(n-1);
  diffapd = createRArr(n-1);
  clock_t t = clock(); 
  //double A[n][n][2] = {{{4,0},{-2,0},{1,0}},{{2,0},{4,0},{1,0}},{{1,0},{1,0},{3,0}}};
  //double A[2][2][2] = {{{1,0},{1,0}},{{-4,0},{-2,0}}};
  double ***A = createMat();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j][0] = (double)(i + j);
            A[i][j][1] =0;
        }
    }
  int k=0;
  if(n>2){
    k = schur(A);
  }
  //printmat(A);
  eigen(A);
  t = clock() - t; 
  double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
  printvec(eig);
  printf("schur() took %f seconds and %d iterations to execute \n", time_taken,k);
}



//funcs
double ***createMat(void){
  double ***A = (double ***)malloc(n*sizeof(double **));
  for(int i=0;i<n;i++){
    A[i] = (double **)malloc(n*sizeof(double *));
    for(int j=0;j<n;j++){
      A[i][j] = (double *)malloc(2*sizeof(double));
    }
  }
  return A;
}
double **createArr(int N){
  double **A = (double **)malloc(N*sizeof(double *));
  for(int i=0;i<N;i++){
    A[i] = (double *)malloc(2*sizeof(double));
  }
  return A;
}
double *createRArr(int N){
  double *A = (double *)malloc(N*sizeof(double));
  return A;
}
void nuller(int N,double *a){
  for(int i=0;i<N;i++){
    a[i] = 0;
  }
}
void matmult(double ***a,double ***b, double ***p){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			nuller(2,t[i][j]);
			for(int k=0;k<n;k++){
			        double temp[2];
				evilequal(temp,eviladd(t[i][j],evilmult(a[i][k],(b[k][j]))));
				evilequal(t[i][j],temp);
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
void matscale(double ***a, double k){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			evilequal(a[i][j],evilscale(a[i][j],k));
		}
	}
}
void transmat(double ***a, double ***at){
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
void matadd(double ***a, double ***b, double ***c){
        //double temp[n][n][2];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			t1[i][j][0] = a[i][j][0]+ b[i][j][0];
			t1[i][j][1] = a[i][j][1]+ b[i][j][1];
                        c[i][j][0] = t1[i][j][0];
			c[i][j][1] = t1[i][j][1];
		}
	}
}
void matsub(double ***a, double ***b, double ***c){
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
void eye(double ***a, double *k){
	for(int i=0; i<n;i++){
		for(int j=0;j<n;j++){
			a[i][j][0] = (i==j)?k[0]:0;
			a[i][j][1] = (i==j)?k[1]:0;
		}
	}
}
void ithcoltoarr(int i, double **qi, double ***q){
	for(int j=0; j<n;j++){
		q[j][i][0] = qi[j][0];
		q[j][i][1] = qi[j][1];
		//evilequal(q[j][i],qi[j])
        }
}

void arrtoithcol(int i, double **qi, double ***q){
	for(int j=0;j<n;j++){
		qi[j][1] = q[j][i][1];
		qi[j][0] = q[j][i][0];
	}
}
double *inner(double **a, double **b){
	//double *sum= (double *)malloc(2*sizeof(double));
	tt[0] = 0; tt[1] = 0;
	for(int i=0;i<n;i++){
		evilequal(tt,eviladd(tt,evilmult(evilcon(a[i]),b[i])));
	}
	return tt;
}
double norm(double **a){
	double sum=0;
	for(int i=0;i<n;i++){
		sum += (a[i][1]*a[i][1]) + a[i][0]*a[i][0];
	}
	return sqrt(sum);
}
void qr(double ***a, double ***q, double ***r){
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

void extractpd(double ***A){
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
void printmat(double ***A){
  for(int i=0;i<n;i++){
    for(int j=0; j<n;j++){
      printf("(%e , %e) ",A[i][j][0],A[i][j][1]);
    }
    printf("||\n");
    printf("\n");
  }
}
void printvec(double **A){
    for(int j=0; j<n;j++){
      printf("(%e , %e)\n",A[j][0],A[j][1]);
    }
}
int isConverged(double ***A){
  extractpd(A);
  pdsub();
  updatepd();
  for(int i=0;i<n;i++){
    if(i<n-1){
      if(fabs(diffpd[i]) > epsilon || fabs(diffapd[i]) > epsilon || fabs(diffbpd[i]) > epsilon){
          return 0;
      }
    }
  }
  return 1;
}

void deflate(double ***A){
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(fabs(sqrt(evilnormsq(A[i][j]))) < delta) evilequal(A[i][j], null);
      else continue;
    }
  }
}
int schur(double ***A){
  double ***Q = createMat(), ***R= createMat(); int i= 0;
  while(i<1000){
    qr(A,Q,R);
    matmult(R,Q,A);
    deflate(A);
    if(isConverged(A)){
      return i+1;
    }
    else i++;
  }
  return 23;
}
/*
void schurshift(double A[n][n][2]){
  double Q[n][n][2]={0}, R[n][n][2]={0}, U[n][n][2];
  for(int i =0;i<800;i++){
    eye(U,A[n-1][n-1]);
    matsub(A,U,A);
    qr(A,Q,R);
    matmult(R,Q,A);
    matadd(A,U,A);
  }
}*/

void root(void){
  //printf("%lf\n",dis[0]);
  evilequal(roots[0], eviladd(b,evilsqrt(evilsub(evilmult(b,b),evilscale(c,4)))));
  evilequal(roots[1], evilsub(b,evilsqrt(evilsub(evilmult(b,b),evilscale(c,4)))));
  evilequal(roots[0],evilscale(roots[0], (double)1/2));
  evilequal(roots[1],evilscale(roots[1], (double)1/2));
}

void eigen(double ***a){
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
