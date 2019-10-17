#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXORD 8
#define MAXPAR 16

#define YES 1
#define NO 0

enum alogithms {EULER, RK2, RK4};

char *prog;

int
main(int argc, char *argv[]){
  double r0, f0, ro, p[MAXPAR], tmax, t, dt, r[MAXORD];
  void euler(int, double*, double, double, double*);
  void rk2(int, double*, double, double, double*);
  void rk4(int, double*, double, double, double*);
  int n, nsteps;
  int al;

  prog = argv[0];



  printf("Welcome to the rabbit/fox population simulator.\n\n");

  printf("\nRatio of initial rabbit population density to critical population density:\n");
  scanf("%le", r0);

  printf("\nRatio of initial fox population density to critical population density:\n");
  scanf("%le", f0);

  printf("\nRatio of fox starvation rate to rabbit growth rate:\n");
  scanf("%le", p);

  printf("\nAlgorithm to be used (euler, rk2, rk4):\n");
  scanf("%s", al);

  if(al == "rk2"){
    al = RK2;
  }
  else if(al == "rk4"){
    al = RK4;
  }
  else if(al == "euler"){
    al = EULER;
  }


  printf("\nMax time:\n")
  scanf("%le", tmax);

  r[0] = r0; // initial rabbit density to critical density ratio
  r[1] = f0; // initial fox density to critical density ratio
  p[0] = ro;  // (growth rabbits)/(growth foxes)

  dt = tmax/nsteps;
  t = 0.0;

  for (n = 0; n < nsteps; n++) {  // N.B.: start counting at n=0 so that ...

    t = n*dt; // ... "t" is the time at the BEGINNING of the interval!

    /*
    * advance solution from t --> t+dt:
    */
    switch (al) {

      case EULER:
        euler(2, r, t, dt, p);
      break;

      case RK2:
        rk2(2, r, t, dt, p);
      break;

      case RK4:
        rk4(2, r, t, dt, p);
      break;

      default: // this shouldn't happen ...
        fprintf(stderr, "%s: invalid algorithm?\n", prog);
        return 5;
      break;

    }

    t += dt;


  }

}


void
rk4(int ord, double *r, double t, double dt, double *p){

  double drdt[MAXORD], rTemp[MAXORD], k1[MAXORD], k2[MAXORD], k3[MAXORD], k4[MAXORD];

  void model(int, double*, double*, double, double*);

  model(ord, r, k1, t, p);

  for(int i = 0; i < ord; i++){
    rTemp[i] = r[i] + dt*k1[i]/2;
  }

  model(ord, rTemp, k2, t+dt/2, p);

  for(int i = 0; i < ord; i++){
    rTemp[i] = r[i] + dt*k2[i]/2;
  }

  model(ord, rTemp, k3, t+dt/2, p);

  for(int i = 0; i < ord; i++){
    rTemp[i] = r[i] + dt*k3[i];
  }

  model(ord, rTemp, k4, t+dt, p);

  for(int i = 0; i < ord; i++){
    r[i] += dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
  }

  return;
}

void
rk2(int ord, double *r, double t, double dt, double *p){

  double drdt[MAXORD], r12[MAXORD];

  void model(int, double*, double*, double, double*);
  model(ord, r, drdt, t, p);

  for(int i = 0; i < ord; i++){
    r12[i] = r[i] + dt*drdt[i]/2;
  }

  model(ord, r12, drdt, t+dt/2, p);

  for(int i = 0; i < ord; i++){
    r[i] += dt*drdt[i];
  }

  return;
}

void
euler(int ord, double *r, double t, double dt, double *p){

  double drdt[MAXORD];

  void model(int, double*, double*, double, double*);
  model(ord, r, drdt, t, p);

  for(int i = 0; i < ord; i++){
    r[i] += dt*drdt[i];
  }

  return;
}

void
model(int ord, double *r, double *drdt, double t, double *p){

  int i;

  drdt[0] = r[0]*(1-r[1]); // dn_r/dtau = n_R*(1 - n_F)
  drdt[1] = p[0]*r[1]*(r[0]-1); // dn_F/dtau = ro*n_F*(n_R - 1)

  return;
}
