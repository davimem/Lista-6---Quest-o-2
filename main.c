#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x, double y, double yl)
{
  return y*y*y - y*yl;
}

double fy(double x, double y, double yl)
{
  return 3*y*y - yl;
}

double fyl(double x, double y, double yl)
{
  return -y;
}


double *NewtonNL(double a, double b, int N, double alfa, double beta)
{
  int i;
  double t, yl, norma;
  double h=(b-a)/(N+1);
  double h2=h*h;
  
  double *D=malloc(sizeof(double)*(N+1));
  double *U=malloc(sizeof(double)*(N+1));
  double *L=malloc(sizeof(double)*(N+1));
  double *R=malloc(sizeof(double)*(N+1));
  double *l=malloc(sizeof(double)*(N+1));
  double *u=malloc(sizeof(double)*(N+1));
  double *z=malloc(sizeof(double)*(N+1));
  double *y=malloc(sizeof(double)*(N+2));
  double *v=malloc(sizeof(double)*(N+2));

  for(i=0; i<N+2; i++)y[i]=alfa + (beta-alfa)/(b-a) * (i*h);
  
  do{
    
    t = a + h;

    double yl = (y[2] - y[0]) / (2 * h);
    D[1]      = 2 + h2 * fy(t, y[1], yl);
    U[1]      = -1 + h / 2 * fyl(t, y[1], yl);
    R[1]      = -(2 * y[1] - y[2] + h2 * f(t, y[1], yl) - alfa);

    for (int i = 2; i < N; i++) {
      t += h;
      yl   = (y[i + 1] - y[i - 1]) / (2 * h);
      D[i] = 2 + h2 * fy(t, y[i], yl);
      U[i] = -1 + h / 2 * fyl(t, y[i], yl);
      L[i] = -1 - (h / 2) * fyl(t, y[i], yl);
      R[i] = -(2 * y[i] - y[i + 1] - y[i - 1] + h2 * f(t, y[i], yl));
    }

    t += h;
    yl   = (y[N] - y[N - 1]) / (2 * h);
    D[N] = 2 + h2 * fy(t, y[N], yl);
    L[N] = -1 - (h / 2) * fyl(t, y[N], yl);
    R[N] = -(2 * y[N] - y[N - 1] - beta + h2 * f(t, y[N], yl));

    memcpy(v, y, (N + 2) * sizeof(double));
  
    l[1] = D[1];
    u[1] = U[1]/D[1];
    z[1] = R[1]/l[1];
    
    for(i=2; i<N; i++) 
    {
      l[i] = D[i] - L[i]*u[i-1];
      u[i] = U[i]/l[i];
      z[i] = (R[i] - L[i]*z[i-1])/l[i];
    }
    
    l[N] = D[N] - L[N]*u[N-1];
    z[N] = (R[N] - L[N]*z[N-1])/l[N];
    
    v[N] = z[N];  
    
    for (int i = N - 1; i > 0; i--)
    {
      v[i] = z[i] - u[i] * v[i + 1];
    }
    
    norma=0;
    for (i=N-1; i>0; i--) 
    {
      if(fabs(v[i])>norma) norma=fabs(v[i]);
      y[i]+=v[i];
    }

  }while((norma>1e-5));
  return y;
}

  
int main(int argc, char **argv)
{
  double a, b, alfa, beta, h, *x;
  int N,i;
  
  a=1;
  b=2;
  alfa=1/2.;
  beta=1/3.;
  N=9;
  h=(b-a)/(N+1);
  
  
  x=calloc(N+2, sizeof(double));

  x=NewtonNL(a, b, N, alfa, beta);
  
  for (i=0; i<N+2; i++) printf("%lf %lf\n", a+i*h, fabs(x[i]-(1/(a+i*h+1))));
  return 0;
}