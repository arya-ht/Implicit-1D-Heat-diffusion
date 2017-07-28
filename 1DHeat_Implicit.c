#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#pragma warning(disable:4996)
// created by Arya HajiTaheri
#define tMax 1
//M_PI works fine in codeblocks
#define dx (M_PI/10)
#define xf M_PI
#define xi 0
#define D 1.0
#define F 0.0
#define n 1.0

#define lambda .4
#define dt (lambda*dx*dx/D)
//init
double u0(double);
double exact(double,double);
void heat_implicit();
int main(void){
    heat_implicit();
    return 0;
}
void heat_implicit(){
    //M and N are limits for dt and dx
    int N = (int)((xf-xi)/dx), M = (int)(tMax/dt);
    double u[M+1][N+1];
    double a = 1 + 2*lambda, b = -lambda, c = -lambda;
    double alpha[N],g[N];

    int i, j, t;
    double k, l;

    for(i=0; i<=N; i++)
        u[0][i] = u0((double)(i*dx));

    for(i=0; i<=M; i++){
        u[i][0] = 0;
        u[i][N] = 0;
    }
    alpha[1] = a;
    for(t=1; t<=M; t++)  {
        g[1] = u[t-1][1];
        for(i=2; i<N; i++) {
            alpha[i] = a - c*b/alpha[i-1];
            g[i] = u[t-1][i] - b/alpha[i-1]*g[i-1];
        }

        u[t][N-1] = g[N-1]/alpha[N-1];

        for(i=N-2; i>0; i--)
            u[t][i] = (g[i] - c*u[t][i+1])/alpha[i];
    }
    for(i=0; i<=M; i++) {
        for(j=0; j<=N; j++)  {
            printf("Estimate: %.8lf\n",u[i][j]);
            printf("Exact: %.8lf\n\n",exact(j*dx,i*dt));
        }
        printf("\n");
    }
    printf("\n\n");
//fileIO
    FILE *ft = fopen("t.csv","w"), *fx = fopen("x.csv","w"),*fu = fopen("u.csv","w"),*fuexact = fopen("u_exact.csv","w"),*f_error = fopen("error.csv","w");
    for (k=0; k<=N; k++) {
        fprintf(ft,"%lf,",k);
    }
    for (k=0; k<=M; k++) {
        fprintf(fx,"%lf\n",k);
    }
    double error = 0.0;
    for(i=0; i<=M; i++) {
        for(j=0; j<=N; j++) {
            fprintf(fu,"%lf,",u[i][j]);
            fprintf(fuexact,"%.8lf,",exact(j*dx,i*dt));
            if(i == 0 || i == M || j == 0 || j == N ) {
                error = 0;
            }
            else {
                error = fabs((u[i][j]-(exact(j*dx,i*dt)))/u[i][j])*100;
            }
            fprintf(f_error,"%.8lf,",error);
        }
        fprintf(fu, "\n");
        fprintf(fuexact, "\n");
        fprintf(f_error, "\n");
    }

    fclose(fx);
    fclose(ft);
    fclose(f_error);
    fclose(fu);
    fclose(fuexact);
}
double u0(double x){  return sin(n*x);}

double exact(double x,double t){  return exp(-n*n*t)*sin(n*x);}

