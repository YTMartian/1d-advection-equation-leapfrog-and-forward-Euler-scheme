/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"


void shockwave(const char *q1_file) {
    printf("shockwave() - IMPLEMENT ME!\n");
}

void linsolve(const char *q3_file) {
    printf("linsolve() - IMPLEMENT ME!\n");
}


#define TYPE double

void interp(const char *q4_file, double xo) {
    //read csv.
    FILE *stream = fopen(q4_file, "r");
    char line[1024];
    char tmp[1024];
    int first_line = 1;
    int N = 0;
    //first read to calculate N.
    while (fgets(line, 1024, stream)) {
        if(first_line == 1) {
            first_line = 0;
            continue;
        }
        N++;
    }

    TYPE *x = (TYPE *)calloc(N, sizeof(TYPE));
    TYPE *a = (TYPE *)calloc(N, sizeof(TYPE));//a_i = y_i
    TYPE *b = (TYPE *)calloc(N - 1, sizeof(TYPE));
    TYPE *c = (TYPE *)calloc(N, sizeof(TYPE));
    TYPE *d = (TYPE *)calloc(N - 1, sizeof(TYPE));
    TYPE *h = (TYPE *)calloc(N - 1, sizeof(TYPE));//h_i = x_i+1 - x_i.
    int n = 0;
    //second read to get data.
    fclose(stream);
    stream = fopen(q4_file, "r");
    first_line = 1;
    while (fgets(line, 1024, stream)) {
        if(first_line == 1) {
            first_line = 0;
            continue;
        }
        strcpy(tmp, line); //copy string.
        const char *tok = NULL;
        int i = 0;
        for (tok = strtok(tmp, ","); tok && *tok; tok = strtok(NULL, ",\n"), i++) {
            if (i == 0) x[n] = atof(tok);
            else if (i == 1) a[n] = atof(tok);
        }
        n++;
    }
    fclose(stream);
    //printf("%d*\n", N);

    //method from Numerical Alalysis Ninth Edition.
    //calculate h_i.
    for(int i = 0; i < N - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }
    //calculate Equation(55) right half part.
    TYPE *A = (TYPE *)calloc(N - 1, sizeof(TYPE));
    for(int i = 0; i < n - 1; i++) {
        A[i] = 3 * (a[i + 1] - a[i]) / h[i] + 3 * (a[i - 1] - a[i]) / h[i - 1];
    }
    //LU Factorization.
    TYPE *L = (TYPE *)calloc(N, sizeof(TYPE));
    TYPE *U = (TYPE *)calloc(N, sizeof(TYPE));
    TYPE *Z = (TYPE *)calloc(N, sizeof(TYPE));
    L[0] = 1;
    U[0] = Z[0] = 0;
    for(int i = 1; i < n - 1; i++) {
        L[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * U[i - 1];
        U[i] = h[i] / L[i];
        Z[i] = (A[i] - h[i - 1] * Z[i - 1]) / L[i];
    }
    L[n - 1] = 1;
    Z[n - 1] = c[n - 1] = 0;
    for(int i = n - 2; i >= 0; i--) {
        c[i] = Z[i] - U[i] * c[i + 1];
        //Equation(54) calculate b_i.
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
        //Equation(52) calculate d_i.
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    /*
    for (int i = 0; i < n - 1; ++i) {
        double fxo = a[i] + b[i] * (xo - x[i]) + c[i] * (xo - x[i]) * (xo - x[i]) + d[i] * (xo - x[i]) * (xo - x[i]) * (xo - x[i]);
        printf("%2d %18.6lf %18.6lf %18.6lf %18.6lf %18.6lf\n", i, a[i], b[i], c[i], d[i], fxo);
    }
    */
    /*
    //out detailed data.
    FILE *fp2 = fopen("out_interp2.csv","w");
    for(int i = 0; i < N - 1; i++) {
        TYPE t = x[i + 1] - x[i];
        TYPE xx = x[i];
        TYPE t2 = t / 10.0;
        for(int j = 0; j < 10; j++) {
            TYPE yy = a[i] + b[i] * (xx - x[i]) + c[i] * (xx - x[i]) * (xx - x[i]) + d[i] * (xx - x[i]) * (xx - x[i]) * (xx - x[i]);
            fprintf(fp2, "%.6lf,%.6lf\n", xx, yy);
            xx += t2;
        }
    }
    fclose(fp2);
    */

    FILE *fp = fopen("out_interp.csv", "w");
    fprintf(fp, "xo,f(xo)\n");
    //to determine the interval xo in.
    for(int i = 0; i < N - 1; i++) {
        if((xo >= x[i] && xo <= x[i + 1]) || (xo <= x[i] && xo >= x[i + 1])) {
            TYPE fxo = a[i] + b[i] * (xo - x[i]) + c[i] * (xo - x[i]) * (xo - x[i]) + d[i] * (xo - x[i]) * (xo - x[i]) * (xo - x[i]);
            fprintf(fp, "%.6lf,%.6lf\n", xo, fxo);
        }
    }
    fclose(fp);
    free(x);
    free(a);
    free(b);
    free(c);
    free(d);
    free(h);
    free(A);
    free(L);
    free(U);
    free(Z);
}


void advection(const char *q5_file) {
    const double PI = 3.14159265358979;
    TYPE u, CFL, t_final;
    int Nx;
    //read csv file.
    FILE *stream = fopen(q5_file, "r");
    char line[1024];
    char tmp[1024];
    int first_line = 1;
    while (fgets(line, 1024, stream)) {//fgets:get one line.
        if(first_line == 1) {
            first_line = 0;
            continue;
        }
        strcpy(tmp, line);//copy string.
        char *tok = NULL;
        int i = 0;
        for (tok = strtok(tmp, ","); tok && *tok; tok = strtok(NULL, ","), i++) {
            if (i == 0) u = atof(tok);
            else if (i == 1) Nx = atoi(tok);
            else if(i == 2)CFL = atof(tok);
            else if(i == 3)t_final = atof(tok);
        }
        //free(tmp);
    }
    //printf("%lf %d %lf %lf\n", u, Nx, CFL, t_final);

    TYPE delta_x = 1.0 / (double)Nx;
    TYPE delta_t = (delta_x * CFL) / u;
    int Nt = (int)(round(t_final / delta_t));

    TYPE *array_x = (TYPE *)calloc(Nx + 1, sizeof(TYPE));

    //printf("delta_x=%lf\n", delta_x);
    for(int i = 0; i < Nx + 1; i++) {
        array_x[i] = i * delta_x;
    }

    TYPE **phi = (TYPE **)calloc(Nt + 1, sizeof(TYPE *));
    for(int i = 0; i < Nt + 1; i++) {
        phi[i] = (TYPE *)calloc(Nx + 1, sizeof(TYPE));
    }
    //initialize phi
    for(int i = 0; i < Nx + 1; i++) {
        if(0 <= array_x[i] && array_x[i] < 0.125) phi[0][i] = 0;
        else if(0.125 <= array_x[i] && array_x[i] <= 0.375) phi[0][i] = 0.5 * (1 - cos(8.0 * PI * (array_x[i] - 0.125)));
        else if(0.375 < array_x[i] && array_x[i] <= 1.0) phi[0][i] = 0;
    }

    //first step: forward Euler scheme
    TYPE half_CFL = CFL / 2.0;
    phi[1][0] = phi[0][0] - half_CFL * (phi[0][1] - phi[0][Nx - 1]);//at the x_0.
    for(int i = 1; i < Nx; i++) {
        phi[1][i] = phi[0][i] - half_CFL * (phi[0][i + 1] - phi[0][i - 1]);
    }
    phi[1][Nx] = phi[0][Nx] - half_CFL * (phi[0][1] - phi[0][Nx - 1]);//at the x_Nx.

    //printf("%d***\n",Nt);
    //subsequent steps: leapfrog scheme
    for(int n = 1; n < Nt; n++) {
        phi[n + 1][0] = phi[n - 1][0] - CFL * (phi[n][1] - phi[n][Nx - 1]);//at the x_0.
        for(int i = 1; i < Nx; i++) {
            phi[n + 1][i] = phi[n - 1][i] - CFL * (phi[n][i + 1] - phi[n][i - 1]);
        }
        phi[n + 1][Nx] = phi[n - 1][Nx] - CFL * (phi[n][1] - phi[n][Nx - 1]);//at the x_Nx.
    }

    //save to csv.
    FILE *fp = fopen("out_advection.csv", "w");
    fprintf(fp, "x,phi\n");
    for(int i = 0; i <= Nx; i++) {
        fprintf(fp, "%.6lf,%.6lf\n", delta_x * i, phi[Nt][i]);
    }
    fclose(fp);

    free(array_x);
    //free(array_t);
    for(int i = 0; i < Nt; i++) {
        free(phi[i]);
    }
    free(phi);
}






