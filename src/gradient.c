#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

// calculate the gradient of a scaller over element in tri3 mesh
int gradient_ele_tri3(int nelem, int *elems, double *ptxyz, int *esure, double *normedge,double *u,double **grad2){
    double *grad = (double *)malloc((size_t)nelem*3*sizeof(double));
    if (grad==NULL){
        fprintf(stderr,"! there is error in memory allocation for grad pointer.\n");
        return -1;
    }
    for (int ele = 0;ele<nelem;ele++){
        // Coordinates of points in a triangular element
		int p1 = elems [3*ele] -1;
		int p2 = elems [3*ele+1] -1;
		int p3 = elems [3*ele+2] -1;
		double cp1[3] = {ptxyz[3*p1], ptxyz[3*p1+1], ptxyz[3*p1+2]};
		double cp2[3] = {ptxyz[3*p2], ptxyz[3*p2+1], ptxyz[3*p2+2]};
		double cp3[3] = {ptxyz[3*p3], ptxyz[3*p3+1], ptxyz[3*p3+2]};
        double a1 = SQUARE (cp1[0]-cp2[0]);
        a1 += SQUARE (cp1[1]-cp2[1]);
        a1 += SQUARE (cp1[2]-cp2[2]);
        a1 = sqrt(a1);
        double a2 = SQUARE (cp2[0]-cp3[0]);
        a2 += SQUARE (cp2[1]-cp3[1]);
        a2 += SQUARE (cp2[2]-cp3[2]);
        a2 = sqrt(a2);
        double a3 = SQUARE (cp3[0]-cp1[0]);
        a3 += SQUARE (cp3[1]-cp1[1]);
        a3 += SQUARE (cp3[2]-cp1[2]);
        a3 = sqrt(a3);
        double u1 = (u[esure[3*ele]] + u[ele])/2;
        double u2 = (u[esure[3*ele+1]] + u[ele])/2;
        double u3 = (u[esure[3*ele+2]] + u[ele])/2;
        double flux_x = u1 * normedge [9*ele] * a1;
        flux_x += u2 * normedge [9*ele+3] * a2;
        flux_x += u3 * normedge [9*ele+6] * a3;
        double flux_y = u1 * normedge [9*ele+1] * a1;
        flux_y += u2 * normedge [9*ele+3+1] * a2;
        flux_y += u3 * normedge [9*ele+6+1] * a3;
        double flux_z = u1 * normedge [9*ele+2] * a1;
        flux_z += u2 * normedge [9*ele+3+2] * a2;
        flux_z += u3 * normedge [9*ele+6+2] * a3;
        grad [3*ele] = flux_x;
        grad [3*ele+1] = flux_y;
        grad [3*ele+2] = flux_z;
    }
    *grad2 = grad;
    return 0; // success signal
}