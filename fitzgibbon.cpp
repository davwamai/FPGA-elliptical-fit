// #include <hls_math.h>
#include "fitzgibbon.h"
#include <cmath> 
#include <iostream>
#define MAX_SIZE 300
#define DIM_1 3
#define DIM_2 3
using namespace std;

void inverse_3x3(double A[3][3], double B[3][3]) {
    double det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) -
                 A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                 A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    double invdet = 1 / det;

    B[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * invdet;
    B[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * invdet;
    B[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * invdet;
    B[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * invdet;
    B[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * invdet;
    B[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * invdet;
    B[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * invdet;
    B[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * invdet;
    B[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * invdet;
}

void solve_cubic(double a, double b, double c, double d, double roots[3]) {
    // convert to depressed cubic t^3 + pt + q = 0 (subst x = t - b/3a)
    double p = c / a - b * b / (3 * a * a);
    double q = 2 * b * b * b / (27 * a * a * a) - b * c / (3 * a * a) + d / a;
    double sub = b / (3 * a);

    double q2 = q / 2;
    double discriminant = q2 * q2 + p * p * p / 27;
    if (discriminant < 0) {
        // three real roots
        double mp3 = -p / 3;
        double mp33 = mp3 * mp3 * mp3;
        double r = sqrt(mp33);
        double t = -q / (2 * r);
        double cosphi = t < -1 ? -1 : t > 1 ? 1 : t;
        double phi = acos(cosphi);
        double crtr = cbrt(r);
        double t1 = 2 * crtr;
        roots[0] = t1 * cos(phi / 3) - sub;
        roots[1] = t1 * cos((phi + M_PI * 2) / 3) - sub;
        roots[2] = t1 * cos((phi + 4 * M_PI) / 3) - sub;
    } else {
        // one real root
        double sd = sqrt(discriminant);
        roots[0] = cbrt(q2 + sd) + cbrt(q2 - sd) - sub;
    }
}

void eigenvalues_3x3(double A[3][3], double eigenvalues[3]) {
    // coefficients of the characteristic equation
    double a = -1;
    double b = A[0][0] + A[1][1] + A[2][2];
    double c = A[0][1]*A[1][0] + A[0][2]*A[2][0] + A[1][2]*A[2][1] - A[0][0]*A[1][1] - A[0][0]*A[2][2] - A[1][1]*A[2][2];
    double d = A[0][0]*A[1][1]*A[2][2] + 2*A[0][1]*A[1][2]*A[2][0] - A[0][0]*A[1][2]*A[1][2] - A[1][1]*A[0][2]*A[0][2] - A[2][2]*A[0][1]*A[0][1];

    // solve the characteristic equation using Cardano's method
    solve_cubic(a, b, c, d, eigenvalues);
}

void eigenvectors_3x3(double A[3][3], double eigenvalues[3], double eigenvectors[3][3]) {
    for (int i = 0; i < 3; i++) {
        double lambda = eigenvalues[i];

        // create the matrix (A - lambda I)
        double B[3][3];
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                B[j][k] = A[j][k];
            }
            B[j][j] -= lambda;
        }

        // solve (A - lambda I)v = 0
        // we can find a solution by setting one of the variables to 1 and solving for the others
        double x = 1;
        double y = -B[0][1] / B[0][0];
        double z = -B[0][2] / B[0][0];

        // normalize the eigenvector
        double norm = sqrt(x*x + y*y + z*z);
        eigenvectors[i][0] = x / norm;
        eigenvectors[i][1] = y / norm;
        eigenvectors[i][2] = z / norm;
    }
}

void transpose_300(double A[300][3], double B[3][300], int rows, int cols) {
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
            B[j][i] = A[i][j];
}

void matmul_300(double A[3][300], double B[300][3], double C[3][3], int rowsA, int colsA, int colsB) {
    for(int i = 0; i < rowsA; i++)
        for(int j = 0; j < colsB; j++)
            for(int k = 0; k < colsA; k++)
                C[i][j] += A[i][k] * B[k][j];
}

void matmul_3x3(double A[3][3], double B[3][3], double C[3][3]) {
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            C[i][j] = 0;
            for(int k = 0; k < 3; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void print_matrix_300x3(double matrix[300][3]) {
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("...\n\n");
}

void print_matrix_3x300(double matrix[3][300]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 10; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("...\n\n");

}

void print_matrix_3x3(double matrix[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("...\n\n");

}

void fitzgibbon(double* x, double* y, int size, double* a) {
    double D1[MAX_SIZE][DIM_1];
    double D2[MAX_SIZE][DIM_2];
    double S1[DIM_1][DIM_1];
    double S2[DIM_1][DIM_2];
    double S3[DIM_2][DIM_2];
    double T[DIM_2][DIM_2];
    double M[DIM_1][DIM_1];
    double evec[DIM_1][DIM_1];
    double eval[DIM_1];
    double D1_t[DIM_1][MAX_SIZE];
    double cond[DIM_1];
    double a1[DIM_1];
    int i;

    for(i = 0; i < size; i++) {
        D1[i][0] = x[i] * x[i];
        D1[i][1] = x[i] * y[i];
        D1[i][2] = y[i] * y[i];

        D2[i][0] = x[i];
        D2[i][1] = y[i];
        D2[i][2] = 1;
    }

    print_matrix_300x3(D1);

    transpose_300(D1, D1_t, 300, 3);
    print_matrix_3x300(D1_t);
 
    matmul_300(D1_t, D1, S1, 3, 300, 3);
    print_matrix_3x3(S1);

    // Compute the inverse of S3 and store the result in T
    // Compute the eigenvectors and eigenvalues of M and store the results in evec and eval

    // for(i = 0; i < DIM_1; i++) {
    //     cond[i] = 4 * evec[i][0] * evec[i][2] - evec[i][1] * evec[i][1];
    // }

    // Compute the values of a1 based on the condition
    // Compute the values of a based on a1 and T
}

