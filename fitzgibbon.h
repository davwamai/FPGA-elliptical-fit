#ifndef FITZGIBBON_H
#define FITZGIBBON_H    
void fitzgibbon(double* x, double* y, int size, double* a);

void inverse_3x3(double A[3][3], double B[3][3]);

void solve_cubic(double a, double b, double c, double d, double roots[3]);

void eigenvalues_3x3(double A[3][3], double eigenvalues[3]);

void eigenvectors_3x3(double A[3][3], double eigenvalues[3], double eigenvectors[3][3]);

void transpose_300(double A[300][3], double B[3][300], int rows, int cols);

void matmul_300(double A[3][300], double B[300][3], double C[3][3], int rowsA, int colsA, int colsB);

void matmul_3x3(double A[3][3], double B[3][3], double C[3][3]);

void print_matrix_300x3(double matrix[300][3]);
void print_matrix_3x300(double matrix[3][300]);
void print_matrix_3x3(double matrix[3][3]);
#endif