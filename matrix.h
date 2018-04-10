#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class matrix
{
public:
    std::vector<std::vector<double>> m;

    matrix(int r, int c, double initialvalue = 0.0);
    matrix();
    matrix(matrix *matrixcopyfrom);
    matrix(const matrix &m); //copy constructor

    void initmatrix(int r, int c, double initialvalue = 0.0);
    void copyfrom(const matrix *matrixcopyfrom);
    void nullthematrix();
    void swoprowsinthematrix(int r1, int r2);
    void deletelastrow();
    void addarow();
    void deletelastcolumn();
    void addacolumn();
    void resize(int r, int c);
    void idthematrix();
    bool issquare();
    void ludecomposition(matrix *l, matrix *u);

};

matrix &operator+(matrix &m1, matrix &m2);
matrix &operator-(matrix &m1, matrix &m2);
matrix &operator*(matrix &m1, matrix &m2);
matrix &operator*(double d, matrix &m2);
matrix &transposevector(std::vector<double> v);
matrix &transpose(matrix &m);
void choleskyLDLT(matrix *a, matrix *l, matrix *d);
double euclideannorm(matrix *m);
void solveLYequalsB(matrix *L, matrix *Y, matrix *B);
void solveUXequalsY(matrix *U, matrix *X, matrix *Y);
matrix &identitymatrix(int size);

#endif // MATRIX_H
