#include <math.h>
#include <cfloat>
#include "matrix.h"
#include "global.h"

matrix::matrix(int r, int c, double initialvalue)
{
    initmatrix(r, c, initialvalue);

    //nullthematrix();
}

matrix::matrix() //Constructor for populating the matrix by the user himself.
{
    m = std::vector<std::vector<double>> ();
}

matrix::matrix(matrix *matrixcopyfrom)
{
    copyfrom(matrixcopyfrom);
}

matrix::matrix(const matrix &m) //copy constructor
{
    copyfrom(&m);
}

void matrix::initmatrix(int r, int c, double initialvalue)
{
    m = std::vector<std::vector<double>> (r, std::vector<double>(c));
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            m[i][j] = initialvalue;
        }
    }
}

void matrix::copyfrom(const matrix *matrixcopyfrom)
{
    int copyfromr = matrixcopyfrom->m.size();
    int copyfromc;
    if (copyfromr == 0) { copyfromc = 0; }
    else { copyfromc = matrixcopyfrom->m[0].size(); }

    initmatrix(copyfromr, copyfromc);

    for (int i = 0; i < copyfromr; i++)
    {
        for (int j = 0; j < copyfromc; j++)
        {
            m[i][j] = matrixcopyfrom->m[i][j];
        }
    }
}

void matrix::nullthematrix()
{
    for (int r = 0; r < m.size(); r++)
    {
         for (int c = 0; c < m[r].size(); c++)
         {
             m[r][c] = 0;
         }
    }
}

void matrix::swoprowsinthematrix(int r1, int r2)
{
    std::vector<double> temprow = std::vector<double>(m[r1].size(), 0);
    for (int c = 0; c < m[r1].size(); c++)
    {
         temprow[c] = m[r1][c];
         m[r1][c] = m[r2][c];
         m[r2][c] = temprow[c];
    }
}

void matrix::deletelastrow()
{
    if (m.size() > 0) { m.pop_back(); }
}

void matrix::addarow()
{
    if (m.size() > 0)
    {
        m.push_back(std::vector<double>(m[0].size(), 0.0));
    }
}

void matrix::deletelastcolumn()
{
    if (m[0].size() > 0)
    {
        for (int r = 0; r < m.size(); r++)
        {
            m[r].pop_back();
        }
    }
}

void matrix::addacolumn()
{
    if (m.size() > 0)
    {
        for (int r = 0; r < m.size(); r++)
        {
            m[r].push_back(0.0);
        }
    }
}

void matrix::resize(int r, int c)
{
    while (m.size() > r)
    {
        deletelastrow();
    }
    while (m.size() < r)
    {
        addarow();
    }
    while (m[0].size() > c)
    {
        deletelastcolumn();
    }
    while (m[0].size() < c)
    {
        addacolumn();
    }
}

void matrix::idthematrix()
{
    nullthematrix();
    for (int i = 0; i < m.size(); i++)
    {
        if (i < m[0].size())
        {
             m[i][i] = 1.0;
        }
    }
}

bool matrix::issquare()
{
    return (m.size() == m[0].size());
}

void matrix::ludecomposition(matrix *l, matrix *u) //LU decomposition with the Doolittle algorithm
{
     int c = m[0].size();
     if (issquare() && l->issquare() && u->issquare() && l->m.size() == u->m.size() && l->m.size() == c)
     {
                l->idthematrix();
                u->idthematrix();
                for (int i = 0; i < c; i++)
                {
                    for (int j = i; j < c; j++)
                    {
                        u->m[i][j] = m[i][j];

                        for (int k = 0; k <= i; k++)
                        {
                            if (k != i)
                            {
                                u->m[i][j] = u->m[i][j] - l->m[i][k] * u->m[k][j];
                            }
                        }
                        u->m[i][j] /= l->m[i][i] + global::Epsilon;
                    }

                    for (int j = i + 1; j < c; j++)
                    {
                        l->m[j][i] = m[j][i];

                        for (int k = 0; k < i; k++)
                        {
                            if (k != i)
                            {
                                l->m[j][i] = l->m[j][i] - l->m[j][k] * u->m[k][i];
                            }
                        }
                        l->m[j][i] /= u->m[i][i] + global::Epsilon;
                    }
                }
      }
}



//the operator functions as defined below are not member functions of the class.
matrix &operator+(matrix &m1, matrix &m2)
{
    int mr = m1.m.size();
    int mc = m1.m[0].size();
    matrix mat(mr, mc);
    for (int r = 0; r < mr; r++)
    {
        for (int c = 0; c < mc; c++)
        {
            mat.m[r][c] = m1.m[r][c] + m2.m[r][c];
        }
    }
    return mat;
}

matrix &operator-(matrix &m1, matrix &m2)
{
     int mr = m1.m.size();
     int mc = m1.m[0].size();
     matrix mat(mr, mc);
     for (int r = 0; r < mr; r++)
     {
         for (int c = 0; c < mc; c++)
         {
             mat.m[r][c] = m1.m[r][c] - m2.m[r][c];
         }
     }
     return mat;
}

matrix &operator*(matrix &m1, matrix &m2)
{
     int m1r = m1.m.size();
     int m1c = m1.m[0].size();

     int m2r = m2.m.size();
     int m2c = m2.m[0].size();
     matrix mat(m1r, m2c);
     for (int r = 0; r < m1r; r++)
     {
          for (int c = 0; c < m2c; c++)
          {
               mat.m[r][c] = 0.0;
               for (int k = 0; k < m1c; k++)
               {
                    mat.m[r][c] += m1.m[r][k] * m2.m[k][c];
               }
          }
     }
     return mat;
}

matrix &operator*(double d, matrix &m2)
{
     int m2r = m2.m.size();
     int m2c = m2.m[0].size();
     matrix mat(m2r, m2c);
     for (int r = 0; r < m2r; r++)
     {
          for (int c = 0; c < m2c; c++)
          {
               mat.m[r][c] = d*m2.m[r][c];
          }
     }
     return mat;
}

matrix &transposevector(std::vector<double> v)
{
     int size = v.size();
     matrix mat(1, size);

     for (int r = 0; r < size; r++)
     {
          mat.m[0][r] = v[r];
     }

     return mat;
}

matrix &transpose(matrix &m)
{
     int mr = m.m.size();
     int mc = m.m[0].size();

     matrix mat(mc, mr);
     for (int r = 0; r < mr; r++)
     {
          for (int c = 0; c < mc; c++)
          {
               mat.m[c][r] = m.m[r][c];
          }
     }
     return mat;
}

void choleskyLDLT(matrix *a, matrix *l, matrix *d)
{
     int ar = a->m.size();
     int ac = a->m[0].size();
     int n = ar; //a must be square.
     double sumdl2 = 0;
     double sumdll = 0;

     //d = identitymatrix(n);
     l->initmatrix(n, n);
     d->initmatrix(n, n);
     matrix *c = new matrix(n, n);

     for (int j = 0; j < n; j++)
     {
           sumdl2 = 0;
           for (int s = 0; s <= j - 1; s++)
           {
               sumdl2 += d->m[s][s] * pow(l->m[j][s], 2);
           }
           c->m[j][j] = a->m[j][j] - sumdl2;
           d->m[j][j] = c->m[j][j];
           for (int i = j; i < n; i++)
           {
               sumdll = 0;
               for (int s = 0; s <= j - 1; s++)
               {
                    sumdll += d->m[s][s] * l->m[i][s]*l->m[j][s];
               }
               c->m[i][j] = a->m[i][j] - sumdll;
               l->m[i][j] = c->m[i][j] / (d->m[j][j] + global::Epsilon);
            }
      }
     delete c;

}

double euclideannorm(matrix *m) //the norm will be for the whole matrix, but will only really be correct if the matrix is a vector.
{
     int mr = m->m.size();
     int mc = m->m[0].size();
     double norm = 0;
     for (int r = 0; r < mr; r++)
     {
          for (int c = 0; c < mc; c++)
          {
                norm += pow(m->m[r][c],2);
          }
     }
     norm = sqrt(norm);
     return norm;
}

void solveLYequalsB(matrix *L, matrix *Y, matrix *B) //solves for the X vector the equation: LX = B, with L being a LU decomoposition matrix.
        // Dimensions have to be n x n for L
{
     if (L->issquare())
     {
          int Lr = L->m.size();
          int Lc = L->m[0].size();
          Y->nullthematrix();
          for (int r = 0; r < Lr; r++)
          {
               double sum = 0;
               for (int c = 0; c < r; c++)
               {
                    sum += L->m[r][c] * Y->m[c][0];
               }
               if (L->m[r][r] == 0) { Y->m[r][0] = DBL_MAX; }
               else { Y->m[r][0] = (B->m[r][0] - sum) / L->m[r][r]; }
          }
      }
}

void solveUXequalsY(matrix *U, matrix *X, matrix *Y)
{
      if (U->issquare())
      {
           int Ur = U->m.size();
           int Uc = U->m[0].size();
           X->nullthematrix();
           for (int r = Ur - 1; r >= 0; r--)
           {
               double sum = 0;
               for (int c = Uc - 1; c > r; c--)
               {
                    sum += U->m[r][c] * X->m[c][0];
               }
               if (U->m[r][r] == 0) { U->m[r][r] = global::Epsilon; }//X.m[r][0] = Double.NaN; }
               else { X->m[r][0] = (Y->m[r][0] - sum) / U->m[r][r]; }
            }
       }
}

void solveAXequalsB(matrix *A, matrix *X, matrix *B)
{
      int size = A->m.size();
      matrix lmatrix(size, size);
      matrix umatrix(size, size);
      matrix ymatrix(size, 1);

      A->ludecomposition(&lmatrix, &umatrix);
      //matrix tempm = lmatrix * umatrix;
      solveLYequalsB(&lmatrix, &ymatrix, B);
      //matrix tempm2 = lmatrix * ymatrix;
      solveUXequalsY(&umatrix, X, &ymatrix);

}

matrix &identitymatrix(int size) //Returns identity matrix.
{
      matrix mat(size, size);
      for (int r = 0; r < size; r++)
      {
           mat.m[r][r] = 1;
      }
      return mat;
}









