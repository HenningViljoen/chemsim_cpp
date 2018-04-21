#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>
#include "global.h"

class complex
{
public:
    double a; //real component.
    double b; //imaginary component.
    //complex *I = nullptr; //; {0, 1};

    complex(double aa, double ab);
    complex(const complex &c);
    complex();
    ~complex();

    complex &operator=(const complex &c);
    static complex &power(complex &x, double y, int n = 0);
    static complex &power(double x, double y);
};

complex &c(double aa);
double rad(complex &c);
double alpha(complex &c);
complex &operator+(complex &c1, complex &c2);
complex &operator+(double d, complex &c);
complex &operator+(complex &c, double d);
complex &operator-(complex &c1, complex &c2);
complex &operator-(double d, complex &c);
complex &operator-(complex &c, double d);
complex &operator*(complex &c1, complex &c2);
complex &operator*(double d, complex &c);
complex &operator*(complex &c, double d);
complex &conj(complex &c);
double realconj(complex &c);
complex &operator/(complex &c, double d);
complex &operator/(double d, complex &c);


#endif // COMPLEX_H









