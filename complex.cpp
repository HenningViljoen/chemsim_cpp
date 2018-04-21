#include "complex.h"

complex::complex(double aa, double ab)
{
     a = aa;
     b = ab;
     //I = new complex(0, 1);
}

complex::complex(const complex &c)
{
    a = c.a;
    b = c.b;
}

complex::complex()
{
    a = 0;
    b = 0;
}

complex::~complex()
{
    //if (I != nullptr) {delete I;}
}

complex &complex::operator=(const complex &c)
{
    complex comp(c.a, c.b);
    return comp;
}

complex &complex::power(complex &x, double y, int n) //the specific root that is looked for will be specified here.  Could be a number between 0 and 1/y.
        //This is assuming that y is smaller than 1.
{
    double r = rad(x);
    double al = alpha(x);
    double newr = pow(r, y);
    double newalpha;
    if (y < 1)
    {
         newalpha = al * y + 2.0 * M_PI * y * n;
    }
    else
    {
         newalpha = al * y;
    }
    complex comp(newr * cos(newalpha), newr * sin(newalpha));
    return comp;
}

complex &complex::power(double x, double y) //This function will have to be renamed when the code base where it is called since it was
                                                //renamed here.
{
    if (x < 0)
    {
        complex comp(0, pow(-x, y));
        return comp;
    }
    else
    {
        complex comp(pow(x, y), 0);
        return comp;
    }
}




//static and operator methods and functions
complex &c(double aa)
{
    complex comp(aa, 0);
    return comp;
}

double rad(complex &c)
{
     return sqrt(c.a * c.a + c.b * c.b);
}

double alpha(complex &c)
{
     if (c.b == 0 && c.a < 0) { return -M_PI; }
     else if (c.a == 0 && c.b > 0) { return M_PI / 2.0; }
     else if (c.a == 0 && c.b < 0) { return -M_PI / 2.0; }
     else { return atan(c.b / c.a); }
}

complex &operator+(complex &c1, complex &c2)
{
    complex comp(c1.a + c2.a, c1.b + c2.b);
    return comp;
}

complex &operator+(double d, complex &c)
{
    complex comp(d + c.a, c.b);
    return comp;
}

complex &operator+(complex &c, double d)
{
    complex comp(d + c.a, c.b);
    return comp;
}

complex &operator-(complex &c1, complex &c2)
{
    complex comp(c1.a - c2.a, c1.b - c2.b);
    return comp;
}

complex &operator-(double d, complex &c)
{
    complex comp(d - c.a, -c.b);
    return comp;
}

complex &operator-(complex &c, double d)
{
    complex comp(c.a - d, c.b);
    return comp;
}

complex &operator*(complex &c1, complex &c2)
{
    complex comp(c1.a * c2.a - c1.b * c2.b, c1.b * c2.a + c1.a * c2.b);
    return comp;
    //(a+bi) (c+di) = (ac-bd) + (bc+ad)i.
}

complex &operator*(double d, complex &c)
{
    complex comp(d * c.a, d * c.b);
    return comp;
}

complex &operator*(complex &c, double d)
{
    complex comp(d * c.a, d * c.b);
    return comp;
}





complex &conj(complex &c)
{
    complex comp(c.a, -c.b);
    return comp;
}

double realconj(complex &c)
{
    return c.a * c.a + c.b * c.b;
}

complex &operator/(complex &c1, complex &c2)
{
     complex num = c1 * conj(c2);
     double den;

     den = realconj(c2);
     complex comp(num.a / den, num.b / den);
     return comp;
}

complex &operator/(complex &c, double d)
{
    complex comp(c.a / d, c.b / d);
    return comp;
}

complex &operator/(double d, complex &c)
{
    complex num = d * conj(c);
    double den;

    den = realconj(c);
    complex comp(num.a / den, num.b / den);
    return comp;
}








