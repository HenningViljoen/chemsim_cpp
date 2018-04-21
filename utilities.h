#ifndef UTILITIES_H
#define UTILITIES_H

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <string>
#include <math.h>

#include "global.h"

#include "point.h"

class utilities
{
public:
    utilities();
    static double calcdirection(double deltay, double deltax);
    static double distance(point &p1, point &p2);
    static double distance(double deltax, double deltay);
};

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

static inline double round_precision(double num, double num_digits)
{
    return ceil(num*pow(10,num_digits))/pow(10,num_digits);
}

double dndt2fps(double adndt);

#endif // UTILITIES_H
