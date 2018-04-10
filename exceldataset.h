#ifndef EXCELDATASET_H
#define EXCELDATASET_H

#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "global.h"


#include "utilities.h"

class exceldataset
{
public:
    std::string filename;
    std::vector<double> data;
    std::list<std::string> values;
    int arraysize;

    exceldataset(std::string afilename);
};

void split_line(std::string& line, std::string delim, std::list<std::string>& values);

#endif // EXCELDATASET_H
