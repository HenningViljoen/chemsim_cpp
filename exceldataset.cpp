#include "exceldataset.h"

bool invalidChar(char c)
{
    return !(c>=0 && c <128);
}

void stripUnicode(std::string & str)
{
    str.erase(std::remove_if(str.begin(),str.end(), invalidChar), str.end());
}

exceldataset::exceldataset(std::string afilename)
{
    std::vector<double> tempdata;
    int temparraysize;
    int inbetweendatapoints = (int)ceil(global::PlantDataSampleT / global::SampleT);
    filename = afilename;



    /*for (int i = 0; i < rowCount; i++)
    {
                    for (int j = 0; j < colCount; j++)
                    {
                        tempdata[rowCount * j + i] = Convert.ToDouble(xlRange.Cells[i + 1, j + 1].Value2);
                        //MessageBox.Show(xlRange.Cells[i, j].Value2.ToString());
                    }
    }*/

    std::ifstream file(filename); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
    std::string astring;

    while ( file.good() )
    {
        std:getline( file, astring, '\r' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
        astring = ltrim(astring);
        stripUnicode(astring);
        //double dtemp = std::stod("334.48");
        //std::string::size_type sz;
        double d = std::stod(astring);
        tempdata.push_back(d);
        /*if (value.find('\n') != std::string::npos) {
                split_line(value, "\n", values);
        } else {
                values.push_back(value);
        }*/
    }

    /*temparraysize = rowCount * colCount;
    inbetweendatapoints = Convert.ToInt32(Math.Ceiling(global.PlantDataSampleT / global.SampleT));
    arraysize = temparraysize * inbetweendatapoints;
    tempdata = new double[temparraysize];*/
    //data = new double[arraysize];
    //std::vector<double> data;

    for (int i = 0; i < tempdata.size(); i++)
    {
         for (int j = 0; j < inbetweendatapoints; j++)
         {
              data.push_back(tempdata[i]);

         }
    }
}

void split_line(std::string& line, std::string delim, std::list<std::string>& values)
{
    size_t pos = 0;
    while ((pos = line.find(delim, (pos + 1))) != std::string::npos) {
        std::string p = line.substr(0, pos);
        values.push_back(p);
        line = line.substr(pos + 1);
    }

    if (!line.empty()) {
        values.push_back(line);
    }
}
