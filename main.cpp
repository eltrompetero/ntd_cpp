//
//  main.cpp
//  
//
//  Created by Eddie on 12/3/19.
//

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include "ntd.hpp"

// write output to csv for communication with Python
template <class T>
bool write(std::ofstream &myfile, std::vector<T> &v) {
    bool success = false;

    if (myfile.is_open()) {
        for (typename std::vector<T>::iterator it=v.begin(); it!=v.end(); ++it) {
            myfile << *it << ",";
        }
        myfile << "\n";
        success = true;
    } else {
        std::cout << "Unable to write file\n";
    }
    return success;
};

// read csv file from Python
// this is for the duration of simulations
bool read(std::string fname, std::vector<int> &t) {
    std::string line;
    int x;
    int i = 0;
    std::string num;
    bool success = false;
    std::ifstream myfile (fname);

    if (myfile.is_open()) {
        getline (myfile,line);
        while ( i<line.size() ) {
            num.clear();
            while (line[i]!=',' & i<line.size()) {
                num += line[i];
                i++;
            }
            i++;
            if (num.size()) {
                std::sscanf(num.data(), "%d", &x);
                t.push_back(x);
            }
        }
        myfile.close();
        success = true;
    } else {
        std::cout << "Unable to read file\n";
    }
    return success;
};

// command line args are
// input file name
// output files are the input file name with ".s" and ".r" suffixes
int main(int argc, const char * argv[]) {
    // parse command line args
    std::string instring;
    // read in input file name
    std::string ifile(argv[1]);
    // read in r
    instring = std::string (argv[2]);
    int r;
    std::sscanf(instring.data(), "%d", &r);
    // read in b
    instring = std::string (argv[3]);
    double b;
    std::sscanf(instring.data(), "%lf", &b);
    // read in theta
    instring = std::string (argv[4]);
    double theta;
    std::sscanf(instring.data(), "%lf", &theta);
    // read in gamma
    instring = std::string (argv[5]);
    double gmma;
    std::sscanf(instring.data(), "%lf", &gmma);
    
    // set output file name
    std::string osfname = ifile+".s";
    std::string orfname = ifile+".r";
    std::ofstream sfile (osfname.data());
    std::ofstream rfile (orfname.data());
    std::vector<int> t;  // to read in sim durations
    
    read(ifile, t);
    ConflictReportsTrajectory crt(r, b, theta, gmma);
    
    for (std::vector<int>::iterator it=t.begin(); it!=t.end(); ++it) {
        std::cout << "t" << *it << "\n";
        crt.grow(*it);
        write(sfile, crt.cumS);
        write(rfile, crt.radius);
    }
    sfile.close();
    rfile.close();
}//end main
