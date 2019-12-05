//
//  ntd.hpp
//  ntd
//
//  Created by Eddie on 12/3/19.
//  Copyright © 2019 Santa Fe Institute. All rights reserved.
//

#ifndef ntd_hpp
#define ntd_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <random>
#include <array>
#include <vector>
#include <math.h>
#include <sstream>
#include <assert.h>

#endif /* ntd_hpp */

class Branch{
public:
    std::string label;
    int len;
    int pos;
    int ancestralLen;
    
    Branch();
    Branch(std::string, int, int=0);
    bool grow(int=1);
};

class NTD{
public:
    double b;
    int r;
    int rngSeed;
    std::vector<Branch> growingBranches;
    std::vector<Branch> deadBranches;
    std::vector<int> cumS;
    std::vector<int> radius;
    std::mt19937 rd;
    
    NTD();
    NTD(int, double, int=-1);
    void clear();
    void grow_random(int, int=10, double mx_rand_factor=-1.0);
    void print();
};

class ConflictReportsTrajectory : public NTD {
public:
    double theta, gamma;
    
    ConflictReportsTrajectory();
    ConflictReportsTrajectory(int, double, double, double, int=-1);
    void grow(int, int=10, double=-1.0);
    void sample(int);
};
