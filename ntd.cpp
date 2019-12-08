//
//  ntd.cpp
//  ntd
//
//  Created by Eddie on 12/3/19.
//  Copyright © 2019 Santa Fe Institute. All rights reserved.
//

#include "ntd.hpp"

int max(std::vector<int> x) {
    int mx = x[0];
    for (int i : x) {
        if (i>mx) {
            mx = i;
        }
    }
    return mx;
};

int max_len(std::vector<Branch> &x) {
    int mx = x[0].pos + x[0].ancestralLen;
    for (std::vector<Branch>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ((it->pos + it->ancestralLen)>mx) {
            mx = it->pos + it->ancestralLen;
        }
    }
    return mx;
};


Branch::Branch() {};

Branch::Branch(std::string newlabel, int newlength, int newancestrallength) {
    label = newlabel;
    len = newlength;
    ancestralLen = newancestrallength;
    pos = 0;
};

bool Branch::grow(int step_size) {
    pos += step_size;
    if (pos>len) {
        pos = len;
        return false;
    }
    return true;
};

void NTD::print() {
    std::cout << "Growing branches\n";
    for (Branch br : growingBranches) {
        std::cout << br.label << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Dead branches\n";
    for (Branch br : deadBranches) {
        std::cout << br.label << "\n";
    }
};


NTD::NTD() {};

NTD::NTD(int newr, double newb, int rng_seed) {
    r = newr;
    b = newb;
    rngSeed = rng_seed;
};

void NTD::clear() {
    growingBranches.clear();
    deadBranches.clear();
    radius.clear();
};

void NTD::grow_random(int n_steps, int record_every, double mx_rand_factor) {
    // Parse args
    assert (n_steps>0);
    assert (record_every>=1);
    if (mx_rand_factor==-1) {
        mx_rand_factor = b;
    }
    assert (mx_rand_factor>0);
    
    int n, randix, i;
    int counter = 1;
    Branch nb, br;
    Branch *gb;
    std::vector<Branch> newBranches;
    std::uniform_real_distribution<double> factorrng(1/mx_rand_factor, mx_rand_factor);
    std::uniform_int_distribution<int> intrng;
    std::ostringstream stringStream;
    
    this->clear();
    
    for (int i=0; i<r; ++i) {
        stringStream << i;
        growingBranches.push_back(Branch(stringStream.str(), b));
        stringStream.str("");
    }
    
    while (counter<n_steps) {
        // iterate through the current set of branches that are growing without considering
        // new branches that are added in this loop
        n = growingBranches.size();
        i = 0;  // counter for number in current generation remaining
        while (i<n) {
            intrng = std::uniform_int_distribution<int> (0,n-1);
            randix = intrng(rd);
            gb = &(growingBranches[randix]);  // have to access this by ref to modify instance
                                              // stored in vector
            if (!gb->grow()) {
                // shouldn't these be extended by 1 when they're initialized
                // difficultly arises from the fact that some may start with 0 length
                for (int j=0; j<r; ++j) {
                    stringStream << gb->label << j;
                    nb = Branch(stringStream.str(),
                                int(pow(b, gb->label.size()+1)*factorrng(rd)),
                                gb->len+gb->ancestralLen);
                    stringStream.str("");
                    growingBranches.push_back(nb);
                    // vector heap allocation might have changed
                    gb = &(growingBranches[randix]);
                }
                deadBranches.push_back(growingBranches[randix]);
                growingBranches.erase(growingBranches.begin()+randix);
                n--;
                i--;
            }//end if
            if ((counter%record_every)==0) {
                radius.push_back(max_len(growingBranches));
            }//end if
            counter++;
            i++;
        }//end while
    }//end while
};

ConflictReportsTrajectory::ConflictReportsTrajectory() {};

ConflictReportsTrajectory::ConflictReportsTrajectory(int newr,
                                                     double newb,
                                                     double newtheta,
                                                     double newgammas,
                                                     double newgammaf,
                                                     int rng_seed) {
    assert (newr>1);
    assert (newb>1);
    assert (newtheta>=0.0);
    assert (newgammas>=-1.0);
    assert (newgammaf>=-1.0);
    
    r = newr;
    b = newb;
    theta = newtheta;
    gammas = newgammas;
    gammaf = newgammaf;
    rngSeed = rng_seed;
};

void ConflictReportsTrajectory::grow(int n_steps, int record_every, double mx_rand_coeff) {
    assert (n_steps>0);
    assert (record_every>=1);
    cumS = std::vector<int> (n_steps, 0);
    cumF = std::vector<int> (n_steps, 0);
    if (mx_rand_coeff==-1) {
        mx_rand_coeff = b;
    }
    assert (mx_rand_coeff>0);
    
    int randix;
    int counter = 0;
    std::uniform_real_distribution<double> factorrng(1/mx_rand_coeff, mx_rand_coeff);
    std::uniform_int_distribution<int> randixrng;
    std::ostringstream stringStream;
    Branch *gb;
    Branch nb;
    
    this->clear();
    
    // create first generation emanating from center
    for (int i=0; i<r; ++i) {
        stringStream << i;
        growingBranches.push_back(Branch(stringStream.str(), b));
        stringStream.str("");
    }
    
    while (counter<n_steps) {
        randixrng = std::uniform_int_distribution<int> (0, int(growingBranches.size())-1);
        randix = randixrng(rd);
        gb = &(growingBranches[randix]);
        // if branch fails to grow, then spawn children and try again
        if (!gb->grow()) {
            // spawn new branches
            for (int j=0; j<r; ++j) {
                stringStream << gb->label << j;
                nb = Branch(stringStream.str(),
                            int(pow(b,gb->label.size()+1) * factorrng(rd)),
                            gb->len+gb->ancestralLen);
                nb.grow();
                growingBranches.push_back(nb);
                stringStream.str("");
                gb = &(growingBranches[randix]);
            }
            deadBranches.push_back(growingBranches[randix]);
            growingBranches.erase(growingBranches.begin()+randix);
        } else {
            for (int j=counter; j<n_steps; ++j) {
                cumS[j] += int(pow(j-counter, 1+gammas) * pow(counter+1, -theta));
                cumS[j] += 1;
                cumF[j] += int(pow(j-counter, 1+gammaf) * pow(counter+1, -theta));
            }
            if ((counter%record_every)==0) {
                radius.push_back(max_len(growingBranches));
            }//end if
            counter++;
        }//end while
    }//end while
};//end grow

