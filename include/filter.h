/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <cmath>

//struct used for PLL state values
typedef struct pllState{
    float integrator;
    float phaseEst;
    float feedbackI;
    float feedbackQ;
    float nco0;
    float trigOffset;
}pllState;

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void impulseResponseLPF2(float, float, unsigned short int, std::vector<float> &, int);
void impulseResponseAPF(std::vector<float> &, std::vector<float> &, std::vector<float> &);
void convolveFIR(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &);
void convolveFIR_DS(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, const int);
void convolveFIR_DSUS(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, const int, const int);
void impulseResponseBPF (std::vector <float> &, float, float, float,  unsigned short int);
void fmPLL(std::vector<float>&, std::vector<float>&, pllState &, float, float, float, float, float);
void RRC(int , int , std::vector<float> &);
void CDR (std::vector <bool> &, std::vector <float>&, std::vector <int> &, int &,  int &, int);
void diff_decode(std::vector<bool>&, std::vector<bool>&, bool);
int find_sampling_pos(std::vector <float>, int);




#endif // DY4_FILTER_H
