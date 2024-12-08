/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


// function to compute the impulse response for a LPF, "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h){
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);
	double norm_Cutoff = Fc/(Fs/2);

	for(int i=0; i<num_taps; i++){
		float c =0;

		if (i == (num_taps-1)/2){
			c = norm_Cutoff;
		}
		else{
			c = norm_Cutoff*sin(PI*norm_Cutoff*(i-(num_taps-1)/2))/(PI*norm_Cutoff*(i-(num_taps-1)/2));
		}
		h[i] = c*pow(sin(i*PI/num_taps), 2);
	}
}

// function to compute the impulse response for a LPF, "h" based on the sinc function and scaled by a factor of U
void impulseResponseLPF2(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h, int U){
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);
	double norm_Cutoff = Fc/(Fs/2);

	for(int i=0; i<num_taps; i++){
		float c =0;
		if (i == (num_taps-1)/2){
			c = norm_Cutoff;
		}
		else{
			c = norm_Cutoff*sin(PI*norm_Cutoff*(i-(num_taps-1)/2))/(PI*norm_Cutoff*(i-(num_taps-1)/2));
		}
		h[i] = U* (c*pow(sin(i*PI/num_taps), 2));
	}
}

//APF used to add delay to match BPF delay
void impulseResponseAPF(std::vector<float> &y, std::vector<float> &x, std::vector<float> &state){

	std::copy(state.begin(), state.end(), y.begin());
	std::copy(x.begin(),x.end()-state.size(),y.begin()+state.size());
	std::copy(x.end()-state.size(),x.end(),state.begin());
}

//Basic convolve
void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &state){

	for(int n=0; n<x.size(); n++){
		y[n] = 0;
		for(int k=0; k<h.size(); k++){
			if (n<k){
				y[n] += state[state.size()+n-k]*h[k];
			}
			else{
				y[n] += x[n-k]*h[k];
			}
		}
	}
	y.resize(x.size());
	for(int a=0; a<state.size(); a++){
		state[a] = x[x.size()-state.size()+a];
	}
}

//Fast convolve using DS and US factors to only compute required outputs
void convolveFIR_DSUS(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &state, const int DS, const int US){

	y.resize(x.size()*US/DS);

	for(int n=0; n<y.size(); n++){
		y[n] = 0;
		int phase = (DS*n)%US;
		int index;
		for(int k=phase; k<h.size(); k+=US){
			index = (int)(n*DS-k)/US;
			if (index<0){
				y[n] += state[state.size()+index]*h[k];
			}
			else{
				y[n] += x[index]*h[k];
			}
		}
	}
	//y.resize(x.size()*US/DS);
	for(int a=0; a<state.size(); a++){
		state[a] = x[x.size()-state.size()+a];
	}
}


//Fast convolve using DS factor to only compute required outputs
void convolveFIR_DS(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &state, const int DS){

	y.resize(x.size()/DS);
	
	for(int n=0; n<x.size(); n++){
		y[n] = 0;
		int index;
		for(int k=0; k<h.size(); k++){
			if ((index<k)){
				y[n] += state[state.size()+index]*h[k];
			}
			else{
				y[n] += x[index]*h[k];
			}
		}
	}

	for(int a=0; a<state.size(); a++){
		state[a] = x[x.size()-state.size()+a];
	}
}

//BPF coefficients to remove frequencies not within the FB to FE range
void impulseResponseBPF (std::vector <float> &h, float FB, float FE, float FS,  unsigned short int N_taps){
	
	h.resize(N_taps,0.0);
	double norm_C = ((FE+FB)/2)/(FS/2);
	double norm_P = (FE-FB)/(FS/2);

	for (int i = 0; i< N_taps;i++){
		if (i == (N_taps - 1)/2){
			h[i] = norm_P;
		}
		else{
			h[i] = norm_P*((sin(PI*(norm_P/2)*(i-(N_taps-1)/2)))/(PI*(norm_P/2)*(i-(N_taps-1)/2)));
		}
		h[i] = h[i]*cos(i*PI*norm_C);
		h[i] = h[i]*pow(sin(i*PI/N_taps),2);
	}
}

//PLL used to synchronize phase
void fmPLL(std::vector<float>&pllIn, std::vector<float>&ncoOut, pllState &state,float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth){

	float Cp = 2.666;//scale factors
	float Ci = 3.555;

	float Kp = normBandwidth*Cp;//proportional gain
	float Ki = normBandwidth*Ci;//integrator gain

	ncoOut.resize(pllIn.size()+1, 0.0);//output array for the NCO

	//initialize internal state
	float integrator = 0;
	float phaseEst = 0;
	float feedbackI = 0;
	float feedbackQ = 0;
	ncoOut[0] = 0;
	float trigOffset = 0;

	for(int k=0; k<pllIn.size(); k++){

		//phase detector
		float errorI = pllIn[k]*(+feedbackI);//complex conjugate of the
		float errorQ = pllIn[k]*(-feedbackQ);//feedback complex exponential

		//four-quadrant arctangent discriminator for phase error detection
		float errorD = std::atan2(errorQ, errorI);

		//loop filter
		integrator += Ki*errorD;

		//loop filter
		phaseEst += Kp*errorD + integrator;

		//internal oscillator
		trigOffset += 1;
		float trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);

		state.integrator = integrator;
		state.phaseEst = phaseEst;
		state.feedbackI = feedbackI;
		state.feedbackQ = feedbackQ;
		state.nco0 = ncoOut.back();
		state.trigOffset = trigOffset;

		//for stereo only the in-phase NCO component should be returned
		//for block processing you should also return the state

		//for RDS add also the quadrature NCO component to the output
	}
}

//Root Raised Cosine coefficients used to combat inter-symbol interference
void RRC(int FS, int N_taps, std::vector<float> &impulseResponseRRC){

/*
	Root raised cosine (RRC) filter

	Fs  		sampling rate at the output of the resampler in the RDS path
				sampling rate must be an integer multipler of 2375
				this integer multiple is the number of samples per symbol

	N_taps  	number of filter taps
*/

	float T_symbol = 1/2375.0;
	float beta = 0.9;

	impulseResponseRRC.resize(N_taps, 0.0);

	for(int k=0; k<N_taps; k++){
		float t = float(k-N_taps/2)/FS;
		if(t == 0.0){
			impulseResponseRRC[k] = 1.0 + beta*((4/PI)-1);
		}
		else if((t == -T_symbol/(4*beta)) || (t == T_symbol/(4*beta))){
			impulseResponseRRC[k] = (beta/sqrt(2))*(((1+2/PI)*(sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
		}
		else{
			impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_symbol) + 4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/(PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
		}
	}
}

//Converts BPSK into digital signal by synchronizing with local min, max to digitize SPS data
void CDR (std::vector <bool> &output, std::vector <float>&preCDR, std::vector <int> &state, int &starting,  int &flag, int SPS)
{
	std::vector<int> temp;
	std::vector<int> indication;
	
	output = {};
	
	int sizepreCDR = preCDR.size();
	int sizestate = state.size();
	int sizeind = indication.size();
	
	indication.resize(sizepreCDR + sizestate);
	fill(indication.begin(), indication.end(), 0);
	
	copy(state.begin(), state.end(), indication.begin());
	state = {};
	
	output.clear();
	
	for (int i = starting; i >sizepreCDR; i+= SPS)
	{
		if (preCDR[i]> 0)
		{
			indication.push_back(1);
		}
		else
		{
			indication.push_back(-1);
		}
	}
	
	if (flag) 
	{
    bool has_left_symbol = (indication.size() % 2 != 0);

    for (auto i = indication.begin(); i < indication.end() - 1; i += 2) 
    {
        int a = *i;
        int b = *(i + 1);

        if (a == b) 
        {
            output.clear();
            for (auto j = indication.begin() + 1; j < indication.end() - 1; j += 2) 
            {
                int c = *j;
                int value = (c == 1) ? 1 : 0;
                output.push_back(value);
            }
            if (has_left_symbol) 
            {
                state.push_back(indication.back());
                has_left_symbol = false;
            }
            // return output, state, start
        } 
        else 
        {
            int value = (a == 1) ? 1 : 0;
            output.push_back(value);
        }
    }

		if (has_left_symbol) 
		{
			state.push_back(indication.back());
		}
	}
	
	// else, the block is considered normal and henceforth is paired directly. 
	
	else{
		for ( auto i = indication.begin(); i< indication.end() -1 ; i +=2){
			int value = (*i == 2) ? 1: 0;
			output.push_back(value);
		}
		if (indication.size() %2 !=0){
			 state.push_back(indication[sizeind-1]);
		}
	}
		
}

//Implement XOR on every Manchester decoded bit
void diff_decode(std::vector<bool>& decoded_data, std::vector<bool>& encoded_data, bool diff_state) {
    decoded_data.resize(encoded_data.size());

    for (size_t i = 0; i < encoded_data.size(); ++i) {
        decoded_data[i] = (diff_state ^ encoded_data[i]);
        diff_state = (diff_state ^ decoded_data[i]);
    }
}

int find_sampling_pos(std::vector <float> data_block, int SPS){

    float mid = 0.0f;
    int i = 0, j = 0, pos = 0;
    std::vector<float> chunk;

    int offset = SPS / 2;

    std::vector<float> array = std::vector<float>(SPS, 0);

    for (int i = 0; i < SPS; i++) 
    {
    int j = i;
    float error = 0.0f;
    
    std::vector<float> chunk(SPS + 1, 0.0f);
    while (j + SPS + 1 < data_block.size()) 
    {
        copy(data_block.begin() + j, data_block.begin() + j + SPS + 1, chunk.begin());
        mid = chunk[offset];
        if (chunk[0] * chunk.back() >= 0) {
            error += std::abs(mid - (chunk[0] + chunk.back()) / 2);
        } else {
            error += std::abs(mid - 0.0f);
        }
        j += SPS;
    }
    array.push_back(error);
    }
    auto position = min_element(array.begin(), array.end()) - array.begin();
    return position;
}
	
