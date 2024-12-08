/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <atomic>
#include <algorithm>

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "mode.h"

//declaration of data struct
data d = chooseMode(0,0);

void RF_Thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& audio_read_offset, std::atomic<int>& RDS_read_offset);
void audio_Thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset);
void RDS_Thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset);

int main(int argc, char *argv[]){

	// validate input and set mode
	if (argc < 2) {
		std::cerr << "Operating in default mode 0" << std::endl;
	}
	else if (argc == 2) {
		if ((atoi(argv[1])) > 3) {
			std::cerr << "Wrong mode" << (atoi(argv[1])) << std::endl;
			exit(1);
		}
		else{
			d = chooseMode(atoi(argv[1]),atoi(argv[2]));
			d.channel = 0;
		}
	}
	else if (argc == 3) {
		if ((atoi(argv[1])) > 3) {
			std::cerr << "Wrong mode" << (atoi(argv[1])) << std::endl;
			exit(1);
		}
		else {
			d = chooseMode(atoi(argv[1]),atoi(argv[2]));
		}
		if(atoi(argv[2]) == 0 || atoi(argv[2]) == 1 || atoi(argv[2]) == 2) {
			if((d.mode==1 || d.mode==3) && atoi(argv[2]) == 2){
				std::cerr << "RDS is not supported for mode " << (atoi(argv[1])) << std::endl;
				exit(1);
			}
		}
		else{
			std::cerr << "Wrong channel " << (argv[2]) << std::endl;
		}
	}
	else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode> " << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;

		exit(1);
	}
	std::cerr << "Operating in mode " << d.mode << std::endl;
	
	std::vector<float> rf_queue(d.rf_block_size*d.q_ele);
	std::atomic<int> write_offset(0);
	std::atomic<int> audio_read_offset(0);
	std::atomic<int> RDS_read_offset(0);
	
	//Threads
	std::thread RF_produce = std::thread(RF_Thread, std::ref(rf_queue), std::ref(write_offset), std::ref(audio_read_offset), std::ref(RDS_read_offset));
	std::thread audio_consume = std::thread(audio_Thread, std::ref(rf_queue), std::ref(write_offset), std::ref(audio_read_offset));
	std::thread RDS_consume = std::thread(RDS_Thread, std::ref(rf_queue), std::ref(write_offset), std::ref(RDS_read_offset));
		
	RF_produce.join();
	audio_consume.join();
	RDS_consume.join();
	
	return 0;
}

void read_stdin_block(unsigned int num_samples, \
unsigned int block_id, \
std::vector<float> &block_data)
{
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		//auto normalized to +-1
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}

//Processes input data and outputs to audio thread and/or RDS thread
void RF_FrontEnd(std::vector<float>&fm_demod, std::vector<float>&iq_block, std::vector<float>&i_filt_state, std::vector<float>&q_filt_state, float&i_state, float&q_state){
	
	std::vector<float> rf_coeff;

	float i_der = 0;
	float q_der = 0;

	std::vector<float> yb_i;
	std::vector<float> yb_q;
	std::vector<float> xb_i;
	std::vector<float> xb_q;

	// coefficients for the front-end low-pass filter
	impulseResponseLPF(d.rf_FS,d.rf_FC,d.rf_taps,rf_coeff);
	
	//slice into I and Q blocks
	for(int i=0; i<iq_block.size(); i++){
		if (i%2 == 0){
			xb_i.push_back(iq_block[i]);
		} else{
			xb_q.push_back(iq_block[i]);
		}
	}

	//apply RF LPF to I and Q samples
	//if(d.mode == 2 || d.mode == 3){
		//convolveFIR_DSUS(yb_i,xb_i,rf_coeff,i_filt_state,d.rf_decim, d.expander);
		//convolveFIR_DSUS(yb_q,xb_q,rf_coeff,q_filt_state,d.rf_decim, d.expander);
	//}
	//else if (d.mode == 0 || d.mode == 1) {
		convolveFIR_DSUS(yb_i,xb_i,rf_coeff,i_filt_state,d.rf_decim, 1);
		convolveFIR_DSUS(yb_q,xb_q,rf_coeff,q_filt_state,d.rf_decim, 1);
	//}
	
	

	//demodulate, combine I and Q samples
	for(int i=0; i<yb_i.size(); i++){
		fm_demod.resize(yb_i.size(), 0.0);
			if(i==0){
				i_der = yb_i[i] - i_state;
				q_der = yb_q[i] - q_state;
			}
			else{
				i_der = yb_i[i] - yb_i[i-1];
				q_der = yb_q[i] - yb_q[i-1];
			}
			if((yb_i[i]!=0)&&(yb_q[i]!=0)){
				fm_demod[i] = (yb_i[i]*q_der - yb_q[i]*i_der) / (pow(yb_i[i],2) + pow(yb_q[i],2));
			}
			else if (yb_i[i] == 0 && yb_q[i] == 0 && i !=0){//avoid divide by zero error
				fm_demod[i] = fm_demod[i-1];
			}
			else{
				fm_demod[i] = 0;
			}
		}
		i_state = yb_i[yb_i.size()-1];//last element
		q_state = yb_q[yb_q.size()-1];//last element
}

void RF_Thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& audio_read_offset, std::atomic<int>& RDS_read_offset){
	std::vector<float> fm_demod;
	std::vector<float> i_filt_state(d.rf_taps, 0.0);
	std::vector<float> q_filt_state(d.rf_taps, 0.0);
	//std::cerr << "Entered RF THREAD" << std::endl;
	float i_state = 1;
	float q_state = 1;

	//Getting inputs
	for (unsigned int block_id=0; ; block_id++){
		std::vector<float> block_data(d.iq_block_size);
		read_stdin_block(d.iq_block_size, block_id, block_data);
		std::cerr << "read block " << block_id << "\n";
		if((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			exit(1);
		}

	RF_FrontEnd(fm_demod, block_data, i_filt_state, q_filt_state, i_state, q_state);

	//copy to queue and wait if necessary
	//while((write_offset.load() >= (audio_read_offset.load()+d.q_ele))  || (write_offset.load() >= (RDS_read_offset.load()+d.q_ele))) {
	while((write_offset.load() >= (audio_read_offset.load()+d.q_ele)) ) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    std::vector<float>::difference_type address_offset = (write_offset.load() % d.q_ele)*d.rf_block_size;
    std::copy_n(fm_demod.begin(), fm_demod.size(), rf_queue.begin()+address_offset);
    write_offset.fetch_add(1);
  }
}

void audio_Mono(std::vector<float>&processed_data, std::vector<float>&fm_demod, std::vector<float>&audio_coeff, std::vector<float>&filter_state, std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset,std::vector<float> &fm_demod_wait,std::vector<float> &fm_demod_state){
	while(write_offset.load() <= read_offset.load()) {
	std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
    std::vector<float>::difference_type address_offset = (read_offset.load() % d.q_ele)*d.rf_block_size;
    std::copy_n(rf_queue.begin()+address_offset, fm_demod.size(),fm_demod.begin());
    read_offset.fetch_add(1);
	//impulseResponseAPF(fm_demod_wait,fm_demod,fm_demod_state);
	//Extract mono data
    convolveFIR_DSUS(processed_data, fm_demod, audio_coeff, filter_state, d.audio_decim, d.expander);
}

void audio_Stereo(std::vector<float>&l_Channel, std::vector<float>&r_Channel,std::vector<float>&processed_data, std::vector<float>&fm_demod, std::vector<float>&pilot_coeff, std::vector<float>&stereo_coeff, std::vector<float>&stereo_filt_coeff, std::vector<float>&pilot_state,std::vector<float>&stereo_state, std::vector<float>&stereo_filt_state,pllState pll_state){
	std::vector<float> processed_pilot;
	std::vector<float> Carrier;
	std::vector<float> processed_stereo;
	std::vector<float> stereo_filt;

	//Carrier Recovery
	convolveFIR_DS(processed_pilot,fm_demod,pilot_coeff,pilot_state,1);
	fmPLL(processed_pilot,Carrier,pll_state,d.pilot_F,d.rf_IF,2.0,0.0,0.01);

	//Extract stereo data
	convolveFIR_DS(processed_stereo,fm_demod,stereo_coeff,stereo_state,1);

	//Mixing
	std::vector<float>mixing_data(processed_stereo.size(),0);
	std::transform(processed_stereo.begin(),processed_stereo.end(),Carrier.begin(),mixing_data.begin(),std::multiplies<float>());
		for (int i = 0;i<mixing_data.size();i++){
			mixing_data[i] *= 2;
		}

	//Rate conversion and final filter
	convolveFIR_DSUS(stereo_filt,mixing_data,stereo_filt_coeff,stereo_filt_state,d.audio_decim,d.expander);
	
	//Stereo combiner
	l_Channel.resize(stereo_filt.size(),0);
	r_Channel.resize(stereo_filt.size(),0);
	std::transform(processed_data.begin(),processed_data.end(),stereo_filt.begin(),l_Channel.begin(),std::plus<float>());
	std::transform(processed_data.begin(),processed_data.end(),stereo_filt.begin(),r_Channel.begin(),std::minus<float>());
}
	
void audio_Thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset){
  	//MONO
	std::vector<float> audio_coeff(d.audio_taps);
  	std::vector<float> demod_state((d.audio_taps-1)/2, 0);
	std::vector<float> processed_data;
	std::vector<float> filter_state(d.audio_taps*d.expander-1,0.0);
	std::vector<float> fm_demod(d.rf_block_size);
	std::vector<float> fm_demod_wait(fm_demod.size(),0);
	std::vector<float> fm_demod_state((d.rf_taps-1)/2,0);
	
	//STEREO
	std::vector<float> pilot_coeff;
	std::vector<float> stereo_coeff;
	std::vector<float> stereo_filt_coeff;
	std::vector<float> pilot_state(d.audio_taps-1,0.0);
	std::vector<float> stereo_state(d.audio_taps-1,0.0);
	std::vector<float> stereo_filt_state(d.audio_taps*d.expander-1,0.0);
	std::vector<float> l_Channel;
	std::vector<float> r_Channel;
	pllState pll_state;
  	pll_state.integrator = 0.0;
  	pll_state.phaseEst = 0.0;
  	pll_state.feedbackI = 1.0;
  	pll_state.feedbackQ = 0.0;
  	pll_state.nco0 = 1.0;
  	pll_state.trigOffset = 0.0;

	//STEREO
	impulseResponseBPF(pilot_coeff, d.p_FB, d.p_FE, d.rf_IF, d.audio_taps);
	impulseResponseBPF(stereo_coeff, d.audio_FB, d.audio_FE, d.rf_IF, d.audio_taps);
	impulseResponseLPF2(d.rf_IF*d.expander,d.audio_FC,d.audio_taps*d.expander,stereo_filt_coeff,d.expander);

	//MONO
  	impulseResponseLPF2(d.rf_IF*d.expander,d.audio_FC,d.audio_taps*d.expander,audio_coeff,d.expander);

	//MONO or STEREO
  	while(1){

	audio_Mono(processed_data,fm_demod,audio_coeff,filter_state,rf_queue,write_offset,read_offset,fm_demod_wait,fm_demod_state);

	if(d.channel==1){
		audio_Stereo(l_Channel,r_Channel,processed_data,fm_demod,pilot_coeff,stereo_coeff,stereo_filt_coeff,pilot_state,stereo_state,stereo_filt_state,pll_state);
	}
    if (d.channel == 0) {
    	std::vector<short int> sample(processed_data.size());
    	for (int k=0; k < processed_data.size(); k++) {
    		if (std::isnan(processed_data[k])) sample[k] = 0;
    		else sample[k] = static_cast<short int> (processed_data[k] * 16384);
      }
      fwrite(&sample[0], sizeof(short int), sample.size(), stdout);
	}
	else if (d.channel == 1){
		std::vector<short int> sample(l_Channel.size()*2);
    	for (int k=0; k < sample.size(); k++) {
			if (k%2 == 0){
        		if (std::isnan(l_Channel[k/2])) sample[k] = 0;
        		else sample[k] = static_cast<short int> (l_Channel[k/2] * 16384);
      		}
			else{
				if (std::isnan(r_Channel[(k-1)/2])) sample[k] = 0;
        		else sample[k] = static_cast<short int> (r_Channel[(k-1)/2] * 16384);
	  		}
		}
		fwrite(&sample[0], sizeof(short int), sample.size(), stdout);
	}
  }
}

void RDS(std::vector<float>&rf_queue, std::vector <float> RDS_state, std::vector<float>&RDS_pilot_state, pllState&RDS_pilot_pll_state,  std::vector<float>&RDS_resampled_state, std::vector<float>&RDS_RRC_state, std::vector <int> &RDS_CDR_state, bool& RDS_diff_state, std::vector<float>&RDS_coeff, std::vector<float>&RDS_pilot_coeff, std::vector<float>&RDS_demod_coef, std::vector<float>&RRC_coeff, std::atomic<int>& write_offset, std::atomic<int>& read_offset){
	
	std::vector<float> fm_demod(d.rf_block_size);

	while(write_offset.load() <= read_offset.load()){
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
	std::vector<float>::difference_type address_offset = (read_offset.load() % d.q_ele)*d.rf_block_size;
	std::copy_n(rf_queue.begin()+address_offset, fm_demod.size(),fm_demod.begin());
	read_offset.fetch_add(1);

	std::vector <float> RDS_data;
	std::vector <float> RDS_delayed;		
	std::vector <float> RDS_delayed_state;
	std::vector <float> RDS_square;
	std::vector <float> RDS_pilot;
	std::vector <float> RDS_pilot_pll;
	std::vector<float>mixing_data;
	std::vector<float>RDS_resampled;
	std::vector<float>RDS_RRC;
	std::vector<bool>RDS_CDR;
	std::vector<bool> RDS_diff_decode;

	//BPF extract RDS channel
	convolveFIR_DSUS(RDS_data, fm_demod, RDS_coeff, RDS_state, 1, 1);

	//APF to sync
	RDS_delayed.resize(RDS_data.size());
	impulseResponseAPF(RDS_delayed_state, RDS_data, RDS_delayed_state);

	//Squaring Nonlinearity to double frequency
	RDS_square.resize(RDS_data.size());
	std::transform(RDS_data.begin(), RDS_data.end(), RDS_data.begin(), RDS_square.begin(), std::multiplies<float>());

	//BPF extract pilot
	convolveFIR_DSUS(RDS_pilot, RDS_square, RDS_pilot_coeff, RDS_pilot_state, 1, 1);

	//PLL to synchronize phase
	fmPLL(RDS_pilot, RDS_pilot_pll, RDS_pilot_pll_state, d.RDS_pilot_F, d.rf_IF, 0.5, 0.0, 0.01);

	//Mixing
	mixing_data.resize(RDS_data.size(),0);
	std::transform(RDS_data.begin(),RDS_data.end(),RDS_pilot_pll.begin(),mixing_data.begin(),std::multiplies<float>());
	for (int i = 0;i<mixing_data.size();i++){
		mixing_data[i] *= 2;
	}

	//LPF and Rational Resampler
	convolveFIR_DSUS(RDS_resampled, mixing_data, RDS_demod_coef, RDS_resampled_state, d.RDS_DS, d.RDS_US);

	//Root Raised Cosine Filter
	convolveFIR_DSUS(RDS_RRC, RDS_resampled, RRC_coeff, RDS_RRC_state, 1, 1);

	//Carrier and Data Recovery
	//find_sampling_pos();
	//CDR(RDS_CDR, RDS_RRC, RDS_CDR_state, 0, 0, d.SPS);

	//Differential Decoder
	//diff_decode(RDS_diff_decode, RDS_CDR, RDS_diff_state);

	//Frame Sync
	//frame_sync();

}

void RDS_Thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset){
/*
	std::vector<float> RDS_coeff;
	std::vector<float> RDS_pilot_coeff;
	std::vector<float> RDS_demod_coeff;
	std::vector<float> RCC_coeff;

	std::vector<float> RDS_state;
	std::vector<float> RDS_pilot_state;
	pllState RDS_pilot_pll_state;
	std::vector<float> RDS_resampled_state;
	std::vector<float> RDS_RRC_state;
	std::vector<int> RDS_CDR_state;
	bool RDS_diff_state;

	if(d.mode == 0 || d.mode == 2){

		std::cerr<<"Run RDS"<<std::endl;

		impulseResponseBPF(RDS_coeff, d.RDS_FB, d.RDS_FE, d.rf_IF, d.RDS_taps);
		impulseResponseBPF(RDS_pilot_coeff, d.RDS_p_FB, d.RDS_p_FE, d.rf_IF, d.RDS_taps);
		impulseResponseLPF2(d.RDS_FS_demod, d.RDS_FC, d.RDS_taps_demod, RDS_demod_coeff, d.RDS_US);
		RRC(d.SPS*2375, d.RDS_taps,RCC_coeff);

		while(1){
			RDS(rf_queue, RDS_state, RDS_pilot_state, RDS_pilot_pll_state, RDS_resampled_state, RDS_RRC_state, RDS_CDR_state, RDS_diff_state, RDS_coeff, RDS_pilot_coeff, RDS_demod_coeff, RCC_coeff, write_offset, read_offset);
		}
	}*/
}
