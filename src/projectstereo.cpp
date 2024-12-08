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

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

void readStdinBlockData(unsigned int num_samples, \
unsigned int block_id,
std::vector<float> &block_data)
{
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		//auto normalized to +-1
		block_data[k] = float(((unisgned char)raw_data[k]-128)/128.0);
	}
}

int main()
{

	/*

	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
	std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	// the function is already provided in fourier.cpp
	//FFT_optimized(slice_data,Xf);
	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp
	//computeVectorMagnitude(Xf,Xmag);
	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//

*/

	int rf_FS = 2.4e6;
	int rf_FC = 100e3;
	int rf_taps = 151;
	int rf_decim = 10;

	int audio_FS = 48e3;
	int audio_FC = 16e3;
	int audio_taps = 151;
	int audio_decim = 5;

	int iq_block_size = rf_decim*2*audio_decim*1024;
	int block_size = rf_decim*audio_decim*1024;


for(unsigned int block_id=0; ; block_id++){
	std::vector<float> block_data(block_size);
	readStdinBlockData(block_size, block_id, block_data);
	if ((std::cin.rdstate()) != 0){
		std:cerr <<"End of input stream reached" << std::end1;
		exit(1);
	}
	std::cerr << "read block" << block_id << std:endl;
}

/*
	// read the raw IQ data from the recorded file
  // IQ data is assumed to be in 8-bits unsigned (and interleaved)
  std::string in_fname = "../data/iq_samples.raw";
  std::ifstream in_file(in_fname, std::ios::binary);

if(!in_file) {
	std::cout << "File " << in_fname << " not found ... exiting\n";
	exit(1);
} else {
	std::cout << "Reading raw audio from \"" << in_fname << "\"\n";
}

// search for end of file to count the number of samples to be read
	in_file.seekg(0, std::ios::end);
	// we assume the Python script has written data in 32-bit floats
	const unsigned int num_samples = in_file.tellg();

    std::vector<uint8_t> raw_data(num_samples, 0.0);

	//	for(int i=0; i<raw_data.size(); i++){
	//	std::cout<<"raw_data before read: "<<raw_data[i]<<"\n";
	//}

	// allocate memory space to store all the samples
	//raw_data.clear(); raw_data.resize(num_samples);
	// back to the beginning of the file to read all samples at once
	in_file.seekg(0, std::ios::beg);
	// do a single read for audio data from the input file stream
	in_file.read(reinterpret_cast<char*>(&raw_data[0]), \
						num_samples);
	// close the input file
	in_file.close();
  std::cout << "Read raw RF data from \"" << in_fname << "\" in unsigned 16-bit format\n";

	//for(int i=0; i<raw_data.size(); i++){
	//	std::cout<<"raw_data after read: "<<raw_data[i]<<"\n";
	//}

//for(int b=0; b<raw_data.size(); b++){
		//raw_data[b] = uint16_t(raw_data[b]);
		//std::cout<<"\nindex: "<<b<<"\traw: "<<(raw_data[b]);
	//}

  // IQ data is normalized between -1 and +1 in 32-bit float format
  std::vector<float> iq_data(raw_data.size(), 0.0);
  for (int i = 0; i < raw_data.size(); i++) {
      iq_data[i] = (float(raw_data[i]) - 128.0f) / 128.0f;
  }
/////////////////////
  	for(int i=0; i<iq_data.size(); i++){
		std::cout<<"iq_data: "<<iq_data[i]<<"\n";
	}
	//std::cout<<"\ninput is correct!\n\n";
/////////////////////////


  std::cout << "Reformatted raw RF data to 32-bit float format (" << iq_data.size() * sizeof(float) << " bytes)\n";


/*
	int num_samples = 0;
	std::vector<char> raw_data(num_samples);
	std::vector<float> iq_data(raw_data.size());

	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(int k=0; k<(int)num_samples; k++){
		//normalize
		iq_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}

*/
	///////////////////
	//int num_samples = 0;
	//std::vector<char> raw_data(num_samples,0.0);
	//std::vector<float> iq_data(raw_data.size(),0.0);

	///////////////////////////


/*
	std::ifstream myData("../data/iq_samples.raw", std::ios::binary);
    uint8_t value;
	float valuef = 0.0;
    int i = 0;
    char buf[sizeof(short)];
    while (myData.read(buf, sizeof(buf)))
    {
        memcpy(&value, buf, sizeof(value));
		//valuef = ((float(value)-128.0f)/128.0f);
        //std::cout << valuef << "\n";
        i++;
    }

    std::cout << "\n" << "Total count: " << i << "\n";

	return 0;
*/

	///////////////////////////




	int m = 0;//block_count
	int num_blocks = iq_data.size()/block_size;

	std::vector<float> rf_coeff;
	std::vector<float> audio_coeff;

	std::vector<float> i_filt;
	std::vector<float> q_filt;
	std::vector<float> i_filt_state(rf_taps, 0.0);
	std::vector<float> q_filt_state(rf_taps, 0.0);

	std::vector<float> i_ds(block_size/rf_decim,0.0);
	std::vector<float> q_ds(block_size/rf_decim,0.0);
	float i_state = 1;
	float q_state = 1;

	std::vector<float> fm_demod(block_size/rf_decim,0.0);
	float i_der = 0;
	float q_der = 0;

	std::vector<float> audio_filter(block_size/rf_decim, 0.0);
	std::vector<float> audio_block(block_size/rf_decim/audio_decim, 0.0);//50 lags..? 10 doesnt??
	std::vector<float> processed_data((block_size/rf_decim/audio_decim)*num_blocks, 0.0);

	std::vector<float> yb_i(block_size,0.0);
	std::vector<float> yb_q(block_size,0.0);
	std::vector<float> xb_i(block_size,0.0);
	std::vector<float> xb_q(block_size,0.0);



  // coefficients for the front-end low-pass filter
	impulseResponseLPF(rf_FS,rf_FC,rf_taps,rf_coeff);
	// coefficients for the modes
	impulseResponseLPF(rf_FS,audio_FC,audio_taps,audio_coeff);

	std::vector<float> filter_state(audio_coeff.size()-1, 0);

	while((m+1)*iq_block_size < iq_data.size()){

		//slice into input block
		int bi=0;//index i and q blocks
		int bq=0;
		for(int a=0; a<iq_block_size; a++){//still checking this part, incoming block is 2* size of out i and q blocks
			if(a%2==0){
			xb_i[bi] = iq_data[(m*iq_block_size)+a];//even
			bi++;
			}else{
			xb_q[bq] = iq_data[(m*iq_block_size)+a];//odd
			bq++;
			}//i,q every other value
		}

		for(int b=0; b<xb_i.size(); b++){
			std::cout<<"\nindex: "<<b<<"\txb_i: "<<xb_i[b];
		}

	//	for(int b=0; b<xb_q.size(); b++){
		//	std::cout<<"\nindex: "<<b<<"\txb_q: "<<xb_q[b];
		//}

		//std::cout<<"before conv\n";

		//LPF
		//for (int i =0; i<10000;i++){
			//std::cout<<"Filt at "<<i<<" x: " <<xb_i[i]<< " h: " <<rf_coeff[i]<< " y: "<<yb_i[i]<< " state: "<<i_filt_state[i]<< "\n";
			//std::cout<<"bruh: "<<i_filt_state[i]<<"\n";
		//}
		//convolveFIR(yb_i,xb_i,rf_coeff,i_filt_state,block_size);

		//convolveFIR(yb_q,xb_q,rf_coeff,q_filt_state,block_size);

		//i
		for(int n=0; n<xb_i.size(); n++){
			yb_i[n] = 0;
			for(int k=0; k<rf_coeff.size(); k++){//h size not filter
				if (n<k){
					yb_i[n] += i_filt_state[i_filt_state.size()+n-k]*rf_coeff[k];
					//std::cout<<"for i\ty index: "<<n<<"\tY: "<<yb_i[n];
					//std::cout<<"\tindex: "<<i_filt_state.size()+n-k<<"\tSTATE: "<<i_filt_state[i_filt_state.size()+n-k]<<"\n\n";
				}
				else{
					yb_i[n] += xb_i[n-k]*rf_coeff[k];
					//std::cout<<"y index: "<<n<<"\tY: "<<yb_i[n]<<"\n";
				}
			}
		}
		yb_i.resize(xb_i.size());

		//std::cout<<"DONE \n";

		for(int a=0; a<rf_coeff.size()-1; a++){
			i_filt_state[a] = xb_i[xb_i.size()-rf_coeff.size()+a+1];
			//std::cout<<"for i\tm: "<<m<<"\ty index: "<<yb_i.size()-rf_coeff.size()+a+1<<"\tY: "<<yb_i[yb_i.size()-rf_coeff.size()+a+1]<<"\tstate index: "<<a<<"\tSTATE: "<<i_filt_state[a]<<"\n\n";
		}

		//q
		for(int n=0; n<xb_q.size(); n++){
			yb_q[n] = 0;
			for(int k=0; k<rf_coeff.size(); k++){
				if (n<k){
					yb_q[n] += q_filt_state[q_filt_state.size()+n-k]*rf_coeff[k];
					//std::cout<<"for q\ty index: "<<n<<"\tY: "<<yb_q[n];
					//std::cout<<"\tindex: "<<q_filt_state.size()+n-k<<"\tSTATE: "<<q_filt_state[q_filt_state.size()+n-k]<<"\n\n";
				}
				else{
					yb_q[n] += xb_q[n-k]*rf_coeff[k];
					//std::cout<<"y index: "<<n<<"\tY: "<<yb_i[n]<<"\n";
				}
			}
		}//y resize to x size????????????????????????????????????????
		yb_q.resize(xb_q.size());

		//std::cout<<"DONE \n";

		for(int a=0; a<rf_coeff.size()-1; a++){
			q_filt_state[a] = xb_q[xb_q.size()-rf_coeff.size()+a+1];
			//std::cout<<"for q\tm: "<<m<<"\ty index: "<<yb_q.size()-rf_coeff.size()+a+1<<"\tY: "<<yb_q[yb_q.size()-rf_coeff.size()+a+1]<<"\tstate index: "<<a<<"\tSTATE: "<<q_filt_state[a]<<"\n\n";
		}

		//if (m<3){
			//for (int i =0; i<10000;i++){
				//std::cout<<"Filt at "<<i<<" x: " <<xb_i[i]<< " h: " <<rf_coeff[i]<< " y: "<<yb_i[i]<< " state: "<<i_filt_state[i]<< "\n";
				//std::cout<<"bruh"<<i_filt_state[i]<<"\n";
			//}

		//for(int b=0; b<yb_i.size(); b++){
		//	std::cout<<"\nyb_i: "<<yb_i[b];
		//}

		//downsampling
		for(int i=0; i<yb_i.size(); i++){
			if (i%rf_decim == 0){

				i_ds[i/rf_decim] = yb_i[i];
				q_ds[i/rf_decim] = yb_q[i];
			}
		}

		//for(int b=0; b<i_ds.size(); b++){
		//	std::cout<<"\ni_ds: "<<i_ds[b];
		//}

		//std::cout<<"\nafter ds"<< std::flush;

		//..
		//std::cout<<"\before demod"<< std::flush;

		//demod
		for(int i=0; i<i_ds.size(); i++){//need to increase fm_demod block size because we recombine i and q?

			if(i==0){
				i_der = i_ds[i] - i_state;
				q_der = q_ds[i] - q_state;
			}
			else{
				i_der = i_ds[i] - i_ds[i-1];
				q_der = q_ds[i] - q_ds[i-1];
			}
			if((i_ds[i]!=0)&&(q_ds[i]!=0)){
				fm_demod[i] = (i_ds[i]*q_der - q_ds[i]*i_der) / (pow(i_ds[i],2) + pow(q_ds[i],2));//set to 0 if i==0&q==0 avoid 1/0
			}
			else if (i_ds[i] == 0 && q_ds[i] == 0 && i !=0){
				fm_demod[i] = fm_demod[i-1];
			}
			else{
				fm_demod[i] = 0;
			}
		}
		i_state = i_ds[i_ds.size()-1];//last element
		q_state = q_ds[q_ds.size()-1];//last element

		//i_ds.back() //alt to ^

		//std::cout<<"\nafter demod"<< std::flush;

		//LPF
		//convolveFIR(audio_filter, fm_demod, audio_coeff, filter_state, audio_block_size);

		//for(int b=0; b<fm_demod.size(); b++){
		//	std::cout<<"\nfm_demod: "<<fm_demod[b];
		//}

		//audio
		for(int n=0; n<fm_demod.size(); n++){
			audio_filter[n] = 0;
			for(int k=0; k<audio_coeff.size(); k++){
				if (n<k){
					audio_filter[n] += filter_state[filter_state.size()+n-k]*audio_coeff[k];
					//std::cout<<"for audio\ty index: "<<n<<"\tY: "<<audio_filter[n];
					//std::cout<<"\tindex: "<<filter_state.size()+n-k<<"\tSTATE: "<<filter_state[filter_state.size()+n-k]<<"\n\n";
				}
				else{
					audio_filter[n] += fm_demod[n-k]*audio_coeff[k];
					//std::cout<<"y index: "<<n<<"\tY: "<<yb_i[n]<<"\n";
				}
			}
		}
		audio_filter.resize(fm_demod.size());

		for(int a=0; a<audio_coeff.size()-1; a++){
			filter_state[a] = fm_demod[(int)fm_demod.size()-audio_coeff.size()+a+1];
			//std::cout<<"for audio\tm: "<<m<<"\ty index: "<<audio_filter.size()-audio_coeff.size()+a+1<<"\tY: "<<audio_filter[audio_filter.size()-rf_coeff.size()+a+1]<<"\tstate index: "<<a<<"\tSTATE: "<<filter_state[a]<<"\n\n";
		}

		//for(int b=0; b<audio_filter.size(); b++){
		//	std::cout<<"\naudio_filter: "<<audio_filter[b];
		//}

		//donwsample
		for(int i=0; i<audio_filter.size(); i++){
			if (i%audio_decim == 0){
				audio_block[i/audio_decim] = audio_filter[i];
			}
		}

		//for(int b=0; b<audio_block.size(); b++){
		//	std::cout<<"\naudio_block: "<<audio_block[b];
		//}

		//std::cout<<"\nafter audio ds"<< std::flush;

		//processed_data.insert(processed_data.end(), audio_block.begin(), audio_block.end());//insert should work
		for(int a=0; a<audio_block.size(); a++){
			processed_data[((block_size/rf_decim/audio_decim)*m)+a] = audio_block[a];
		}

		for(int b=0; b<audio_block.size()*(m+1); b++){
			std::cout<<"\processed_data: "<<processed_data[b];
		}

		m++;

		std::cerr<<"\nm: "<<m<<" of: "<<iq_data.size()/iq_block_size;
	}

	std::cerr<<"\ndone blocks\n";
	
	/*
	
	// if you wish to write some binary files, see below example
	//
	 const std::string out_fname = "../data/outdata.bin";
	 writeBinData(out_fname, processed_data);

	std::cout<<"\nwrote to bin file\n";
	//
	// output files can be imported, for example, in Python
	// for additional analysis or alternative forms of visualization

	// naturally, you can comment the line below once you are comfortable to run GNU plot
	//std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	*/
	
	//std::vector<float> processed_data(block_size);
	
	std::vector<short int> audio_data(block_size);
	for(unsigned int k=0; k<processed_data.size(); k++){
		if (std::isnan(processed_data[k])) sample =0;
		
		else sample = static_cast<short int>(processed_data[k] * 16384);
		
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
	}

	return 0;
}
