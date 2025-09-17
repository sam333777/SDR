
#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h" 
#include "logfunc.h"
#include "rdsSupport.h"
#include <cstring>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>

#define Queue_Elems 5

int mode, rf_Fs, rf_Fc, rf_taps, rf_decim, intermediate_Fs, audio_Fs, audio_decim, audio_taps, audio_Fc, audio_up_decim, BLOCK_SIZE, channel;
void RF_Thread(std::queue<std::vector<real>> &audio_queue, std::queue<std::vector<real>> &rds_queue, std::mutex& audio_mutex, std::mutex& rds_mutex, std::condition_variable& audio_var, std::condition_variable& rds_var);
void audio_Thread(std::queue<std::vector<real>> &audio_queue, std::mutex &audio_mutex, std::condition_variable& audio_var);
void rds_Thread(std::queue<std::vector<real>> &rds_queue, std::mutex &rds_mutex, std::condition_variable& rds_var);

int main(int argc, char* argv[])
{

	int mode = 0;

	if (argc < 2)
	{
		std::cerr<< "Operating in default mode 0 and mono mode" << std::endl;
	}
	else if (argc == 2)
	{
		mode = atoi(argv[1]);
		if (mode > 3)
		{
			std::cerr<< "Wrong mode" << mode << std::endl;
			exit(1);
		}
	}
	else if (argc == 3)
	{
		mode = atoi(argv[1]);
		if (mode > 3)
		{
			std::cerr<< "Wrong mode" << mode << std::endl;
			exit(1);
		}

		if (strcmp(argv[2], "m") == 0)
		{
			channel = 0;
		}
		else if (strcmp(argv[2], "s") == 0)
		{
			channel = 1;
		}
		else if (strcmp(argv[2], "r") == 0)
		{
			channel = 2;
		}
		else
		{
			std::cerr<< "Wrong type input, operating in default mono mode" << std::endl;
		}
	}
	else
	{
		std::cerr<< "Usage: " << argv[0] << std::endl;
		std::cerr<< "or: " << std::endl;
		std::cerr<< "Usage: " << argv[0] << "<mode>" << std::endl;
		std::cerr<< "\t\t <mode> is a value from 0 to 3" << argv[0] << std::endl;
		exit(1);
	}

	switch(mode){
		case 0:
			rf_Fs = 2.4e6;
			rf_Fc = 100e3;
			rf_taps = 101;
			rf_decim = 10;

			intermediate_Fs = 240e3;

			audio_Fs = 48e3;
			audio_decim = 5;
			audio_up_decim = 1;
			audio_taps = 101;
			audio_Fc = 16e3;
			break;

		case 1:
			rf_Fs = 1.92e6;
			rf_Fc = 100e3;
			rf_taps = 101;
			rf_decim = 12;

			intermediate_Fs = 160e3;

			audio_Fs = 40e3;
			audio_decim = 4;
			audio_up_decim = 1;
			audio_taps = 101;
			audio_Fc = 16e3;
			if (channel == 2)
			{
				channel = 1;
			}
			break;

		case 2:
			rf_Fs = 2.4e6;
			rf_Fc = 100e3;
			rf_taps = 101;
			rf_decim = 10;
		
			intermediate_Fs = 240e3;
		
			audio_Fs = 44.1e3;
			audio_up_decim = 147;
			audio_decim = 800;
			audio_taps = 101;
			audio_Fc = 16e3;
			break;
		case 3:
			rf_Fs = 1.44e6;
			rf_Fc = 100e3;
			rf_taps = 101;
			rf_decim = 5;
		
			intermediate_Fs = 288e3;
		
			audio_Fs = 44.1e3;
			audio_decim = 320;
			audio_up_decim = 49;
			audio_taps = 101;
			audio_Fc = 16e3;
			if (channel == 2)
			{
				channel = 1;
			}
			break;

		default:
			rf_Fs = 2.4e6;
			rf_Fc = 100e3;
			rf_taps = 101;
			rf_decim = 10;

			intermediate_Fs = 240e3;

			audio_Fs = 48e3;
			audio_decim = 5;
			audio_taps = 101;
			audio_Fc = 16e3;
			break;
	}

	BLOCK_SIZE = (rf_Fs / 25) * 2;

	std::queue<std::vector<real>> audio_queue;
	std::queue<std::vector<real>> rds_queue;
	std::mutex audio_mutex;
	std::mutex rds_mutex;
	std::condition_variable audio_var;
	std::condition_variable rds_var;

	std::thread rf_thread = std::thread(RF_Thread, std::ref(audio_queue), std::ref(rds_queue), std::ref(audio_mutex), std::ref(rds_mutex), std::ref(audio_var), std::ref(rds_var));
	std::thread audio_thread = std::thread(audio_Thread, std::ref(audio_queue), std::ref(audio_mutex), std::ref(audio_var));
	//std::thread rds_thread = std::thread(rds_Thread, std::ref(rds_queue), std::ref(rds_mutex), std::ref(rds_var));

	rf_thread.join();
	audio_thread.join();
	//rds_thread.join();
	return 0;
}

void RF_Thread(std::queue<std::vector<real>> &audio_queue, std::queue<std::vector<real>> &rds_queue, std::mutex& audio_mutex, std::mutex& rds_mutex, std::condition_variable& audio_var, std::condition_variable& rds_var)
{
	real phase_diff = 0;

	std::vector<real> stateILpf100K;
	std::vector<real> stateQLpf100K;

	stateILpf100K.resize(rf_taps, 0.0);
	stateQLpf100K.resize(rf_taps, 0.0);

	std::vector<real> frontEndFilter;
	frontEndFilter.resize(rf_taps,0.0);

	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, frontEndFilter);

	for (unsigned int block_id = 0; ; block_id++)
	{
		std::vector<real> block_data(BLOCK_SIZE);
		readStdinBlockData(BLOCK_SIZE, block_id, block_data);
		//std::cerr << "Read block" << block_id << std::endl;

		std::vector<real> I_Data(BLOCK_SIZE/2);
		std::vector<real> Q_Data(BLOCK_SIZE/2);

		int k = 0;

		for (int i = 0; i < BLOCK_SIZE; i = i + 2)
		{
			I_Data[k] = block_data[i];
			Q_Data[k] = block_data[i+1];
			k += 1;
		}
		
		std::vector<real> I_Filt_Data(BLOCK_SIZE/(2*rf_decim));
		std::vector<real> Q_Filt_Data(BLOCK_SIZE/(2*rf_decim));

		downSampler(I_Filt_Data, I_Data, frontEndFilter, rf_decim, stateILpf100K);
		downSampler(Q_Filt_Data, Q_Data, frontEndFilter, rf_decim, stateQLpf100K);

		std::vector<real> FM_Demod;
		FM_Demod.resize(BLOCK_SIZE/(2*rf_decim),0.0);
		lightweight_fmDemodArctan(FM_Demod, phase_diff, I_Filt_Data, Q_Filt_Data);
		
		std::unique_lock<std::mutex> audio_lock(audio_mutex);
		while (audio_queue.size() >= Queue_Elems)
		{
			audio_var.wait(audio_lock);
		}
		audio_queue.push(FM_Demod);
		audio_var.notify_one();
		audio_lock.unlock();

		std::unique_lock<std::mutex> rds_lock(rds_mutex);
		while (rds_queue.size() >= Queue_Elems)
		{
			rds_var.wait(rds_lock);
		}
		rds_queue.push(FM_Demod);
		rds_var.notify_one();
		rds_lock.unlock();
	}
	std::cerr << "RF Thread" << std::endl;
}

void audio_Thread(std::queue<std::vector<real>> &audio_queue, std::mutex &audio_mutex, std::condition_variable& audio_var)
{
	std::vector<real> audioFiltState(rf_taps - 1, 0.0);;
	std::vector<real> stereo_pilot_state(audio_taps - 1, 0.0);
	std::vector<real> stereo_audio_state(audio_taps - 1, 0.0);
	std::vector<real> stereo_filt_state(audio_taps - 1, 0.0);
	std::vector<real> audio_filt_state(audio_taps - 1, 0.0);
	std::vector<real> stereo_pilot_coeff(audio_taps, 0.0);
	std::vector<real> stereo_audio_coeff(audio_taps, 0.0);
	std::vector<real> audio_coeff(audio_taps, 0.0);
	real prevI = 1.0;
	real prevQ = 0.0;
	real integrator = 0.0;
	real eastmatePhase = 0.0;
	real trigOffset = 0;
	real firstNco = 0.0;

	std::vector<real> ncoOut;
	std::vector<real> stereo_pilot_filt(BLOCK_SIZE* audio_up_decim /(rf_decim*2), 0.0);
	std::vector<real> stereo_audio_filt(BLOCK_SIZE* audio_up_decim /(rf_decim*2), 0.0);
	

	std::vector<real> audio_filt_state_delayed((audio_taps - 1) / 2, 0.0); 
	std::vector<real> audio_filt_delayed(BLOCK_SIZE * audio_up_decim / (2 * rf_decim), 0.0);
	std::vector<real> stereo_audio_filt_final(BLOCK_SIZE * audio_up_decim / (2 * rf_decim), 0.0); 
	std::vector<real> audio_filt(BLOCK_SIZE * audio_up_decim / (2 * rf_decim), 0.0); 
	std::vector<real> mono_audio_block(BLOCK_SIZE * audio_up_decim / (2 * rf_decim * audio_decim), 0.0); 
	std::vector<real> stereo_audio_block(BLOCK_SIZE * audio_up_decim / (2 * rf_decim * audio_decim), 0.0); 

    std::vector<real> intermediateFilter(audio_taps * audio_up_decim - 1, 0.0);	
	std::vector<real> Stereo_Resampling_Coeff;

	impulseResponseBPF(intermediate_Fs, 18.5e3, 19.5e3, audio_taps, stereo_pilot_coeff);
	impulseResponseBPF(intermediate_Fs, 22e3, 54e3, audio_taps, stereo_audio_coeff);
	impulseResponseLPFwithResamplerFactor(intermediate_Fs, audio_Fc, audio_taps, Stereo_Resampling_Coeff, audio_up_decim);
	impulseResponseLPFwithResamplerFactor(intermediate_Fs, audio_Fc, audio_taps, intermediateFilter, audio_up_decim);


	for (unsigned int block_id = 0; ; block_id++)
	{
		
		std::vector<real> FM_Demod;
		FM_Demod.resize(BLOCK_SIZE/(2*rf_decim),0.0);

		// while (audio_queue.size() < 1) 
		// {
		// 	std::cerr << "Waiting for data" << std::endl;
		// }
		std::unique_lock<std::mutex> audio_lock(audio_mutex);
		//while (audio_queue.size() >= Queue_Elems)
		while (audio_queue.empty())
		{
			audio_var.wait(audio_lock);
		}
		FM_Demod = audio_queue.front();
		audio_queue.pop();
		audio_var.notify_one();
		audio_lock.unlock();
		
		if (channel == 1 || channel == 2) { // Stereo branch

			
			downSampler(stereo_pilot_filt, FM_Demod, stereo_pilot_coeff, 1, stereo_pilot_state);
			downSampler(stereo_audio_filt, FM_Demod, stereo_audio_coeff, 1, stereo_audio_state);

			allPassFilter(audio_filt_delayed, FM_Demod, audio_filt_state_delayed);
			resampler_optimized(FM_Demod, audio_filt_delayed, intermediateFilter, audio_decim, audio_up_decim, audioFiltState);
			
			// PLL on the pilot-filtered signal
			fmPllwithState(stereo_pilot_filt, 19e3, intermediate_Fs, 2.0, 0.0, 0.01, prevI, prevQ, integrator, eastmatePhase, trigOffset, ncoOut, firstNco);
			// Mixer
			for (size_t i = 0; i < ncoOut.size(); i++) {
				stereo_audio_filt[i] = 2.0 * ncoOut[i] * stereo_audio_filt[i];
			}

			//filter stero
			resampler_optimized(stereo_audio_filt_final, stereo_audio_filt, intermediateFilter, audio_decim, audio_up_decim, stereo_filt_state);
			// filter mono
			resampler_optimized(audio_filt, audio_filt_delayed, intermediateFilter, audio_decim, audio_up_decim, audioFiltState);
		
			resampler_optimized(stereo_audio_filt_final, stereo_audio_filt, Stereo_Resampling_Coeff, audio_decim, audio_up_decim, stereo_filt_state);

			// separate channels
			//size_t numOut = (mono_audio_block.size() > 0) ? mono_audio_block.size() : 0;
			std::vector<real> stereo_left(mono_audio_block.size(), 0.0);
			std::vector<real> stereo_right(mono_audio_block.size(), 0.0);
			for (int i = 0; i < mono_audio_block.size(); i++) {
				
				stereo_left[i]  = (audio_filt[i] + stereo_audio_filt_final[i]);
				stereo_right[i] = (stereo_audio_filt_final[i] - audio_filt[i]);
				//std::cerr << "L: "<< stereo_left[i]<< std::endl;
				//std::cerr << "R: "<< stereo_right[i]<< std::endl;
			}
			
			std::vector<real> audio_data;
			// // concatinatino 
			for (int i = 0; i < stereo_left.size(); i++) {
				
				audio_data.push_back(stereo_left[i]);
				audio_data.push_back(stereo_right[i]);
			}
			writeStoutBlockData(audio_data.size(), audio_data);
			//std::cerr << "Wrote block" << audio_data.size() << std::endl;
		} 
		else //mono
		{
			std::vector<real> final_audio(BLOCK_SIZE*audio_up_decim/(rf_decim*audio_decim*2));
			resampler_optimized(final_audio, FM_Demod, intermediateFilter, audio_decim, audio_up_decim, audioFiltState);

			
			// for (auto element: final_audio){   //debug loop to see contents of a vector
			// 	std::cerr << element << " ";
			// }
			
			writeStoutBlockData(final_audio.size(), final_audio);
		}

		if ((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			break;
		}
	}
	//std::cerr << "audio Thread" << std::endl;
}

/*
void rds_Thread(std::queue<std::vector<real>> &rds_queue, std::mutex &rds_mutex, std::condition_variable& rds_var)
{	
	
	std::vector<real> rds_signal_coeff(audio_taps, 0.0);
	std::vector<real> rds_signal_state(audio_taps - 1, 0.0);
	std::vector<real> rds_signal_filt(BLOCK_SIZE* audio_up_decim /(rf_decim*2), 0.0);

	std::vector<real> rds_signal_delayed_state(audio_taps/2 - 1, 0.0);
	std::vector<real> rds_signal_delayed_filt(BLOCK_SIZE* audio_up_decim /(rf_decim*2), 0.0);

	std::vector<real> rds_squaring(BLOCK_SIZE* audio_up_decim /(rf_decim*2), 0.0);

	std::vector<real> rds_pilot_coeff(audio_taps, 0.0);
	std::vector<real> rds_pilot_state(audio_taps - 1, 0.0);
	std::vector<real> rds_pilot_filt(BLOCK_SIZE* audio_up_decim /(rf_decim*2), 0.0);
	
	std::vector<real> rrc_i;
	std::vector<real> rrc_i_state(audio_taps-1, 0.0);
	std::vector<real> all_recovered_i;

	real prevI = 1.0;
	real prevQ = 0.0;
	real integrator = 0.0;
	real eastmatePhase = 0.0;
	real trigOffset = 0;
	real firstNco = 0.0;

	std::vector<real> ncoOut;

	int decodeIdentifier;
	int decodeIdentifier_prev;
	bool synced;
	bool ps_ready = false;
	std::string p_service;
	std::vector<bool> slot_state;
	int bit_state;
	int symbol_state;
	int block_count;

	std::vector<bool> decoded_bitstream;
	std::vector<bool> all_decoded_bitstream;

	std::vector<std::string> ps_name_segments(4,"");
	std::string group_type;


	impulseResponseBPF(intermediate_Fs, 54e3, 60e3, audio_taps, rds_signal_coeff);
	impulseResponseBPF(intermediate_Fs, 113.5e3, 114.5e3, audio_taps, rds_pilot_coeff);

	
	int U, D, SPS;

	if (mode == 0)
	{
		U = 19;
		D = 128;
		SPS = 15;
	}
	else
	{
		U = 247;
		D = 960;
		SPS = 26;
	}

	std::vector<real> rds_mixed_i(BLOCK_SIZE* U /(D*rf_decim*2), 0.0);
	std::vector<real> rds_mixed_i_state(audio_taps - 1, 0.0);

	std::vector<real> InPhaseFilter(audio_taps * U - 1, 0.0);	

	impulseResponseLPFwithResamplerFactor(intermediate_Fs, 3e3, audio_taps, InPhaseFilter, U);

	int rrc_Fs = SPS*2375;
	int rrc_taps = 101;
	int offset;

	std::vector<real> I_resampler_state(rrc_taps - 1, 0.0);
	std::vector<real> Q_resampler_state(rrc_taps - 1, 0.0);

	std::vector<real> rrc_coeff(rrc_taps, 0.0);

	impulseResponseRootRaisedCosine(rrc_coeff, rrc_Fs, rrc_taps);

	for (unsigned int block_id = 0; ; block_id++)
	{
		
		std::vector<real> FM_Demod;
		FM_Demod.resize(BLOCK_SIZE/(2*rf_decim),0.0);

		std::unique_lock<std::mutex> rds_lock(rds_mutex);
		while (rds_queue.empty())
		{
			rds_var.wait(rds_lock);
		}
		FM_Demod = rds_queue.front();
		rds_queue.pop();
		rds_var.notify_one();
		rds_lock.unlock();
		
		if (channel == 2)
		{
			downSampler(rds_signal_filt, FM_Demod, rds_signal_coeff, 1, rds_signal_state); //isolates the rds signal from 54 to 60 kHz
			for (int i = 0; i < rds_signal_filt.size(); i++){
				rds_pilot_filt[i] = rds_signal_filt[i]*rds_signal_filt[i] * 2;  //performs squaring non-linearlity
			}
			downSampler(rds_pilot_filt, rds_pilot_filt, rds_pilot_coeff, 1, rds_pilot_state); //isolates the pilot tone by filtering aroung 114kHz

			// PLL on the pilot-filtered signal
			//std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#####################################" << std::endl;
			fmPllwithState(rds_pilot_filt, 114e3, intermediate_Fs, 0.5, 0.0, 0.001, prevI, prevQ, integrator, eastmatePhase, trigOffset, ncoOut, firstNco); //locks the pilot tone using a pll, no phase shift needed

			allPassFilter(rds_signal_delayed_filt, rds_signal_filt, rds_signal_delayed_state); // all pass on the orignal filtered signal to align the phase with the pilot tone
			
			// Mixers
			for (size_t i = 0; i < ncoOut.size(); i++) {
				rds_mixed_i[i] = 2.0 * ncoOut[i] * rds_signal_delayed_filt[i];  //performs a mixing of the signal and the pilot
			}

			resampler_optimized(rds_mixed_i, rds_mixed_i, InPhaseFilter, D, U, rds_mixed_i_state); //resamples the signal for symbol extraction

			downSampler(rrc_i, rds_mixed_i, rrc_coeff, 1, rrc_i_state); //applies the root raised filter

			if(block_count >= 10 && block_count % 10 == 0){ //starts rds recovery at block 10
				std::vector<real> distances(SPS, 0.0);
				int maxVal;
				for (int i = 0; i < SPS; i++){ //loops through all the elements in the range of SPS to find the max value
					real totalDist = 0;
					real numElems = 0;

					for (int j = 0; j < rrc_i.size()-1; j+=SPS){ //grabs every value in this range
						totalDist += std::abs(rrc_i[j]);
						numElems += 1;
					}

					if(numElems > 0){ //if there are elements in the array, finds the average value
						distances[i] = totalDist/numElems;
					} else {
						distances[i] = 0;
					}

					maxVal = 0;
				}

				for (int i = 0; i < distances.size()-1; i++){ //compares all of the average values for each sps range in the finds the highest
					if (distances[i] > maxVal){
						maxVal = distances[i];
						offset = i % SPS; //modulates by sps to align with the first wave
					}
				}
			}

			std::vector<real> recovered_i;
			int symb_count = offset;

			while(symb_count < rrc_i.size()-1){ //for debugging plots
				if (symb_count > 0){
					recovered_i.push_back(rrc_i[symb_count]);
					for (int i = 0; i < SPS-1; i++) {	
						recovered_i.push_back(0);
					}
				}
				symb_count += SPS;
			}

			all_recovered_i.insert(all_recovered_i.end(), recovered_i.begin(), recovered_i.end()); //for debugging plots

			std::vector<real> pre_cdr;

			pre_cdr.assign(rrc_i.begin() + offset-1, rrc_i.end()); //starts the recovery in the wave at the offset

			if (!synced && block_count % 15 == 0) //if the clock is not aligned, realigns the signal
			{
				rrc_i.assign(rrc_i.begin() + SPS, rrc_i.end());
			} 

			differential_decode(rrc_i, SPS, symbol_state, bit_state, symb_count); //performs the manchester decoding and differential decoding
			//all_decoded_bitstream.insert(all_decoded_bitstream.end(), decoded_bitstream.begin(), decoded_bitstream.end());
			frame_sync(decoded_bitstream, slot_state, synced, p_service, decodeIdentifier_prev, decodeIdentifier, ps_name_segments, group_type, ps_ready); //performs the frame syncronization
			if (ps_ready){ //if the ps has been recovered, will print out what has been recovered
				std::cerr << "Program Services:" << p_service;
				ps_name_segments.clear(); //clears the screen for the next time its ready
			}
		}
	}
}
*/