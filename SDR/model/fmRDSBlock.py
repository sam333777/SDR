

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
from fmRRC import impulseResponseRootRaisedCosine
from fmPll import fmPllRds
from fmRDSFunc import differential_decode, frame_sync
# for take-home add your functions

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

intermediate_Fs = 240e3

audio_Fs = 48e3
audio_decim = 5
audio_taps = 101
audio_Fc = 16e3
# add other settings for audio, like filter taps, ...
offset = 0

def my_firwin(n_taps, fc, fs):
    norm_cutoff = fc / (fs / 2)
    middle_index = (n_taps - 1) / 2
    taps = np.zeros(n_taps)

    for i in range(n_taps):
        if i == middle_index:
            taps[i] = norm_cutoff
        else:
            taps[i] = (
                norm_cutoff
                * np.sin(np.pi * norm_cutoff * (i - middle_index))
                / (np.pi * norm_cutoff * (i - middle_index))
            )
        # Apply Hann window
        taps[i] *= 0.5 - 0.5 * np.cos(2 * np.pi * i / (n_taps - 1))

    return taps

def my_block_process(block, impulse_response, state):
    block_filtered = np.zeros(len(block))
    for i in range(len(block)):
        for j in range(len(impulse_response)):
            if i - j >= 0:
                block_filtered[i] += block[i-j] * impulse_response[j]
            else:
                #use the state from the previous block
                block_filtered[i] += state[i - j + len(state)] * impulse_response[j]
    #update the state
    updated_state = np.concatenate((state[len(block):], block))[-len(state):]
    return block_filtered, updated_state

def delayBlock(input_block, state_block):
	output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
	state_block = input_block[-len(state_block):]

	return output_block, state_block

def fastArctan(I, Q, prev_phase = 0.0):
	#
	# the default prev_phase phase is assumed to be zero, however
	# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))
	
	for k in range(len(I)):
		if k > 0:
			dQ = Q[k] - Q[k - 1] 
			dI = I[k] - I[k - 1]

		else: 
			dQ = 0
			dI = 0

		num = (I[k] * dQ) - (Q[k] * dI)
		den = (I[k]**2) + (Q[k]**2)
		if den != 0:
			temp = num/den
		else:
			temp = 0

		current_phase = prev_phase + temp

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = temp

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

	# for RDS add also the quadrature NCO component to the output
if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples3.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	'''
	# IQ data is normalized between -1 and +1 in 64-bit double format
	iq_data = (np.float64(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 64-bit double format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
	'''

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window=('hann'))

	# coefficients for the filter to extract RDS signal
	# updated according to the coeff and FM signal frequency spectrium based on RDS lecture slides page 5
	rds_signal_coeff = signal.firwin(rf_taps, [54e3 / (intermediate_Fs / 2), 60e3 / (intermediate_Fs / 2)], window=('hann'), pass_zero="bandpass")	# updated according to the coeff and FM signal frequency spectrium based on RDS lecture slides page 6
	# rds_signal_coeff = signal.firwin(rf_taps, [54e3, 60e3], fs=240e3, window=('hann'), pass_zero="bandpass")	# updated according to the coeff and FM signal frequency spectrium based on RDS lecture slides page 6
	rds_pilot_coeff = signal.firwin(audio_taps, [(113.5e3 / (intermediate_Fs / 2)), (114.5e3 / (intermediate_Fs / 2)) ], window=('hann'), pass_zero="bandpass")


	#RSD channel & pilot state initialization
	rds_signal_state = np.zeros(rf_taps - 1)
	rds_pilot_state = np.zeros(rf_taps - 1)
	
	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = (int)(rf_Fs/25) * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps - 1)
	state_q_lpf_100k = np.zeros(rf_taps - 1)
	state_phase = 0

	#RDS Demod coeff
	#RDS PLL
	freq_center = 114e3
	pllRds_state = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]

	#LPF state
	rds_pilot_lpf_coeff = signal.firwin(rf_taps, 3e3/(intermediate_Fs/2), window=('hann'))
	rds_pilot_state_i = np.zeros(rf_taps-1)
	rds_pilot_state_q = np.zeros(rf_taps-1)
	#Modify Up and Down sample rate accordingly, these default values are for the Mode 0 and 2
	mode = 1
	if (mode == 0):
		U = 19
		D = 128
		SPS = 15
	else:
		U = 247
		D = 960
		SPS = 26
	
	
	#Modify Sanples per Symbol accordingly, this default value is for the Mode 0
	
	rrc_Fs = SPS*2375 #from page 26 of the RDS lecture slides
	rrc_taps = 101
	rrc_state_i = np.zeros(rrc_taps - 1)
	rrc_state_q = np.zeros(rrc_taps - 1)

	rrc_filter = impulseResponseRootRaisedCosine(rrc_Fs, rrc_taps)

	#RDS Resampler
	resampler_filter_coeff = signal.firwin(rf_taps, (rrc_Fs/2)/(2.4e5*U/2), window=('hann'))
	resampler_state_i = np.zeros(rf_taps - 1)
	resampler_state_q = np.zeros(rf_taps - 1)
	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file

	# Initialize lists to accumulate recovered_i and recovered_q (No need for implementation)
	all_recovered_i = []
	all_recovered_q = []
	all_decoded_bitstream = []

	rds_signal_state_delayed = np.zeros((int)((audio_taps -1 ) / 2))

	#State saving for Differential Decoding
	symbol_state = -1
	bit_state = 0
	symb_count = 0

	#State saving for Sync
	slot_state = np.zeros(26, dtype=int)
	synced = False
	p_service = ""
	di_prev = -1
	di = -1
	ps_name_segments = [""] * 4
	group_type = -1


	while (block_count + 1) * block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[block_count * block_size:(block_count + 1) * block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

		#RDS signal
		rds_signal, rds_signal_state = signal.lfilter(rds_signal_coeff, 1.0, fm_demod, zi = rds_signal_state)
		#RDS Pilot
		rds_squaring = np.square(rds_signal)*2

		rds_pilot, rds_pilot_state = signal.lfilter(rds_pilot_coeff, 1.0, rds_squaring, zi = rds_pilot_state)
		#plt.plot(rds_pilot, c='r')
		#plt.show()

		#RDS PLL
		rds_pll_inphase, rds_pll_quadrature, pllRds_state = fmPllRds(rds_pilot, freq_center, intermediate_Fs, pllRds_state, 0.5, 0.0, 0.001)
		#All-Pass Filter
		rds_signal_delayed, rds_signal_state_delayed = delayBlock(rds_signal, rds_signal_state_delayed)

		#Mixers
		rds_mixed_i = np.zeros(len(rds_pll_inphase))
		rds_mixed_q = np.zeros(len(rds_pll_quadrature))

		for i in range(len(rds_pll_inphase)-1):
			rds_mixed_i[i] = 2 * rds_signal_delayed[i] * rds_pll_inphase[i]
			rds_mixed_q[i] = 2 * rds_signal_delayed[i] * rds_pll_quadrature[i]


		#Low-Pass Filter
		rds_pilot_i, rds_pilot_state_i = signal.lfilter(rds_pilot_lpf_coeff, 1.0, rds_mixed_i, zi = rds_pilot_state_i)
		rds_pilot_q, rds_pilot_state_q = signal.lfilter(rds_pilot_lpf_coeff, 1.0, rds_mixed_q, zi = rds_pilot_state_q)
		rds_pilot_i_upsampled = np.zeros(len(rds_pilot_i)*U)
		rds_pilot_q_upsampled = np.zeros(len(rds_pilot_q)*U)

		#Rational Resampler
		#Upsampler
		for i in range (len(rds_pilot_i)):
			rds_pilot_i_upsampled[i*U] = rds_pilot_i[i]
			rds_pilot_q_upsampled[i*U] = rds_pilot_q[i]

		#Downsampler
		rds_resample_filt_i,resample_state_i= signal.lfilter(resampler_filter_coeff,1.0, rds_pilot_i_upsampled,zi=resampler_state_i)
		rds_resample_filt_q,resample_state_q= signal.lfilter(resampler_filter_coeff,1.0, rds_pilot_q_upsampled,zi=resampler_state_q)
		
		#The final step of a rational resampler interpolation by factor U, filtering, and then decimation by factor D
		# converts the sample rate of the processed RDS signal from intermediate frequency to 57kHz suitable for the RDS demodulation.
		resampled_i = rds_resample_filt_i[::D]*U
		resampled_q = rds_resample_filt_q[::D]*U
		
		#RRC
		rrc_i, rrc_state_i = signal.lfilter(rrc_filter, 1.0, resampled_i, zi = rrc_state_i)
		rrc_q, rrc_state_q = signal.lfilter(rrc_filter, 1.0, resampled_q, zi = rrc_state_q)
		# to save runtime, select the range of blocks to log data
		# this includes both saving binary files and plotting PSD

		#Data Recovery (Constellation) Slide 43
		#             array[start:stop:step]
		#For debugging only
		# recovered_i = rrc_i[offset-1::SPS]
		# recovered_q = rrc_q[offset-1::SPS]
		#For debugging only
		
		#Determines the optimal starting offset (i.e., best sampling position) 
		# initially by identifying the sample with the highest amplitude within the first symbol period
		# max = np.max(rrc_i[0:SPS])
		# if block_count == 1:
		# 	for i in range (SPS):
		# 		if ((rrc_i[i] > rrc_i[i+1]) and (rrc_i[i]>rrc_i[i-1])) or ((rrc_i[i] < rrc_i[i+1]) and (rrc_i[i] < rrc_i[i-1])):
		# 			offset = i
		# 			break
		if (block_count >= 15) and (block_count % 15 == 0):
			distances = np.zeros(SPS)

			# Calculate the average distance of each sample from the origin
			for i in range(SPS):
				totalDist = 0.0
				numElems = 0

				# Calculate the total distance of each sample from the origin
				for j in range(i, len(rrc_i), SPS):
					totalDist += abs(rrc_i[j])
					numElems += 1
				# Avoid zero division
				if numElems > 0:  
					# Calculate the average distance of each sample from the origin
					distances[i] = totalDist / numElems
				else:
					distances[i] = 0

 		# Gets the index of the max distance
				maxval = 0

			for i in range(len(distances)):
				if distances[i] > maxval:
					maxval = distances[i]
					offset = i % SPS
		print("offset", offset)

		#No need for debugging only
		recovered_i = []
		recovered_q = []
		#No need for debugging only
		symb_count  = offset
		while(symb_count < len(rrc_i)):
			if symb_count > 0:
				recovered_i.append(rrc_i[symb_count])
				recovered_q.append(rrc_q[symb_count])
				for i in range(SPS-1):
					recovered_i.append(0)
					recovered_q.append(0)
			
			symb_count += SPS

		# Accumulate recovered_i and recovered_q (Debugging only)
		all_recovered_i.extend(recovered_i)
		all_recovered_q.extend(recovered_q)
		# Accumulate recovered_i and recovered_q (Debugging only)

		pre_cdr = np.array(rrc_i[offset-1:])

		if not synced and block_count % 15 == 0:
			rrc_i = rrc_i[SPS:]
			rrc_q = rrc_q[SPS:]

		print("Synced", synced)

		decoded_bitstream, symbol_state, bit_state = differential_decode(rrc_i, SPS, symbol_state, bit_state, offset)
		all_decoded_bitstream.extend(decoded_bitstream)
		slot_state, synced, p_service, di_prev, di, ps_name_segments, group_type, ps_ready = frame_sync(decoded_bitstream, slot_state, synced, p_service, di_prev, di, ps_name_segments, group_type)
		if ps_ready:
			print("Program Service:", p_service)
			ps_name_segments = [""] * 4

		

		if block_count == 15:
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			'''
			# create binary file where each sample is a 64-bit double
			fm_demod.astype('float64').tofile(fm_demod_fname)
			'''
			fig, (p_adjust1) = plt.subplots(nrows=1)
			fig, (rrc) = plt.subplots(nrows=1)
			fig.subplots_adjust(hspace = 1.0)
			p_adjust1.scatter(recovered_i, recovered_q, s=10)
			p_adjust1.set_ylim(-1.25, 1.25)
			rrc.plot(rrc_i[0:512], c='r')
			plt.plot(recovered_i[0:512], c='g')
			rrc.plot(rrc_q[0:512], c='b')
			fig, (ax0, ax1, ax2, ax4) = plt.subplots(nrows=4)
			fig.subplots_adjust(hspace = 1.0)

			# PSD after FM demodulation
			ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs / rf_decim) / 1e3)
			ax0.set_ylabel('PSD (db/Hz)')
			ax0.set_title('FM_Demod')

			# save PSD plots
			ax1.psd(rds_pilot, NFFT=512, Fs=(rf_Fs / rf_decim) / 1e3)
			ax1.set_ylabel('PSD (db/Hz)')
			ax1.set_title('RDS_Signal')

			ax2.psd(rds_pll_inphase, NFFT=512, Fs=(rf_Fs / rf_decim) / 1e3)
			ax2.set_ylabel('PSD (db/Hz)')
			ax2.set_title('RDS_After_PLL')

			ax4.psd(rrc_q, NFFT=512, Fs=(rrc_Fs / 1e3))
			ax4.set_ylabel('PSD (db/Hz)')
			ax4.set_title('RRC_Filter')

			# fig, (p_adjust1) = plt.subplots(nrows=1)
			# p_adjust1.set_ylim(-1.25, 1.25)
			# plt.plot(rds_pilot[0:512], c='r')
			
			# plt.plot(rds_pll_inphase[0:512], c='b')
			# fig, (ax0, ax1, ax2, ax4) = plt.subplots(nrows=4)
			# fig.subplots_adjust(hspace = 1.0)

			plt.show()
			# save figure to file
			fig.savefig("../data/fmRDSBlock" + str(block_count) + ".png")


		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmSBlock.wav"
	#wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data / 2) * 32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
    # Plot all recovered_i and recovered_q
	fig, (p_adjust1) = plt.subplots(nrows=1)
	fig.subplots_adjust(hspace=1.0)
	p_adjust1.scatter(all_recovered_i, all_recovered_q, s=10)
	p_adjust1.set_ylim(-1.25, 1.25)
	plt.show()
	# uncomment assuming you wish to show some plots
	# plt.show()
